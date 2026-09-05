"""
    Cumulative

Accumulators shared by every pair solve of a pairwise run, guarded by `lock`:
the cumulative and maximum current grids of a raster (`cum_curr`,
`max_curr`; `max_curr` is empty unless `write_max_cur_maps`), or the
cumulative node and branch currents of a network (`cum_node_curr`,
`cum_branch_curr`, with `coords` naming each branch and `branch_index`
mapping an edge in either orientation to its slot).
"""
struct Cumulative{T,V}
    cum_curr::Matrix{T}
    max_curr::Matrix{T}
    cum_branch_curr::Vector{T}
    cum_node_curr::Vector{T}
    coords::Vector{Tuple{V,V}}
    # Edge (in either orientation) -> position in coords / cum_branch_curr,
    # so accumulating a pair's branch currents is O(E), not O(E^2).
    branch_index::Dict{Tuple{V,V},Int}
    lock::ReentrantLock
end

"""
    GraphProblem

A pairwise problem: the graph Laplacian `G`, its connected components `cc`,
the focal nodes as graph indices (`points`) with the user's ids alongside
(`user_points`), the pairs of user ids not to solve (`exclude_pairs`), the
[`Geometry`](@ref) that maps nodes back to output coordinates, the
[`Cumulative`](@ref) accumulators and the solver. Raster and network runs
differ only in `geometry`. Built by [`build_problem`](@ref) and consumed by
[`solve`](@ref).
"""
struct GraphProblem{T,V,W,Geom<:Geometry}
    G::SparseMatrixCSC{T,V}
    cc::Vector{Vector{V}}
    points::Vector{V}
    user_points::Vector{V}
    exclude_pairs::Vector{Tuple{V,V}}
    geometry::Geom
    cum::Cumulative{T}
    solver::W
end

"""
    ComponentData

One connected component as the output writers see it: its node list, its
matrix and its geometry restricted to the component's local numbering.
"""
struct ComponentData{T,V,Geom<:Geometry}
    cc::Vector{V}
    matrix::SparseMatrixCSC{T,V}
    geometry::Geom
end

struct Output{T,V}
    points::Vector{V}
    voltages::Vector{T}
    orig_pts::Tuple{V,V}
    comp_idx::Tuple{V,V}
    resistance::T
    col::V
    cum::Cumulative{T}
end

# The pairwise kernel hands `write_cur_maps` / `write_volt_maps` an `Output`
# per pair; the advanced kernel hands them the voltage vector itself.
_voltages(output::Output) = output.voltages
_voltages(voltages::AbstractVector) = voltages

struct Shortcut{T}
    get_shortcut_resistances::Bool
    voltmatrix::Matrix{T}
    shortcut_res::Matrix{T}
end

"""
    Solver

Abstract supertype of the solver backends. `prepare!` and `solve_pairs!`
dispatch on it; everything else in the pairwise driver is shared.
"""
abstract type Solver end

"""
    CholmodSolver(bs)

Direct solver: sparse Cholesky factorization through SuiteSparse CHOLMOD
(`solver = cholmod`), solving `bs` right-hand sides per batch. Available
without extra packages, in double and single precision.
"""
struct CholmodSolver <: Solver
    bs::Int
end

"""
    AMGSolver()

Iterative solver: conjugate gradient preconditioned by smoothed-aggregation
algebraic multigrid (`solver = cg+amg`, the default). One preconditioner is
built per connected component and shared, with a private workspace per task,
by the parallel pair solves.
"""
struct AMGSolver <: Solver
end

"""
    PardisoSolver(bs)

Direct solver through Intel MKL Pardiso (`solver = pardiso`), `bs` right-hand
sides per batch. Requires the Pardiso.jl package extension (`using Pardiso`)
and double precision.
"""
struct PardisoSolver <: Solver
    bs::Int
end

"""
    AccelerateSolver(bs)

Direct solver through Apple's Accelerate sparse Cholesky
(`solver = accelerate`), `bs` right-hand sides per batch. Requires the
AppleAccelerate.jl package extension (`using AppleAccelerate`) on macOS;
double and single precision.
"""
struct AccelerateSolver <: Solver
    bs::Int
end
"""
Core kernel of Circuitscape - used to solve several pairs

Input:
* data::GraphProblem
"""
function single_ground_all_pairs(prob::GraphProblem{T,V,W}, cfg, log = true) where {T,V,W}
    solve(prob, prob.solver, cfg, log)
end

"""
    get_solver(cfg) -> Solver

The solver backend selected by `cfg.solver`, with `cfg.cholmod_batch_size`
as the batch size of the direct solvers. Logs which solver is used.
"""
function get_solver(cfg)
    s = cfg.solver
    if s == st_cg_amg
        @info("Solver used: AMG accelerated by CG")
        return AMGSolver()
    elseif s == st_cholmod
        @info("Solver used: CHOLMOD")
        bs = cfg.cholmod_batch_size
        return CholmodSolver(bs)
    elseif s == st_pardiso
        @info("Solver used: Pardiso")
        bs = cfg.cholmod_batch_size
        return PardisoSolver(bs)
    elseif s == st_accelerate
        @info("Solver used: Apple Accelerate")
        bs = cfg.cholmod_batch_size
        return AccelerateSolver(bs)
    else
        error("Unknown solver: $s")
    end
end

const DirectSolver = Union{CholmodSolver, PardisoSolver, AccelerateSolver}

"""
    PairJob

One pair of focal nodes to solve within a connected component: the two graph
nodes (`src_node`, `dst_node`), their local indices in the component
(`comp_i`, `comp_j`), and `src_indices` / `dst_indices`, the positions in
`points` that map to each node; a focal region contributes several positions,
all sharing one solve. Produced by [`pair_jobs`](@ref).
"""
struct PairJob{V}
    src_node::V
    dst_node::V
    comp_i::V
    comp_j::V
    src_indices::Vector{Int}
    dst_indices::Vector{Int}
end

# Components from `connected_components` list their nodes in ascending order,
# so a node's local index and its membership are binary searches: O(log n)
# each and allocation-free, where a linear scan or a per-component hash table
# would be O(n) on a component of millions of nodes.
function component_index(comp, node)
    i = searchsortedfirst(comp, node)
    (i <= length(comp) && comp[i] == node) ? i : nothing
end

in_component(comp, node) = component_index(comp, node) !== nothing

# Distinct focal nodes of one component, in first-appearance order of `points`.
component_points(comp, points) = unique(filter(x -> in_component(comp, x), points))

# Positions in `points` that map to each node; a focal region contributes
# several positions to one node.
function point_positions(points::Vector{V}) where V
    positions = Dict{V,Vector{Int}}()
    for (i, p) in enumerate(points)
        push!(get!(() -> Int[], positions, p), i)
    end
    positions
end

"""
    pair_jobs(csub, comp, points, orig_pts, exclude, first_source_only;
              positions = point_positions(points))

The node pairs to solve in one component, grouped by source node, in the
order the log numbers them. A pair is skipped when every combination of its
focal indices is excluded. With `first_source_only` (the resistance shortcut)
only pairs anchored at `csub[1]` are produced. A group is returned for every
source considered, even when it has no pairs, so callers can iterate sources.
`comp` must be sorted ascending, as `connected_components` returns it.
"""
function pair_jobs(csub::Vector{V}, comp, points, orig_pts, exclude, first_source_only;
                   positions = point_positions(points)) where V
    issorted(comp) || throw(ArgumentError("component node list must be sorted ascending"))
    excluded = Set(exclude)
    groups = Vector{Vector{PairJob{V}}}()
    nsrc = first_source_only ? min(1, length(csub)) : length(csub)
    for si in 1:nsrc
        src_node = csub[si]
        comp_i = component_index(comp, src_node)
        comp_i === nothing && error("Node $src_node not found in component")
        src_indices = positions[src_node]
        jobs = PairJob{V}[]
        for di in si+1:length(csub)
            dst_node = csub[di]
            dst_node == src_node && continue
            comp_j = component_index(comp, dst_node)
            comp_j === nothing && error("Node $dst_node not found in component")
            dst_indices = positions[dst_node]
            any(((orig_pts[c_i], orig_pts[c_j]) ∉ excluded
                 for c_i in src_indices, c_j in dst_indices)) || continue
            push!(jobs, PairJob{V}(src_node, dst_node, V(comp_i), V(comp_j),
                                   src_indices, dst_indices))
        end
        push!(groups, jobs)
    end
    groups
end

"""
    solve(prob, solver, cfg, log)

Pairwise driver shared by every solver. Enumerates the pairs per connected
component, hands them to `solve_pairs!` for the solver-specific linear
algebra, and stores resistances / postprocesses maps in `handle_pair`. Only
`prepare!` and `solve_pairs!` differ between the iterative and direct paths.
"""
function solve(prob::GraphProblem{T,V}, solver::Solver, cfg, log)::Matrix{T} where {T,V}

    # Data
    a = prob.G
    cc = prob.cc
    points = prob.points
    exclude = prob.exclude_pairs
    orig_pts = prob.user_points
    geometry = prob.geometry
    cum = prob.cum

    # Output options
    write_volt_maps = cfg.write_volt_maps
    write_cur_maps = cfg.write_cur_maps
    write_cum_cur_map_only = cfg.write_cum_cur_map_only
    write_max_cur_maps = cfg.write_max_cur_maps

    numpoints = size(points, 1)

    @info("Graph has $(size(a,1)) nodes, $numpoints focal points and $(length(cc)) connected components")

    num_pairs, pair_numbers = get_num_pairs(cc, points, exclude, orig_pts)
    log && @info("Total number of pair solves = $num_pairs")

    # Initialize pairwise resistance
    resistances = -1 * ones(T, numpoints, numpoints)
    voltmatrix = zeros(T, size(resistances))
    shortcut_res = -1 * ones(T, size(resistances))

    # Get a vector of connected components
    comps = getindex.([a], cc, cc)

    get_shortcut_resistances = false
    if is_raster(geometry) && !write_volt_maps && !write_cur_maps &&
            !write_cum_cur_map_only && !write_max_cur_maps &&
            isempty(exclude)
        get_shortcut_resistances = true
        log && @info("Triggering resistance calculation shortcut")
        num_pairs, pair_numbers = get_num_pairs_shortcut(cc, points, exclude, orig_pts)
        log && @info("Total number of pair solves has been reduced to $num_pairs")
    end
    shortcut = Shortcut(get_shortcut_resistances, voltmatrix, shortcut_res)

    pair_label(job) = haskey(pair_numbers, (job.src_node, job.dst_node)) ?
        "Solving pair $(pair_numbers[(job.src_node, job.dst_node)]) of $num_pairs" :
        "Solving pair $(job.src_node)-$(job.dst_node)"

    positions = point_positions(points)

    for (cid, comp) in enumerate(cc)

        # Subset of points relevant to CC
        csub = component_points(comp, points)
        isempty(csub) && continue

        # Conductance matrix corresponding to CC; the solver may regularize
        # it in place before it is captured for output writing.
        matrix = comps[cid]
        handle = prepare!(solver, matrix)

        # Geometry of this CC in local numbering - for output writing
        local_geometry = @timeit CSTIMER[] "construct local nodemap" restrict(geometry, comp)

        component_data = ComponentData(comp, matrix, local_geometry)

        groups = pair_jobs(csub, comp, points, orig_pts, exclude, get_shortcut_resistances;
                           positions)

        # Focal indices that share a node have zero resistance between them
        for si in eachindex(groups)
            smash_repeats!(resistances, positions[csub[si]])
        end

        # Store one solve's result for every non-excluded index combination
        # and postprocess its maps. Runs on worker threads; each (c_i, c_j)
        # is owned by exactly one job, so the writes never overlap.
        function handle_pair(job, voltages)
            resistance = voltages[job.comp_j] - voltages[job.comp_i]
            for c_i in job.src_indices, c_j in job.dst_indices
                (orig_pts[c_i], orig_pts[c_j]) in exclude && continue
                resistances[c_i, c_j] = resistance
                resistances[c_j, c_i] = resistance
                output = Output(points, voltages, (orig_pts[c_i], orig_pts[c_j]),
                                (job.comp_i, job.comp_j), resistance, V(c_j), cum)
                postprocess(output, component_data, shortcut, cfg)
            end
        end

        solve_pairs!(handle_pair, handle, solver, matrix, groups, cfg, log, pair_label)

        if get_shortcut_resistances
            idx = first(positions[csub[1]])
            update_shortcut_resistances!(idx, shortcut, resistances, points, comp)
        end
    end

    if get_shortcut_resistances
        resistances = shortcut.shortcut_res
    end

    for i = 1:size(resistances,1)
        resistances[i,i] = 0
    end

    # Pad it with the user points
    r = vcat(vcat(0,orig_pts)', hcat(orig_pts, resistances))

    # Save resistances
    save_resistances(r, cfg)

    r
end

"""
    prepare!(solver, matrix)

Per-component setup, run once before the pairs of a component are solved.
Returns whatever `solve_pairs!` needs for that solver: the AMG preconditioner
for [`AMGSolver`](@ref) (after regularizing `matrix` in place), or the
factorization for a direct solver.
"""
function prepare!(::AMGSolver, matrix)
    # Regularization step
    matrix.nzval .+= eps(eltype(matrix)) * norm(matrix.nzval)

    # Construct preconditioner *once* for every CC.
    # CG requires a symmetric positive definite preconditioner, so the coarse
    # solve must be a symmetric operator. AlgebraicMultigrid's default coarse
    # solver (QR) is not, and its factorization is also not safe to share
    # across the threaded pair solves below; the pseudo-inverse is both.
    @timeit CSTIMER[] "construct preconditioner" aspreconditioner(
        smoothed_aggregation(matrix;
                             coarse_solver = AlgebraicMultigrid.Pinv,
                             presmoother = AlgebraicMultigrid.GaussSeidel(),
                             postsmoother = AlgebraicMultigrid.GaussSeidel()))
end

prepare!(solver::DirectSolver, matrix) =
    @timeit CSTIMER[] "construct cholesky factor" construct_cholesky_factor(matrix, solver)

"""
    solve_pairs!(handle_pair, handle, solver, matrix, groups, cfg, log, pair_label)

Solve every job in `groups` and call `handle_pair(job, voltages)` with the
voltages referenced to the source node. The iterative path solves one pair at
a time in equal-sized chunks of pairs, each task with its own preconditioner
workspace; the direct path factorizes once and solves batches
of right-hand sides, parallelizing the postprocessing of each batch.
"""
function solve_pairs!(handle_pair, P, ::AMGSolver, matrix::SparseMatrixCSC{T},
                      groups, cfg, log, pair_label) where T

    function task(jobs, timer)
        @timeit timer "task" begin
        # Each task needs its own workspace (scratch vectors are mutable)
        ml = P.ml
        local_P = aspreconditioner(AlgebraicMultigrid.MultiLevel(
            ml.levels, ml.final_A, ml.coarse_solver,
            ml.presmoother, ml.postsmoother, deepcopy(ml.workspace)))

        for job in jobs
            current = zeros(T, size(matrix, 1))
            current[job.comp_i] = -1
            current[job.comp_j] = 1

            log && @debug(pair_label(job))
            voltages = @timeit timer "solve linear system" solve_linear_system(matrix, current, local_P; tol = residual_tolerance(cfg, T))
            voltages .-= voltages[job.comp_i]
            @timeit timer "postprocess" handle_pair(job, voltages)
        end
        end # @timeit task
        timer
    end

    # Every pair costs about the same (same matrix, same preconditioner), so
    # equal-count chunks are equal-work chunks. Splitting by source node gave
    # a triangular load: task 1 solved n-1 pairs, task n none, so the wall
    # time was bounded by the first source rather than pairs/threads. A few
    # chunks per thread leaves the scheduler slack for stragglers while still
    # amortizing the workspace copy each task makes.
    jobs = collect(Iterators.flatten(groups))
    isempty(jobs) && return nothing
    nchunks = cfg.parallelize ? min(length(jobs), 4 * Threads.nthreads()) : 1
    chunks = collect(Iterators.partition(jobs, cld(length(jobs), nchunks)))
    timers = [TimerOutput() for _ in chunks]
    @timeit CSTIMER[] "solve and accumulate pairs" if cfg.parallelize
        fetch.(map(i -> Threads.@spawn(task(chunks[i], timers[i])), eachindex(chunks)))
    else
        foreach(i -> task(chunks[i], timers[i]), eachindex(chunks))
    end
    foreach(t -> merge!(CSTIMER[], t), timers)
    nothing
end

function solve_pairs!(handle_pair, factor, solver::DirectSolver, matrix::SparseMatrixCSC{T},
                      groups, cfg, log, pair_label) where T

    jobs = collect(Iterators.flatten(groups))
    n = size(matrix, 1)

    for batch in Iterators.partition(jobs, solver.bs)
        rhs = zeros(T, n, length(batch))
        for (col, job) in enumerate(batch)
            rhs[job.comp_i, col] = -1
            rhs[job.comp_j, col] = 1
        end

        log && @debug("Solving a batch of $(length(batch)) pairs: " * pair_label(first(batch)) *
                      " to " * pair_label(last(batch)))
        lhs = solve_linear_system(factor, matrix, rhs; tol = residual_tolerance(cfg, T))

        function post(col)
            timer = TimerOutput()
            job = batch[col]
            voltages = lhs[:, col]
            voltages .-= voltages[job.comp_i]
            @timeit timer "postprocess" handle_pair(job, voltages)
            timer
        end

        timers = @timeit CSTIMER[] "postprocess pairs" if cfg.parallelize
            fetch.(map(col -> Threads.@spawn(post(col)), 1:length(batch)))
        else
            map(post, 1:length(batch))
        end
        foreach(t -> merge!(CSTIMER[], t), timers)
    end
    nothing
end

# TODO: In the pardiso case, we're not really constructing the factor
# So can we make this consistent?
function construct_cholesky_factor(matrix, ::CholmodSolver)
    T = eltype(matrix)
    factor = cholesky(matrix + sparse(T(10)*eps(T)*I,size(matrix)...))
    factor
end


"""
    get_num_pairs(ccs, fp, exclude_pairs, user_points = fp) -> (n, numbering)

Number of pairs the pairwise driver will solve without the resistance
shortcut: for every connected component, the pairs of distinct focal nodes in
it whose user ids are not in `exclude_pairs`. `numbering` maps each
`(node, node)` pair to its position in that count, which is how the log
labels `Solving pair k of n`.
"""
function get_num_pairs(ccs, fp::Vector{V}, exclude_pairs, user_points::Vector{V}=fp) where V

    num = 0
    d = Dict{Tuple{V,V}, V}()
    # Map graph node indices to user point IDs for exclude comparison
    g2u = Dict(fp[i] => user_points[i] for i in 1:length(fp))

    for (i,cc) in enumerate(ccs)
        sub_fp = component_points(cc, fp)
        l = lastindex(sub_fp)
        for ii = 1:l
            pt1 = sub_fp[ii]
            for jj = ii+1:l
                pt2 = sub_fp[jj]
                if (get(g2u, pt1, pt1), get(g2u, pt2, pt2)) in exclude_pairs
                    continue
                else
                    num += 1
                    d[(pt1, pt2)] = num
                end
            end
        end
    end
    num, d
end

function get_num_pairs_shortcut(ccs, fp::Vector{V}, exclude_pairs, user_points::Vector{V}=fp) where V

    num = 0
    d = Dict{Tuple{V,V}, V}()
    g2u = Dict(fp[i] => user_points[i] for i in 1:length(fp))

    for (i,cc) in enumerate(ccs)
        sub_fp = component_points(cc, fp)
        l = lastindex(sub_fp)
        l == 0 && continue
        for ii = 1:1
            pt1 = sub_fp[ii]
            for jj = ii+1:l
                pt2 = sub_fp[jj]
                if (get(g2u, pt1, pt1), get(g2u, pt2, pt2)) in exclude_pairs
                    continue
                else
                    num += 1
                    d[(pt1, pt2)] = num
                end
            end
        end
    end
    num, d
end
function smash_repeats!(ret, I)
    for i = 1:size(I,1)
        for j = i+1:size(I,1)
            push!(ret, (I[i], I[j], 0))
        end
    end
end

function smash_repeats!(resistances::Matrix{T}, I) where T
    for i = 1:size(I,1)
        for j = i+1:size(I,1)
            resistances[I[i], I[j]] = 0
            resistances[I[j], I[i]] = 0
        end
    end
end

"""
Calculate laplacian of the adjacency matrix of a graph
"""
function laplacian!(G::SparseMatrixCSC{T,V}) where {T,V}
    n = size(G, 1)
    s = Vector{eltype(G)}(undef,n)
    for i = 1:n
        s[i] = sum_off_diag(G, i)
        for j in nzrange(G, i)
            if i == G.rowval[j]
                G.nzval[j] = 0
            else
                G.nzval[j] = -G.nzval[j]
            end
        end
    end
    r = V(1):V(n)
    S = sparse(r, r, s)
    G + S
end

function sum_off_diag(G, i)
     sum = zero(eltype(G))
     for j in nzrange(G, i)
         if G.rowval[j] != i
             sum += G.nzval[j]
         end
     end
     sum
 end

"""
    residual_tolerance(cfg, T)

Largest acceptable true relative residual `‖Gv − b‖ / ‖b‖` for a solve in
precision `T`. `cfg.residual_tolerance` when the user set one; otherwise
`1e-4` in double and `1e-3` in single, since Float32 round-off puts `1e-4`
at the edge of what CG can reach (the fixed `1e-4` gate broke every single
precision run).
"""
residual_tolerance(cfg, ::Type{T}) where {T} =
    cfg.residual_tolerance > 0 ? cfg.residual_tolerance :
    T === Float32 ? TOL_SINGLE : TOL_DOUBLE

"""
    solve_linear_system(G, curr, M; tol)

Preconditioned CG. Krylov's stopping test is on the *preconditioned* residual
`sqrt(r' M⁻¹ r)`, while the result is accepted on the true 2-norm relative
residual. When the AMG preconditioner is a poor fit for a component the two
diverge: CG reports success in its own norm while the true residual is still
above `tol` (issue #470). Rather than abort a multi-hour job on one such
solve, tighten the stopping tolerances and retry; solves that pass first time
are unaffected.
"""
function solve_linear_system(
            G::SparseMatrixCSC{T,V},
            curr::Vector{T}, M; tol = TOL_DOUBLE)::Vector{T} where {T,V}
    # Krylov's own default for atol; in Float32 it is ~3e-4 and is what
    # actually stops CG, so it has to be tightened along with rtol.
    atol = sqrt(eps(T))
    rtol = T === Float32 ? T(1e-5) : T(1e-6)
    local v, stats, residual
    for attempt in 1:3
        v, stats = Krylov.cg(G, curr, M=M, ldiv=true, atol=atol, rtol=rtol, itmax=100_000)
        residual = norm(G*v .- curr) / norm(curr)
        residual < tol && return v
        attempt < 3 && @warn("CG stopped after $(stats.niter) iterations ($(stats.status)) " *
                             "but the relative residual $residual exceeds $tol; " *
                             "retrying with rtol = $(rtol / 100)")
        atol /= 100
        rtol /= 100
    end
    error("CG solver did not converge: relative residual $residual exceeds tolerance $tol " *
          "after 3 attempts (last rtol = $(rtol * 100), $(stats.niter) iterations, $(stats.status))")
end

"""
    refine_columns!(lhs, factor, matrix, rhs, tol, name)

Check every column of a direct solve against `tol`, applying iterative
refinement with the existing factor where needed. In single precision a
Cholesky solve of an ill-conditioned Laplacian lands at 1e-3..1e-2, well short
of what CG reaches with its 1e-5 rtol; each step costs one triangular solve
and recovers several digits. Double precision residuals are far below the
target and take no steps. Shared by the CHOLMOD and Accelerate backends.
"""
function refine_columns!(lhs, factor, matrix, rhs, tol, name)
    target = min(tol, 1e-5)
    for col = 1:size(rhs, 2)
        x = view(lhs, :, col)
        b = view(rhs, :, col)
        residual = norm(matrix*x .- b) / norm(b)
        steps = 0
        while residual >= target && steps < 2
            x .+= factor \ (b .- matrix*x)
            residual = norm(matrix*x .- b) / norm(b)
            steps += 1
        end
        residual < tol || error("$name solver residual $residual exceeds tolerance $tol for column $col after $steps refinement steps")
    end
    lhs
end

function solve_linear_system(factor::SuiteSparse.CHOLMOD.Factor, matrix, rhs; tol = TOL_DOUBLE)
    refine_columns!(factor \ rhs, factor, matrix, rhs, tol, "CHOLMOD")
end

function postprocess(output, component_data, shortcut, cfg)


    orig_pts = output.orig_pts

    # Shortcut flags and data
    get_shortcut_resistances = shortcut.get_shortcut_resistances

    if get_shortcut_resistances
        update_voltmatrix!(shortcut, output, component_data)
        return nothing
    end

    # Nothing below produces a result unless some map output was requested;
    # without this, every pair still computed its node currents and a
    # full-grid current map only to discard them (~1 GB allocated per pair
    # at 1M cells, a fifth of the run time).
    (cfg.write_volt_maps || cfg.write_cur_maps ||
     cfg.write_cum_cur_map_only || cfg.write_max_cur_maps) || return nothing

    name = "_$(orig_pts[1])_$(orig_pts[2])"

    if cfg.write_volt_maps
        write_volt_maps(name, output, component_data, cfg)
    end

    # TODO: Even though this function is called write_cur_maps
    # actually writing the calculated maps depends on some options.
    write_cur_maps(name, output, component_data, [-9999.], cfg)
    nothing
end

function update_voltmatrix!(shortcut, output, component_data)

    # Data
    voltmatrix = shortcut.voltmatrix
    c = output.points
    cc = component_data.cc
    voltages = output.voltages
    r = output.resistance
    j = output.col

    for i = 2:size(c, 1)
        ind = component_index(cc, c[i])
        if ind !== nothing
            voltageAtPoint = voltages[ind]
            voltageAtPoint = 1 - (voltageAtPoint/r)
            voltmatrix[i,j] = voltageAtPoint
        end
    end
end


function update_shortcut_resistances!(anchor, sc, resistances, points, comp)

    # Data
    voltmatrix = sc.voltmatrix
    shortcut = sc.shortcut_res

    check = map(x -> in_component(comp, x), points)
    l = size(resistances, 1)
    for pointx = 1:l
        if check[pointx]
            R1x = resistances[anchor, pointx]
            if R1x != -1
                shortcut[pointx, anchor] = shortcut[anchor, pointx] = R1x
                for point2 = pointx:l
                    if check[point2]
                        R12 = resistances[anchor, point2]
                        if R12 != -1
                            if R1x != RESISTANCE_INVALID
                                shortcut[anchor, point2] = shortcut[point2, anchor] = R12
                                Vx = voltmatrix[pointx, point2]
                                R2x = 2*R12*Vx + R1x - R12
                                if shortcut[point2, pointx] != RESISTANCE_INVALID
                                    shortcut[point2, pointx] = shortcut[pointx, point2] = R2x
                                end
                            else
                                shortcut[pointx, :] = shortcut[:, pointx] = RESISTANCE_INVALID
                            end
                        end
                    end
                end
            end
        end
    end
end
