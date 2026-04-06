struct Cumulative{T,V}
    cum_curr::Matrix{T}
    max_curr::Matrix{T}
	cum_branch_curr::Vector{T}
	cum_node_curr::Vector{T}
	coords::Vector{Tuple{V,V}}
    lock::ReentrantLock
end

struct GraphProblem{T,V,W}
    G::SparseMatrixCSC{T,V}
    cc::Vector{Vector{V}}
    points::Vector{V}
    user_points::Vector{V}
    exclude_pairs::Vector{Tuple{V,V}}
    nodemap::Matrix{V}
    polymap::Matrix{V}
    hbmeta::RasterMeta
    cellmap::Matrix{T}
    cum::Cumulative{T}
    solver::W
end

struct ComponentData{T,V}
    cc::Vector{V}
    matrix::SparseMatrixCSC{T,V}
    local_nodemap::Matrix{V}
    hbmeta::RasterMeta
    cellmap::Matrix{T}
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

struct Shortcut{T}
    get_shortcut_resistances::Bool
    voltmatrix::Matrix{T}
    shortcut_res::Matrix{T}
end

abstract type Solver end

struct CholmodSolver <: Solver
    bs::Int
end

struct AMGSolver <: Solver
end

struct PardisoSolver <: Solver
    bs::Int
end
"""
Core kernel of Circuitscape - used to solve several pairs

Input:
* data::GraphProblem
"""
function single_ground_all_pairs(prob::GraphProblem{T,V,W}, flags, cfg, log = true) where {T,V,W}
    solve(prob, prob.solver, flags, cfg, log)
end

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
    else
        error("Unknown solver: $s")
    end
end

function solve(prob::GraphProblem{T,V}, ::AMGSolver, flags, cfg, log)::Matrix{T} where {T,V}

    # Data
    a = prob.G
    cc = prob.cc
    points = prob.points
    exclude = prob.exclude_pairs
    nodemap = prob.nodemap
    polymap = prob.polymap
    orig_pts = prob.user_points
    hbmeta = prob.hbmeta
    cellmap = prob.cellmap

    # Flags
    outputflags = flags.outputflags
    is_raster = flags.is_raster
    write_volt_maps = outputflags.write_volt_maps
    write_cur_maps = outputflags.write_cur_maps
    write_cum_cur_map_only = outputflags.write_cum_cur_map_only
    write_max_cur_maps = outputflags.write_max_cur_maps

    # Get number of focal points
    numpoints = size(points, 1)

    # Cumulative currents

    cum = prob.cum

    @info("Graph has $(size(a,1)) nodes, $numpoints focal points and $(length(cc)) connected components")

    num_pairs, pair_numbers = get_num_pairs(cc, points, exclude, orig_pts)
    log && @info("Total number of pair solves = $num_pairs")

    # Initialize pairwise resistance
    resistances = -1 * ones(T, numpoints, numpoints)::Matrix{T}
    voltmatrix = zeros(T, size(resistances))::Matrix{T}
    shortcut_res = deepcopy(resistances)::Matrix{T}

    # Get a vector of connected components
    comps = getindex.([a], cc, cc)

    get_shortcut_resistances = false
    if is_raster && !write_volt_maps && !write_cur_maps &&
            !write_cum_cur_map_only && !write_max_cur_maps &&
            isempty(exclude)
        get_shortcut_resistances = true
        @info("Triggering resistance calculation shortcut")
        num_pairs, pair_numbers = get_num_pairs_shortcut(cc, points, exclude, orig_pts)
        @info("Total number of pair solves has been reduced to $num_pairs ")
    end
    shortcut = Shortcut(get_shortcut_resistances, voltmatrix, shortcut_res)

    for (cid, comp) in enumerate(cc)

        # Subset of points relevant to CC
        csub = filter(x -> x in comp, points) |> unique

        if isempty(csub)
            continue
        end

        # Conductance matrix corresponding to CC
        matrix = comps[cid]

        # Regularization step
        matrix.nzval .+= eps(eltype(matrix)) * norm(matrix.nzval)

        # Construct preconditioner *once* for every CC
        t1 = @elapsed P = aspreconditioner(smoothed_aggregation(matrix))
        @info("Time taken to construct preconditioner = $t1 seconds")

        # Get local nodemap for CC - useful for output writing
        t2 = @elapsed local_nodemap = construct_local_node_map(nodemap, comp, polymap)
        @info("Time taken to construct local nodemap = $t2 seconds")

        component_data = ComponentData(comp, matrix, local_nodemap, hbmeta, cellmap)

        function solve_pairs_for_point(point_idx)
            # Each task needs its own workspace (scratch vectors are mutable)
            ml = P.ml
            local_P = aspreconditioner(AlgebraicMultigrid.MultiLevel(
                ml.levels, ml.final_A, ml.coarse_solver,
                ml.presmoother, ml.postsmoother, deepcopy(ml.workspace)))

            results = Vector{Tuple{V,V,T}}()

            src_node = csub[point_idx]
            comp_i = findfirst(isequal(src_node), comp)
            comp_i === nothing && error("Node $src_node not found in component")
            comp_i = V(comp_i)
            src_indices = findall(x -> x == src_node, points)
            smash_repeats!(results, src_indices)

            # Iteration space through all possible pairs
            pair_range = point_idx+1:size(csub, 1)
            if Threads.nthreads() > 1
                for pair_idx in pair_range
                    dst_node = csub[pair_idx]
                    haskey(pair_numbers, (src_node, dst_node)) && @debug("Scheduling pair $(pair_numbers[(src_node, dst_node)]) of $num_pairs to be solved")
                end
            end

            # Loop through all possible pairs
            for pair_idx in pair_range

                dst_node = csub[pair_idx]
                comp_j = findfirst(isequal(dst_node), comp)
                comp_j === nothing && error("Node $dst_node not found in component")
                comp_j = V(comp_j)
                dst_indices = findall(x -> x == dst_node, points)

                if src_node == dst_node
                    continue
                end

                # Check if all index combinations are excluded
                any_included = false
                for c_i in src_indices, c_j in dst_indices
                    if (orig_pts[c_i], orig_pts[c_j]) ∉ exclude
                        any_included = true
                        break
                    end
                end
                !any_included && continue

                # Solve once per (src_node, dst_node) pair
                current = zeros(T, size(matrix, 1))
                current[comp_i] = -1
                current[comp_j] = 1

                log && haskey(pair_numbers, (src_node, dst_node)) && @debug("Solving pair $(pair_numbers[(src_node, dst_node)]) of $num_pairs")
                solve_time = @elapsed voltages = solve_linear_system(matrix, current, local_P)
                @debug("Time taken to solve linear system = $solve_time seconds")

                voltages .= voltages .- voltages[comp_i]
                resistance = voltages[comp_j] - voltages[comp_i]

                # Store result for all non-excluded index combinations
                for c_i in src_indices
                    for c_j in dst_indices
                        if (orig_pts[c_i], orig_pts[c_j]) in exclude
                            continue
                        end
                        push!(results, (c_i, c_j, resistance))
                        if get_shortcut_resistances
                            resistances[c_i, c_j] = resistance
                            resistances[c_j, c_i] = resistance
                        end
                        output = Output(points, voltages, (orig_pts[c_i], orig_pts[c_j]),
                                        (comp_i, comp_j), resistance, V(c_j), cum)
                        postprocess(output, component_data, flags, shortcut, cfg)
                    end
                end
            end

        results
        end

        if get_shortcut_resistances
            idx = findfirst(isequal(csub[1]), points)
            idx === nothing && error("Focal point $(csub[1]) not found in points list")
            solve_pairs_for_point(1)
            update_shortcut_resistances!(idx, shortcut, resistances, points, comp)
        else
            is_parallel = cfg.parallelize
            if is_parallel
                all_results = fetch.(map(pt -> Threads.@spawn(solve_pairs_for_point(pt)), 1:size(csub,1)))
            else
                all_results = map(solve_pairs_for_point, 1:size(csub,1))
            end

            # Set all resistances
            for task_result in all_results
                for (ci, cj, rv) in task_result
                    resistances[ci, cj] = rv
                    resistances[cj, ci] = rv
                end
            end
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

struct CholmodNode{T}
    cc_idx::Tuple{T,T}
    points_idx::Tuple{T,T}
end

function solve(prob::GraphProblem{T,V}, solver::Union{CholmodSolver, PardisoSolver}, flags,
                                  cfg, log) where {T,V}

    # Data
    a = prob.G
    cc = prob.cc
    points = prob.points
    exclude = prob.exclude_pairs
    nodemap = prob.nodemap
    polymap = prob.polymap
    orig_pts = prob.user_points
    hbmeta = prob.hbmeta
    cellmap = prob.cellmap

    # Flags
    outputflags = flags.outputflags
    is_raster = flags.is_raster
    write_volt_maps = outputflags.write_volt_maps
    write_cur_maps = outputflags.write_cur_maps
    write_cum_cur_map_only = outputflags.write_cum_cur_map_only
    write_max_cur_maps = outputflags.write_max_cur_maps

    # Cumulative current map
    cum = prob.cum

    # Batchsize
    batch_size = solver.bs

    # Get number of focal points
    numpoints = size(points, 1)

    @info("Graph has $(size(a,1)) nodes, $numpoints focal points and $(length(cc)) connected components")

    num_pairs, _ = get_num_pairs(cc, points, exclude, orig_pts)
    log && @info("Total number of pair solves = $num_pairs")

    # Initialize pairwise resistance
    resistances = -1 * ones(eltype(a), numpoints, numpoints)
    voltmatrix = zeros(eltype(a), size(resistances))
    shortcut_res = -1 * ones(eltype(a), size(resistances))

    # Get a vector of connected components
    comps = getindex.([a], cc, cc)

    get_shortcut_resistances = false
    if is_raster && !write_volt_maps && !write_cur_maps &&
            !write_cum_cur_map_only  && !write_max_cur_maps &&
            isempty(exclude)
        get_shortcut_resistances = true
        @info("Triggering resistance calculation shortcut")
        num_pairs, _ = get_num_pairs_shortcut(cc, points, exclude, orig_pts)
        @info("Total number of pair solves has been reduced to $num_pairs ")
    end
    shortcut = Shortcut(get_shortcut_resistances, voltmatrix, shortcut_res)

    for (cid, comp) in enumerate(cc)

        # Subset of points relevant to CC
        csub = filter(x -> x in comp, points) |> unique

        if isempty(csub)
            continue
        end

        # Conductance matrix corresponding to CC
        matrix = comps[cid]

        t = @elapsed factor = construct_cholesky_factor(matrix, solver)

        # Get local nodemap for CC - useful for output writing
        t2 = @elapsed local_nodemap = construct_local_node_map(nodemap, comp, polymap)
        @info("Time taken to construct local nodemap = $t2 seconds")

        component_data = ComponentData(comp, matrix, local_nodemap, hbmeta, cellmap)

        cholmod_batch = CholmodNode[]

        # Build batch of pairs to solve
        function build_cholmod_batch(point_idx)

            src_node = csub[point_idx]
            comp_i_raw = findfirst(isequal(src_node), comp)
            comp_i_raw === nothing && error("Node $src_node not found in component")
            comp_i = V(comp_i_raw)
            src_indices = findall(x -> x == src_node, points)
            smash_repeats!(resistances, src_indices)

            # Iteration space through all possible pairs
            pair_range = point_idx+1:size(csub, 1)

            # Loop through all possible pairs
            for pair_idx in pair_range

                dst_node = csub[pair_idx]
                comp_j_raw = findfirst(isequal(dst_node), comp)
                comp_j_raw === nothing && error("Node $dst_node not found in component")
                comp_j = V(comp_j_raw)
                dst_indices = findall(x -> x == dst_node, points)

                if src_node == dst_node
                    continue
                end

                # Forget excluded pairs
                for c_i in src_indices, c_j in dst_indices
                    if (orig_pts[c_i], orig_pts[c_j]) in exclude
                        continue
                    else
                        push!(cholmod_batch,
                          CholmodNode((comp_i, comp_j), (V(c_i), V(c_j))))
                    end
                end
            end
        end

        function postprocess_pair(batch_idx, batch_range, lhs)
            batch_pos = batch_range[batch_idx]
            output = Output(points, lhs[:,batch_idx],
                (orig_pts[cholmod_batch[batch_pos].points_idx[1]],
                orig_pts[cholmod_batch[batch_pos].points_idx[2]]),
                cholmod_batch[batch_pos].cc_idx,
                lhs[cholmod_batch[batch_pos].cc_idx[2], batch_idx] - lhs[cholmod_batch[batch_pos].cc_idx[1], batch_idx],
                V(cholmod_batch[batch_pos].points_idx[2]), cum)
            postprocess(output, component_data, flags, shortcut, cfg)
        end
        if get_shortcut_resistances
            idx = findfirst(isequal(csub[1]), points)
            idx === nothing && error("Focal point $(csub[1]) not found in points list")
            build_cholmod_batch(1)
        else
            build_cholmod_batch.(1:size(csub, 1))
        end

        num_batched_pairs = length(cholmod_batch)

        for st in 1:batch_size:num_batched_pairs

            batch_range = st + batch_size <= num_batched_pairs ?
                            (st:(st+batch_size-1)) : (st:num_batched_pairs)

            @debug("Solving points $(batch_range.start) to $(batch_range.stop)")

            rhs = zeros(eltype(matrix), size(matrix, 1), length(batch_range))

            for (col, batch_pos) in enumerate(batch_range)
                node = cholmod_batch[batch_pos]
                rhs[node.cc_idx[1], col] = -1
                rhs[node.cc_idx[2], col] = 1
            end

            lhs = solve_linear_system(factor, matrix, rhs)

            # Normalisation step
            for (col, batch_pos) in enumerate(batch_range)
                ref_node = cholmod_batch[batch_pos].cc_idx[1]
                ref_volt = lhs[ref_node, col]
                for row = 1:size(matrix, 1)
                    lhs[row, col] = lhs[row, col] - ref_volt
                end
            end

            is_parallel = cfg.parallelize
            if is_parallel
                fetch.(map(bi -> Threads.@spawn(postprocess_pair(bi, batch_range, lhs)), 1:length(batch_range)))
            else
                map(bi -> postprocess_pair(bi, batch_range, lhs), 1:length(batch_range))
            end

            for (col, batch_pos) in enumerate(batch_range)
                coords = cholmod_batch[batch_pos].points_idx
                resistance = lhs[cholmod_batch[batch_pos].cc_idx[2], col] -
                            lhs[cholmod_batch[batch_pos].cc_idx[1], col]
                resistances[coords...] = resistance
                resistances[reverse(coords)...] = resistance
            end
        end

        if get_shortcut_resistances
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

# TODO: In the pardiso case, we're not really constructing the factor
# So can we make this consistent?
function construct_cholesky_factor(matrix, ::CholmodSolver)
    T = eltype(matrix)
    t = @elapsed factor = cholesky(matrix + sparse(T(10)*eps(T)*I,size(matrix)...))
    @info("Time taken to construct cholesky factor = $t")
    factor
end


"""
Returns all possible pairs to solve.

Input:
* ccs::Vector{Vector{Int}} - vector of connected components
* fp::Vector{Int} - vector of focal points
* exclude_pairs::Vector{Tuple{Int,Int}} - vector of point pairs (tuples) to exclude

Output:
* n - total number of pairs
"""
function get_num_pairs(ccs, fp::Vector{V}, exclude_pairs, user_points::Vector{V}=fp) where V

    num = 0
    d = Dict{Tuple{V,V}, V}()
    # Map graph node indices to user point IDs for exclude comparison
    g2u = Dict(fp[i] => user_points[i] for i in 1:length(fp))

    for (i,cc) in enumerate(ccs)
        sub_fp = filter(x -> x in cc, fp) |> unique
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
        sub_fp = filter(x -> x in cc, fp) |> unique
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

function solve_linear_system(
            G::SparseMatrixCSC{T,V},
            curr::Vector{T}, M)::Vector{T} where {T,V}
    v, stats = Krylov.cg(G, curr, M=M, ldiv=true, rtol=T(1e-6), itmax=100_000)
    residual = norm(G*v .- curr) / norm(curr)
    residual < 1e-4 || error("CG solver did not converge: relative residual $residual exceeds tolerance 1e-4")
    v
end


function solve_linear_system(factor::SuiteSparse.CHOLMOD.Factor, matrix, rhs)
    lhs = factor \ rhs
    for col = 1:size(rhs, 2)
        residual = norm(matrix*lhs[:,col] .- rhs[:,col]) / norm(rhs[:,col])
        residual < 1e-4 || error("CHOLMOD solver residual $residual exceeds tolerance 1e-4 for column $col")
    end
    lhs
end

function postprocess(output, component_data, flags, shortcut, cfg)


    voltages = output.voltages
    matrix = component_data.matrix
    local_nodemap = component_data.local_nodemap
    hbmeta = component_data.hbmeta
    orig_pts = output.orig_pts

    # Shortcut flags and data
    get_shortcut_resistances = shortcut.get_shortcut_resistances

    if get_shortcut_resistances
        update_voltmatrix!(shortcut, output, component_data)
        return nothing
    end


    name = "_$(orig_pts[1])_$(orig_pts[2])"

    if flags.outputflags.write_volt_maps
        t = @elapsed write_volt_maps(name, output, component_data, flags, cfg)
        @debug("Time taken to write voltage maps = $t seconds")
    end

    # TODO: Even though this function is called write_cur_maps
    # actually writing the calculated maps depends on some flags.
    t = @elapsed write_cur_maps(name, output, component_data,
                                [-9999.], flags, cfg)
    @debug("Time taken to calculate current maps = $t seconds")
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
        ind = findfirst(isequal(c[i]), cc)
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

    check = map(x -> x in comp, points)
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
