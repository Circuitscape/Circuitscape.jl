"""
    AdvancedProblem

An advanced-mode problem: the graph Laplacian `G`, its connected components
`cc`, the [`Geometry`](@ref) mapping nodes to output coordinates, the source
and ground vectors over the nodes (`sources`, `grounds`, with the finite
grounds separated out in `finitegrounds` for the solver) and the solver.
Built by [`build_problem`](@ref) in advanced mode and, one focal node at a
time, by the one-to-all / all-to-one driver; solved by
[`advanced_kernel`](@ref).
"""
struct AdvancedProblem{T,V,W,Geom<:Geometry}
    G::SparseMatrixCSC{T,V}
    cc::Vector{Vector{V}}
    geometry::Geom
    sources::Vector{T}
    grounds::Vector{T}
    finitegrounds::Vector{T}
    solver::W
end

"""
    sources_and_grounds(geometry, source_map, ground_map, G, cfg, conflict_policy = policy(cfg);
                        cc = nothing)

Source and ground vectors over the graph nodes from the user's source and
ground inputs, with conflicts between the two resolved by `conflict_policy`.
For a raster the inputs are grids read through the nodemap; for a network
they are `(node, value)` lists. Given the connected components `cc`, a
component that had sources and grounds but keeps none of one after the
conflicts are resolved is an error, as in Python (`compute.resolve_conflicts`
is run per component there). One-to-all / all-to-one do not pass `cc`: a focal
node whose solve is left with nothing to solve reports `-1`, as when it is the
only focal node, rather than ending the whole run.
"""
function sources_and_grounds(geometry::Geometry, source_map, ground_map, G, cfg,
                             conflict_policy = policy(cfg); cc = nothing)

    # Initialize sources and grounds
    sources = zeros(eltype(G), size(G, 1))
    grounds = zeros(eltype(G), size(G, 1))

    fill_sources_and_grounds!(sources, grounds, source_map, ground_map, cfg, geometry)
    had_sources = map(!iszero, sources)
    had_grounds = map(!iszero, grounds)

    sources, grounds, finitegrounds =
        resolve_conflicts(sources, grounds, conflict_policy)

    cc === nothing ||
        check_conflicts_resolved(cc, had_sources, had_grounds, sources, grounds)

    sources, grounds, finitegrounds
end

# Python resolves conflicts per component and raises when a component that
# had both sources and grounds is left without one of them; a component that
# never had both is skipped by the kernel instead.
function check_conflicts_resolved(cc, had_sources, had_grounds, sources, grounds)
    for c in cc
        any(i -> had_sources[i], c) && any(i -> had_grounds[i], c) || continue
        any(i -> !iszero(sources[i]), c) || throw(ErrorException(
            "All sources conflicted with grounds and were removed. There is nothing to solve."))
        any(i -> !iszero(grounds[i]), c) || throw(ErrorException(
            "All grounds conflicted with sources and were removed. There is nothing to solve."))
    end
    nothing
end

function fill_sources_and_grounds!(sources, grounds, source_map, ground_map, cfg,
                                   geometry::RasterGeometry{T,V}) where {T,V}
    nodemap = geometry.nodemap
    (i1, j1, v1) = begin _I = findall(!iszero, source_map); getindex.(_I, 1), getindex.(_I, 2), source_map[_I] end
    (i2, j2, v2) = begin _I = findall(!iszero, ground_map); getindex.(_I, 1), getindex.(_I, 2), ground_map[_I] end
    for i = 1:size(i1, 1)
        v = V(nodemap[i1[i], j1[i]])
        if v != 0
            sources[v] += v1[i]
        end
    end
    for i = 1:size(i2, 1)
        v = V(nodemap[i2[i], j2[i]])
        if v != 0
            grounds[v] += v2[i]
        end
    end
    nothing
end

function fill_sources_and_grounds!(sources, grounds, source_map, ground_map, cfg,
                                   ::NetworkGeometry{V}) where {V}
    # Invert into a copy: `ground_map` is the caller's data and may be used
    # to build the problem again.
    ground_list = cfg.ground_file_is_resistances ?
                  hcat(ground_map[:,1], 1 ./ ground_map[:,2]) : ground_map
    assign_network_nodes!(sources, source_map, "source", V)
    assign_network_nodes!(grounds, ground_list, "ground", V)
    nothing
end

# Set `values[node] = value` for each `(node, value)` row. A node id that is
# not in the graph is skipped with a warning, as Circuitscape 4 matches the
# ids against the graph's node names and silently ignores the rest
# (compute.py, advanced_module).
function assign_network_nodes!(values, list, what, ::Type{V}) where {V}
    n = length(values)
    for r in axes(list, 1)
        node = list[r, 1]
        if 1 <= node <= n && isinteger(node)
            values[V(node)] = list[r, 2]
        else
            id = isinteger(node) ? Int(node) : node
            @warn("Ignoring $what at node $id: not a node of the graph (nodes are 1 to $n)")
        end
    end
    values
end

# Resolve the nodes that are both a source and a ground by `policy`, then
# drop any direct ground still under a source. Python's node-level pass
# (`compute.resolve_conflicts`) differs on purpose in two ways that only show
# where its read-time cell pass did not already resolve the conflict (network
# inputs, cells merged into one node by polygons): its `rmvall` removes only
# the source, and it takes the finite grounds before resolving, so a ground
# removed by `rmvgnd` stays on the solver diagonal. Here `rmvall` removes both,
# as its name and the Circuitscape 4 GUI ("Remove both source and ground")
# say, and the finite grounds are the grounds that remain.
function resolve_conflicts(sources::Vector{T},
                            grounds::Vector{T}, policy) where T

    conflicts = map((s, g) -> !iszero(s) && !iszero(g), sources, grounds)

    if policy == :rmvsrc || policy == :rmvall
        sources[conflicts] .= 0
    end
    if policy == :rmvgnd || policy == :rmvall
        grounds[conflicts] .= 0
    end

    # A direct ground under a source (of either sign) is removed, as in Python.
    infconflicts = map((s, g) -> !iszero(s) && g == T(Inf), sources, grounds)
    grounds[infconflicts] .= 0

    finitegrounds = map(x -> x < T(Inf) ? x : T(0), grounds)
    if !any(!iszero, finitegrounds)
        finitegrounds = T[-9999]
    end

    sources, grounds, finitegrounds
end

"""
    advanced_kernel(prob, cfg; check_node, name, accumulate_currents)
        -> (voltages, current_map, solver_called)

Solve every connected component that has both a source and a ground, and
write the voltage and current maps the options ask for, suffixed `name`.
`voltages` is over all graph nodes (zero in components not solved); the
current map is the per-cell sum over components for a raster (empty for a
network), filled when `accumulate_currents` (by default when either current
map option is set). One-to-all / all-to-one solve
one focal node at a time and pass it as `check_node` to skip every other
component.
"""
function advanced_kernel(prob::AdvancedProblem{T,V}, cfg;
                         check_node = V(-1), name = "",
                         accumulate_currents = cfg.write_cur_maps || cfg.write_cum_cur_map_only) where {T,V}

    # Data
    G = prob.G
    geometry = prob.geometry
    sources = prob.sources
    grounds = prob.grounds
    finitegrounds = prob.finitegrounds
    cc = prob.cc

    # Options
    write_v_maps = cfg.write_volt_maps
    write_c_maps = cfg.write_cur_maps
    write_cum_cur_map_only = cfg.write_cum_cur_map_only

    f_local = Vector{eltype(G)}()
    solver_called = false
    voltages = zeros(eltype(G), size(G, 1))
    outvolt = write_v_maps ? alloc_map(geometry) : zeros(0, 0)
    outcurr = alloc_map(geometry)

    for c in cc

        if check_node != -1 && !(check_node in c)
            continue
        end

        # A component spanning the whole graph is the graph: no local copies.
        whole = length(c) == size(G, 1)
        a_local = whole ? G : G[c,c]
        s_local = sources[c]
        g_local = grounds[c]

        # A component is solved when it has a source and a ground; count
        # them as Python does, so sources of +1 and -1 are not skipped.
        if !any(!iszero, s_local) || !any(!iszero, g_local)
            continue
        end

        if finitegrounds != [-9999.]
            f_local = finitegrounds[c]
        else
            f_local = finitegrounds
        end

        voltages[c] .+= multiple_solver(cfg, prob.solver, a_local, s_local, g_local, f_local)
        solver_called = true

        local_geometry = whole && is_raster(geometry) ? geometry : restrict(geometry, c)
        write_v_maps && accum_voltages!(outvolt, voltages[c], local_geometry)
        accumulate_currents && accum_currents!(outcurr, a_local, voltages[c], f_local, local_geometry)
    end

    write_v_maps && write_advanced_volt_map(name, voltages, outvolt, geometry, cfg)

    # Issue 342: in one-to-all / all-to-one every focal node is active on each
    # solve (one as source, the rest as grounds), so all of them are zeroed.
    if (is_onetoall(cfg) || is_alltoone(cfg)) && cfg.set_focal_node_currents_to_zero
        active = V[n for n in eachindex(sources) if sources[n] != 0 || grounds[n] != 0]
        zero_focal_nodes!(outcurr, geometry, Set(active))
    end

    # One-to-all / all-to-one write a map per focal node with
    # `write_cur_maps` unless `write_cum_cur_map_only`, as Circuitscape 4 does;
    # their cumulative map is written by the driver. Advanced mode has a single
    # map, which is its cumulative map, so either option writes it (as in 6.x;
    # Circuitscape 4 ignores `write_cum_cur_map_only` here and writes nothing).
    per_solve = if is_onetoall(cfg) || is_alltoone(cfg)
        write_c_maps && !write_cum_cur_map_only
    else
        write_c_maps || write_cum_cur_map_only
    end
    per_solve && write_advanced_cur_map(name, voltages, outcurr, G, finitegrounds, geometry, cfg)

    voltages, outcurr, solver_called
end

# The value an advanced-mode run returns: every node's voltage for a network,
# the voltage grid for a raster (`[-1]` when nothing could be solved).
advanced_result(voltages, solver_called, geometry::NetworkGeometry) =
    [geometry.nodes voltages]

function advanced_result(voltages::Vector{T}, solver_called, geometry::RasterGeometry) where T
    solver_called || return fill(T(-1), 1, 1)
    _create_voltage_map(voltages, geometry.nodemap, geometry.hbmeta)
end

# Zero the current of every cell whose node is in `nodes`; a network has no
# current map to zero.
zero_focal_nodes!(outcurr, geometry::RasterGeometry, nodes) =
    zero_focal_cells!(outcurr, geometry.nodemap, nodes)
zero_focal_nodes!(outcurr, ::NetworkGeometry, nodes) = outcurr

# Advanced mode writes one voltage and one current map per run from the
# whole-graph voltages: a grid for a raster, node (and branch) lists for a
# network. The output options (log transform, null cells to nodata) are
# applied by the grid writers; the current grid is copied first because the
# raw `outcurr` is also returned, and one-to-all accumulates it.
write_advanced_volt_map(name, voltages, outvolt, geometry::RasterGeometry, cfg) =
    write_voltage_grid(outvolt, name, cfg, geometry)
write_advanced_volt_map(name, voltages, outvolt, geometry::NetworkGeometry, cfg) =
    write_volt_maps(geometry, name, voltages, cfg)

write_advanced_cur_map(name, voltages, outcurr, G, finitegrounds, geometry::RasterGeometry, cfg) =
    write_current_grid(copy(outcurr), name, cfg, geometry)
write_advanced_cur_map(name, voltages, outcurr, G, finitegrounds, geometry::NetworkGeometry, cfg) =
    write_cur_maps(name, voltages, ComponentData(geometry.nodes, G, geometry), finitegrounds, cfg)


# `A + Diagonal(d)`, zeros dropped as by sparse `+`. Written into `A` when the
# caller owns it (`inplace`) and every diagonal entry is already stored.
function add_diagonal!(A::SparseMatrixCSC{T,V}, d, inplace) where {T,V}
    n = size(A, 1)
    rowval, nzval = rowvals(A), nonzeros(A)
    if inplace && all(i -> any(j -> rowval[j] == i, nzrange(A, i)), 1:n)
        for i = 1:n, j in nzrange(A, i)
            rowval[j] == i && (nzval[j] += d[i])
        end
        return dropzeros!(A)
    end
    # Diagonal built with the index type of `A`: `spdiagm` is always `Int`
    # indexed and would promote a 32-bit matrix to 64-bit indices.
    r = V(1):V(n)
    A + sparse(r, r, d, n, n)
end

function multiple_solver(cfg, solver, a::SparseMatrixCSC{T,V}, sources, grounds, finitegrounds) where {T,V}

    # `a` is never written to: it is copied at most once, by dropping the
    # direct grounds or adding the finite ones, and otherwise solved as is.
    infgrounds = findall(x -> x == Inf, grounds)
    deleteat!(sources, infgrounds)
    asolve = a
    if !isempty(infgrounds)
        r = collect(V(1):V(size(a, 1)))
        deleteat!(r, infgrounds)
        asolve = a[r, r]
    end
    if finitegrounds[1] != -9999
        fg = isempty(infgrounds) ? finitegrounds : deleteat!(copy(finitegrounds), infgrounds)
        asolve = add_diagonal!(asolve, fg, asolve !== a)
    end

    volt = multiple_solve(cfg, solver, asolve, sources)

    # Replace the inf with 0
    voltages = zeros(eltype(a), length(volt) + length(infgrounds))
    k = 1
    for i = 1:size(voltages, 1)
        if i in infgrounds
            voltages[i] = 0
        else
            #voltages[i] = volt[1][k]
            voltages[i] = volt[k]
            k += 1
        end
    end
    voltages
end

function multiple_solve(cfg, s::AMGSolver, matrix::SparseMatrixCSC{T,V}, sources::Vector{T}) where {T,V}
    # Pin the coarse solver and smoothers for the same reason as in core.jl:
    # CG needs a symmetric positive definite preconditioner.
    M = aspreconditioner(smoothed_aggregation(matrix;
                                              coarse_solver = AlgebraicMultigrid.Pinv,
                                              presmoother = AlgebraicMultigrid.GaussSeidel(),
                                              postsmoother = AlgebraicMultigrid.GaussSeidel()))
    solve_linear_system(matrix, sources, M; tol = residual_tolerance(cfg, T))
end

function multiple_solve(cfg, s::CholmodSolver, matrix::SparseMatrixCSC{T,V}, sources::Vector{T}) where {T,V}
    factor = construct_cholesky_factor(matrix, s)
    solve_linear_system(factor, matrix, sources; tol = residual_tolerance(cfg, T))
end

function multiple_solve(cfg, s::PardisoSolver, matrix::SparseMatrixCSC{T,V}, sources::Vector{T}) where {T,V}
    factor = construct_cholesky_factor(matrix, s)
    solve_linear_system(factor, matrix, sources; tol = residual_tolerance(cfg, T))
end

function multiple_solve(cfg, s::AccelerateSolver, matrix::SparseMatrixCSC{T,V}, sources::Vector{T}) where {T,V}
    factor = construct_cholesky_factor(matrix, s)
    solve_linear_system(factor, matrix, sources; tol = residual_tolerance(cfg, T))
end

function multiple_solve(cfg, s::AMGXSolver, matrix::SparseMatrixCSC{T,V}, sources::Vector{T}) where {T,V}
    handle = construct_amgx_solver(matrix, s)
    try
        solve_linear_system(handle, matrix, sources; tol = residual_tolerance(cfg, T))
    finally
        release!(handle)
    end
end
