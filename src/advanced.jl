# An advanced-mode problem: the graph Laplacian, its connected components,
# the source and ground vectors over the nodes (finite grounds separated out
# for the solver) and where each node lands in the output.
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
    sources_and_grounds(geometry, source_map, ground_map, G, cfg, conflict_policy = policy(cfg))

Source and ground vectors over the graph nodes from the user's source and
ground inputs, with conflicts between the two resolved by `conflict_policy`.
For a raster the inputs are grids read through the nodemap; for a network
they are `(node, value)` lists.
"""
function sources_and_grounds(geometry::Geometry, source_map, ground_map, G, cfg,
                             conflict_policy = policy(cfg))

    # Initialize sources and grounds
    sources = zeros(eltype(G), size(G, 1))
    grounds = zeros(eltype(G), size(G, 1))

    fill_sources_and_grounds!(sources, grounds, source_map, ground_map, cfg, geometry)

    sources, grounds, finitegrounds =
        resolve_conflicts(sources, grounds, conflict_policy)
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
    if cfg.ground_file_is_resistances
        ground_map[:,2] = 1 ./ ground_map[:,2]
    end
    sources[V.(source_map[:,1])] = source_map[:,2]
    grounds[V.(ground_map[:,1])] = ground_map[:,2]
    nothing
end

function resolve_conflicts(sources::Vector{T},
                            grounds::Vector{T}, policy) where T

    l = size(sources, 1)

    finitegrounds = map(x -> x < T(Inf) ? x : T(0.), grounds)
    if count(x -> x != 0, finitegrounds) == 0
        finitegrounds = T.([-9999.])
    end

    conflicts = falses(l)
    for i = 1:l
        conflicts[i] = sources[i] != 0 && grounds[i] != 0
    end

    if any(conflicts)
        if policy == :rmvsrc
            sources[findall(x->x!=0,conflicts)] .= 0
        elseif policy == :rmvgnd
            grounds[findall(x->x!=0,conflicts)] .= 0
        elseif policy == :rmvall
            sources[findall(x->x!=0,conflicts)] .= 0
        end
    end

    infgrounds = map(x -> x == Inf, grounds)
    infconflicts = map((x,y) -> x > 0 && y > 0, infgrounds, sources)
    grounds[infconflicts] .= 0

    sources, grounds, finitegrounds
end

"""
    advanced_kernel(prob, cfg; check_node, name, accumulate_currents)
        -> (voltages, current_map, solver_called)

Solve every connected component that has both a source and a ground, and
write the voltage and current maps the options ask for, suffixed `name`.
`voltages` is over all graph nodes (zero in components not solved); the
current map is the per-cell sum over components for a raster (empty for a
network), filled when `accumulate_currents`. One-to-all / all-to-one solve
one focal node at a time and pass it as `check_node` to skip every other
component.
"""
function advanced_kernel(prob::AdvancedProblem{T,V}, cfg;
                         check_node = V(-1), name = "",
                         accumulate_currents = cfg.write_cur_maps) where {T,V}

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

        a_local = G[c,c]
        s_local = sources[c]
        g_local = grounds[c]

        if sum(s_local) == 0 || sum(g_local) == 0
            continue
        end

        if finitegrounds != [-9999.]
            f_local = finitegrounds[c]
        else
            f_local = finitegrounds
        end

        voltages[c] .+= multiple_solver(cfg, prob.solver, a_local, s_local, g_local, f_local)
        solver_called = true

        local_geometry = restrict(geometry, c)
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

    if write_c_maps || write_cum_cur_map_only
        write_advanced_cur_map(name, voltages, outcurr, G, finitegrounds, geometry, cfg)
    end

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
# network.
write_advanced_volt_map(name, voltages, outvolt, geometry::RasterGeometry, cfg) =
    write_grid(outvolt, name, cfg, geometry.hbmeta, geometry.cellmap, voltage = true)
write_advanced_volt_map(name, voltages, outvolt, geometry::NetworkGeometry, cfg) =
    write_volt_maps(geometry, name, voltages, cfg)

write_advanced_cur_map(name, voltages, outcurr, G, finitegrounds, geometry::RasterGeometry, cfg) =
    write_grid(outcurr, name, cfg, geometry.hbmeta, geometry.cellmap)
write_advanced_cur_map(name, voltages, outcurr, G, finitegrounds, geometry::NetworkGeometry, cfg) =
    write_cur_maps(name, voltages, ComponentData(geometry.nodes, G, geometry), finitegrounds, cfg)


function multiple_solver(cfg, solver, a::SparseMatrixCSC{T,V}, sources, grounds, finitegrounds) where {T,V}

    asolve = deepcopy(a)
    if finitegrounds[1] != -9999
        # asolve = a + spdiagm(finitegrounds, 0, size(a, 1), size(a, 1))
        asolve = a + spdiagm(0 => finitegrounds)
    end

    infgrounds = findall(x -> x == Inf, grounds)
    deleteat!(sources, infgrounds)
    dst_del = V[]
    append!(dst_del, infgrounds)
    r = collect(1:size(a, 1))
    deleteat!(r, dst_del)
    asolve = asolve[r, r]

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
