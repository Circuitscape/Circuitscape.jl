struct AdvancedProblem{T,V,W,Geom<:Geometry}
    G::SparseMatrixCSC{T,V}
    cc::Vector{Vector{V}}
    geometry::Geom
    sources::Vector{T}
    grounds::Vector{T}
    source_map::Matrix{T} # Needed for one to all mode
    finitegrounds::Vector{T}
    check_node::V
    src::V
    solver::W
end

function raster_advanced(T, V, cfg)::Matrix{T}

    # Load raster data
    rasterdata = load_raster_data(T, V, cfg)

    # Generate advanced
    advanced_data = compute_advanced_data(rasterdata, cfg)

    # Send to main kernel
    v, _ = advanced_kernel(advanced_data, cfg)

    v
end


function compute_advanced_data(data::RasterData{T,V}, cfg)::AdvancedProblem{T,V} where {T,V}

    # Data
    cellmap = data.cellmap
    polymap = data.polymap
    points_rc = data.points_rc
    included_pairs = data.included_pairs
    hbmeta = data.hbmeta
    source_map = data.source_map


    # Options
    avg_res = cfg.connect_using_avg_resistances
    four_neighbors = cfg.connect_four_neighbors_only

    # Nodemap and graph construction
    nodemap = construct_node_map(cellmap, polymap)
    A = construct_graph(cellmap, nodemap, avg_res, four_neighbors)
    G = laplacian!(A)

    # Connected Components
    cc = connected_components(G)

    geometry = RasterGeometry(nodemap, polymap, hbmeta, cellmap)

    # Advanced mode specific stuff
    sources, grounds, finitegrounds =
            get_sources_and_grounds(data, cfg, G, geometry)

    solver = get_solver(cfg)

    AdvancedProblem(G, cc, geometry,
                sources, grounds, source_map,
                finitegrounds, V(-1), V(0), solver)

end

function get_sources_and_grounds(data, cfg, G, geometry)

    # Data
    source_map = data.source_map
    ground_map = data.ground_map

    _get_sources_and_grounds(source_map, ground_map, cfg, G, geometry)
end

"""
    _get_sources_and_grounds(source_map, ground_map, cfg, G, geometry, override_policy = :none)

Source and ground vectors over the graph nodes from the user's source and
ground inputs, with conflicts between the two resolved by `policy(cfg)` (or
`override_policy`). For a raster the inputs are grids read through the
nodemap; for a network they are `(node, value)` lists.
"""
function _get_sources_and_grounds(source_map, ground_map,
                                  cfg, G, geometry::Geometry, override_policy = :none)
    conflict_policy = override_policy == :none ? policy(cfg) : override_policy

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

function advanced_kernel(prob::AdvancedProblem{T,V,S}, cfg)::Tuple{Matrix{T},Matrix{T}} where {T,V,S}

    # Data
    G = prob.G
    geometry = prob.geometry
    sources = prob.sources
    grounds = prob.grounds
    finitegrounds = prob.finitegrounds
    cc = prob.cc
    src = prob.src
    check_node = prob.check_node
    source_map = prob.source_map # Need it for one to all mode

    # Options
    alltoone = is_alltoone(cfg)
    onetoall = is_onetoall(cfg)
    write_v_maps = cfg.write_volt_maps
    write_c_maps = cfg.write_cur_maps
    write_cum_cur_map_only = cfg.write_cum_cur_map_only

    f_local = Vector{eltype(G)}()
    solver_called = false
	voltages = zeros(eltype(G), size(G, 1))
    outvolt = alloc_map(geometry)
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
        write_c_maps && accum_currents!(outcurr, a_local, voltages[c], f_local, local_geometry)
    end

    name = src == 0 ? "" : "_$(V(src))"
    write_v_maps && write_advanced_volt_map(name, voltages, outvolt, geometry, cfg)

    # Issue 342: in one-to-all / all-to-one every focal node is active on each
    # solve (one as source, the rest as grounds), so all of them are zeroed.
    if (onetoall || alltoone) && cfg.set_focal_node_currents_to_zero
        active = V[n for n in eachindex(sources) if sources[n] != 0 || grounds[n] != 0]
        zero_focal_nodes!(outcurr, geometry, Set(active))
    end

    if write_c_maps || write_cum_cur_map_only
        write_advanced_cur_map(name, voltages, outcurr, G, finitegrounds, geometry, cfg)
    end

    advanced_result(voltages, outcurr, solver_called, source_map, geometry, cfg)
end

# The value an advanced-mode run returns. A network run returns every node's
# voltage; a raster run returns the voltage grid, or, in one-to-all /
# all-to-one, the effective resistance of the focal node solved.
function advanced_result(voltages::Vector{T}, outcurr, solver_called, source_map,
                         geometry::NetworkGeometry, cfg) where T
    v = [geometry.nodes voltages]
    v, outcurr
end

function advanced_result(voltages::Vector{T}, outcurr, solver_called, source_map,
                         geometry::RasterGeometry, cfg) where T
    if !solver_called
        ret = Matrix{T}(undef,1,1)
        ret[1] = -1
        return ret, outcurr
    end

    if is_onetoall(cfg)
        volt = _create_voltage_map(voltages, geometry.nodemap, geometry.hbmeta)
        idx = findall(x->x!=0,source_map)
        val = volt[idx] ./ source_map[idx]
        if val[1] ≈ 0
            ret = Matrix{T}(undef,1,1)
            ret[1] = -1
            return ret, outcurr
        else
            ret = Matrix{T}(undef,length(val),1)
            ret[:,1] = val
            return ret, outcurr
        end
    elseif is_alltoone(cfg)
        ret = Matrix{T}(undef,1,1)
        ret[1] = 0
        return ret, outcurr
    end

    _create_voltage_map(voltages, geometry.nodemap, geometry.hbmeta), outcurr
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
