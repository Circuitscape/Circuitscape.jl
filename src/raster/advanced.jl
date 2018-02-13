struct AdvancedData{T,V}
    G::SparseMatrixCSC{T,V}
    cc::Vector{Vector{V}}
    nodemap::Matrix{V}
    polymap::Matrix{V}
    hbmeta::RasterMeta    
    sources::Vector{T}
    grounds::Vector{T}
    finite_grounds::Vector{T}
    check_node::V
    src::V
end

function raster_advanced(T, cfg)
    
    # Load raster data
    rasterdata = load_raster_data(T, cfg)

    # Get flags
    flags = get_raster_flags(cfg)

    # Generate advanced 
    advanced_data = compute_advanced_data(rasterdata, flags)

    # Send to main kernel
    advanced_kernel(advanced_data, flags, cfg)
end

function compute_advanced_data(data, flags)

    # Data
    cellmap = data.cellmap
    polymap = data.polymap
    points_rc = data.points_rc
    included_pairs = data.included_pairs    
    hbmeta = data.hbmeta
    
    # Flags
    avg_res = flags.avg_res
    four_neighbors = flags.four_neighbors
    policy = flags.policy

    # Nodemap and graph construction
    nodemap = construct_node_map(cellmap, polymap)
    G = construct_graph(cellmap, nodemap, avg_res, four_neighbors)
    G = laplacian(G)

    # Connected Components
    cc = connected_components(SimpleWeightedGraph(G))

    # Advanced mode specific stuff
    sources, grounds, finite_grounds = 
            get_sources_and_grounds(data, flags, G, nodemap)

    AdvancedData(G, cc, nodemap, polymap, hbmeta,
                sources, grounds, finite_grounds, -1, 0)
end

function get_sources_and_grounds(data, flags, G, nodemap)

    # Data 
    source_map = data.source_map
    ground_map = data.ground_map

    # Flags
    is_raster = flags.is_raster
    grnd_file_is_res = flags.grnd_file_is_res
    policy = flags.policy

    # Initialize sources and grounds
    sources = zeros(eltype(G), size(G, 1))
    grounds = zeros(eltype(G), size(G, 1))

    if is_raster
        (i1, j1, v1) = findnz(source_map)
        (i2, j2, v2) = findnz(ground_map)
        for i = 1:size(i1, 1)
            v = Int(nodemap[i1[i], j1[i]])
            if v != 0
                sources[v] += v1[i]
            end
        end
        for i = 1:size(i2, 1)
            v = Int(nodemap[i2[i], j2[i]])
            if v != 0
                grounds[v] += v2[i]
            end
        end
    else
        if grnd_file_is_res
            ground_map[:,2] = 1 ./ ground_map[:,2]
        end
        sources[Int.(source_map[:,1])] = source_map[:,2]
        grounds[Int.(ground_map[:,1])] = ground_map[:,2]
    end
    sources, grounds, finitegrounds = 
        resolve_conflicts(sources, grounds, policy)
end

function resolve_conflicts(sources, grounds, policy)
    
    finitegrounds = similar(sources)
    l = size(sources, 1)

    finitegrounds = map(x -> x < Inf ? x : 0., grounds)
    if count(x -> x != 0, finitegrounds) == 0
        finitegrounds = [-9999.]
    end

    conflicts = falses(l)
    for i = 1:l
        conflicts[i] = sources[i] != 0 && grounds[i] != 0
    end

    if any(conflicts)
        if policy == :rmvsrc
            sources[find(conflicts)] = 0
        elseif policy == :rmvgnd
            grounds[find(conflicts)] = 0
        elseif policy == :rmvall
            sources[find(conflicts)] = 0
        end
    end

    infgrounds = map(x -> x == Inf, grounds)
    infconflicts = map((x,y) -> x > 0 && y > 0, infgrounds, sources)
    grounds[infconflicts] = 0
    
    sources, grounds, finitegrounds
end

function advanced_kernel(data, flags, cfg)

    # Data 
    G = data.G
    nodemap = data.nodemap
    polymap = data.polymap
    hbmeta = data.hbmeta
    sources = data.sources
    grounds = data.grounds
    finitegrounds = data.finite_grounds
    cc = data.cc
    src = data.src
    check_node = data.check_node

    # Flags
    is_raster = flags.is_raster
    is_alltoone = flags.is_alltoone
    is_onetoall = flags.is_onetoall
    write_volt_maps = flags.outputflags.write_volt_maps
    write_cur_maps = flags.outputflags.write_cur_maps

    volt = zeros(eltype(G), size(nodemap))
    ind = find(nodemap)
    f_local = Vector{eltype(G)}()
    solver_called = false
    voltages = Vector{eltype(G)}()
    outvolt = alloc_map(hbmeta)
    outcurr = alloc_map(hbmeta)

    for c in cc

        if check_node != -1 && !(check_node in c)
            continue
        end

        # a_local = laplacian(a[c, c])
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

        voltages = multiple_solver(cfg, a_local, s_local, g_local, f_local)
        solver_called = true

        if write_volt_maps && is_raster
            local_nodemap = construct_local_node_map(nodemap, c, polymap)
            accum_voltages!(outvolt, voltages, local_nodemap, hbmeta)
        end
        if write_cur_maps && is_raster
            local_nodemap = construct_local_node_map(nodemap, c, polymap)
            accum_currents!(outcurr, voltages, cfg, a_local, voltages, f_local, local_nodemap, hbmeta)
        end

        for i in eachindex(volt)
            if i in ind
                val = Int(nodemap[i])
                if val in c
                    idx = findfirst(x -> x == val, c)
                    volt[i] = voltages[idx]
                end
            end
        end
    end

    name = src == 0 ? "" : "_$(Int(src))"
    if write_volt_maps
        if !is_raster
            write_volt_maps(name, voltages, collect(1:size(a,1)), cfg)
        else
            write_aagrid(outvolt, name, cfg, hbmeta, voltage = true)
        end
    end

    if write_cur_maps
        if !is_raster
            write_cur_maps(laplacian(a), voltages, finitegrounds, collect(1:size(a,1)), name, cfg)
        else
            write_aagrid(outcurr, name, cfg, hbmeta)
        end
    end

    if !is_raster
        v = [collect(1:size(G, 1))  voltages]
        return v
    end

    scenario = cfg["scenario"]
    if !solver_called
        return [-1.]
    end

    if is_onetoall
        idx = find(source_map)
        val = volt[idx] ./ source_map[idx]
        if val[1] ≈ 0
            return [-1.]
        else
            return val
        end
    elseif is_alltoone
        return [0.]
    end

    return volt
end

function multiple_solver(cfg, a, sources, grounds, finitegrounds)
    
    asolve = deepcopy(a)
    if finitegrounds[1] != -9999
        asolve = a + spdiagm(finitegrounds, 0, size(a, 1), size(a, 1))
    end

    infgrounds = find(x -> x == Inf, grounds)
    deleteat!(sources, infgrounds)
    dst_del = Int[]
    append!(dst_del, infgrounds)
    r = collect(1:size(a, 1))
    deleteat!(r, dst_del)
    asolve = asolve[r, r]

    M = aspreconditioner(smoothed_aggregation(asolve))
    volt = solve_linear_system(cfg, asolve, sources, M)

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
    