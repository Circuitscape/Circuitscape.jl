struct RasterFlags
    is_raster::Bool 
    is_pairwise::Bool
    is_advanced::Bool 
    four_neighbors::Bool
    avg_res::Bool
    solver::String
    outputflags::OutputFlags
end

function raster_pairwise(T, cfg)

    # Get input
    rasterdata = load_raster_data(T, cfg)
    
    # Get compute flags
    flags = get_raster_flags(cfg)
    
    # Compute graph data based on compute flags
    graphdata = compute_graph_data(rasterdata, flags)
    
    # Send to main kernel
    single_ground_all_pairs(graphdata, flags, cfg)
end

function get_raster_flags(cfg)
    
    # Computation flags
    is_raster = true
    is_pairwise = cfg["scenario"] in PAIRWISE
    is_advanced = cfg["scenario"] in ADVANCED
    four_neighbors = cfg["connect_four_neighbors_only"] in truelist
    avg_res = cfg["connect_using_avg_resistances"] in truelist
    solver = cfg["solver"]

    # Output flags
    write_volt_maps = cfg["write_volt_maps"] in truelist
    write_cur_maps = cfg["write_cur_maps"] in truelist
    write_cum_cur_maps_only = cfg["write_cum_cur_map_only"] in truelist
    write_max_cur_maps = cfg["write_max_cur_maps"] in truelist
    set_null_currents_to_nodata = cfg["set_null_currents_to_nodata"] in truelist
    set_null_voltages_to_nodata = cfg["set_null_voltages_to_nodata"] in truelist
    compress_grids = cfg["compress_grids"] in truelist
    log_transform_maps = cfg["log_transform_maps"] in truelist

    o = OutputFlags(write_volt_maps, write_cur_maps,
                    write_cum_cur_maps_only, write_max_cur_maps,
                    set_null_currents_to_nodata, set_null_voltages_to_nodata,
                    compress_grids, log_transform_maps)
    
    RasterFlags(is_raster, is_pairwise, is_advanced, 
                four_neighbors, avg_res, solver, o)
end

function compute_graph_data(rasterdata, flags)
    
    points_rc = rasterdata.points_rc

    pt_file_contains_polygons = 
            length(points_rc[1]) != length(unique(points_rc[3]))

    if !pt_file_contains_polygons
        graphdata = _pairwise_no_polygons(rasterdata, flags)
    end

    graphdata
end

function _pairwise_no_polygons(data, flags)

    # Data
    cellmap = data.cellmap
    polymap = data.polymap
    points_rc = data.points_rc
    included_pairs = data.included_pairs    
    hbmeta = data.hbmeta
    
    # Flags
    avg_res = flags.avg_res
    four_neighbors = flags.four_neighbors

    # Nodemap and graph construction
    nodemap = construct_node_map(cellmap, polymap)
    G = construct_graph(cellmap, nodemap, avg_res, four_neighbors)
    G = laplacian(G)

    # Connected Components
    cc = connected_components(SimpleWeightedGraph(G))    

    # Generate exclude pairs array
    if !isempty(included_pairs)
        exclude_pairs = generate_exclude_pairs(points_rc, included_pairs)
    else
        exclude_pairs = Tuple{Int,Int}[]
    end

    points = zeros(Int, length(points_rc[3]))
    for (i,v) in enumerate(zip(points_rc[1], points_rc[2]))
        points[i] = nodemap[v...]
    end

    GraphData(G, cc, points, points_rc[3], 
                exclude_pairs, nodemap, polymap, hbmeta)
end
Base.isempty(t::IncludeExcludePairs) = t.mode == :undef

function generate_exclude_pairs(points_rc, included_pairs)

    exclude_pairs_array = Tuple{Int,Int}[]
    mat = included_pairs.include_pairs
    mode = included_pairs.mode == :include ? 0 : 1    

        prune_points!(points_rc, included_pairs.point_ids)
        for j = 1:size(mat, 2)
            for i = 1:size(mat, 1)
                if mat[i,j] == mode
                    push!(exclude_pairs_array, (i,j))
                end
            end
        end

    exclude_pairs_array
end

function construct_node_map(gmap, polymap)

    nodemap = zeros(Int, size(gmap))
    
    if isempty(polymap)
        ind = gmap .> 0
        nodemap[ind] = 1:sum(ind)
        return nodemap
    end

    d = Dict{Int, Vector{Int}}()
    #=for i in unique(polymap)
        d[i] = find(x -> x == i, polymap)
    end=#
    m, n = size(polymap)
    I, J, V = findnz(polymap)
    for (i,v) in enumerate(V)
        if v != -9999
            if haskey(d, v)
                push!(d[v], sub2ind((m, n), I[i], J[i]))
            else
                d[v] = [sub2ind((m, n), I[i], J[i])]
            end
        end
    end
    d[0] = find(x -> x == 0, polymap)
    k = 1
    for i in find(gmap)
        if i in d[0]
            nodemap[i] = k
            k += 1
        else
            for key in keys(d)
                if i in d[key]
                    if i == first(d[key])
                        nodemap[i] = k
                        k += 1
                    else
                        f = first(d[key])
                        if polymap[f] != 0 && nodemap[f] == 0
                            nodemap[i] = k
                            nodemap[f] = k
                            k += 1
                        else
                            nodemap[i] = nodemap[f]
                        end
                    end
                end
            end
        end
    end
    for i in eachindex(polymap)
        if polymap[i] != 0 && nodemap[i] == 0
            val::Float64 = polymap[i]
            index::Int64 = findfirst(x -> x == val, polymap)
            nodemap[i] = nodemap[index]
        end
    end

    nodemap
end

function construct_graph(gmap, nodemap, avg_res, four_neighbors)
    f1 = avg_res ? res_avg : cond_avg
    f2 = avg_res ? weirder_avg : weird_avg
    I = Int64[]
    J = Int64[]
    V = Float64[]
    for j = 1:size(gmap, 2)
        for i = 1:size(gmap, 1)
            if nodemap[i,j] == 0
                continue
            else
                # Horizontal neighbour
                if j != size(gmap, 2) && nodemap[i,j+1] != 0
                    push!(I, nodemap[i,j])
                    push!(J, nodemap[i,j+1])
                    push!(V, f1(gmap[i,j], gmap[i,j+1]))
                end

                # Vertical neighbour
                if i != size(gmap, 1) && nodemap[i+1, j] != 0
                    push!(I, nodemap[i,j])
                    push!(J, nodemap[i+1,j])
                    push!(V, f1(gmap[i,j], gmap[i+1,j]))
                end

                if !four_neighbors
                    # Diagonal neighbour
                    if i != size(gmap, 1) && j != size(gmap, 2) && nodemap[i+1, j+1] != 0
                        push!(I, nodemap[i,j])
                        push!(J, nodemap[i+1,j+1])
                        push!(V, f2(gmap[i,j], gmap[i+1,j+1]))
                    end

                    if i != 1 && j != size(gmap, 2) && nodemap[i-1, j+1] != 0
                        push!(I, nodemap[i,j])
                        push!(J, nodemap[i-1,j+1])
                        push!(V, f2(gmap[i,j], gmap[i-1,j+1]))
                    end
                end
            end
        end
    end
    m = maximum(nodemap)
    a = sparse(I,J,V, m, m)
    a = a + a'
    a
end

res_avg(x, y) = 1 / ((1/x + 1/y) / 2)
cond_avg(x, y) = (x + y) / 2
weird_avg(x,y) = (x + y) / (2*√2)
weirder_avg(x, y) = 1 / (√2 * (1/x + 1/y) / 2)