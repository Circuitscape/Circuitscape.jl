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

    pt_file_contains_polygons = length(rasterdata.points_rc[1]) != 
                                length(unique(rasterdata.points_rc[3]))
    
    if pt_file_contains_polygons
        _pt_file_polygons_path(rasterdata, flags, cfg)
    else
        _pt_file_no_polygons_path(rasterdata, flags, cfg)
    end
    # Compute graph data based on compute flags
    # graphdata = compute_graph_data(rasterdata, flags)
    
    # Send to main kernel
    # single_ground_all_pairs(graphdata, flags, cfg)
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

function _pt_file_no_polygons_path(rasterdata, flags, cfg)

    graphdata = compute_graph_data_no_polygons(rasterdata, flags)
    single_ground_all_pairs(graphdata, flags, cfg)
end

function _pt_file_polygons_path(rasterdata, flags, cfg)

    # get unique list of points
    # for every point pair do
        # construct new polymap
        # construct new nodemap
        # construct new graph
        # solve for two points
    # end

    # Data
    gmap = rasterdata.cellmap
    polymap = rasterdata.polymap
    points_rc = rasterdata.points_rc
    avg_res = flags.avg_res
    four_neighbors = flags.four_neighbors

    pts = unique(points_rc[3])
    resistances = -1 * ones(length(pts), length(pts))

    for i = 1:size(pts, 1)
        pt1 = pts[i]
        for j = i+1:size(pts, 1)
            pt2 = pts[j]
            #=newpoly = create_new_polymap(gmap, polymap, points_rc, pt1, pt2)
            nodemap = construct_node_map(gmap, newpoly)
            a = construct_graph(gmap, nodemap, avg_res, four_neighbors)
            x,y = 0,0
            x = find(x -> x == pt1, points_rc[3])[1]
            y = find(x -> x == pt2, points_rc[3])[1]
            c1 = nodemap[points_rc[1][x], points_rc[2][x]]
            c2 = nodemap[points_rc[1][y], points_rc[2][y]]
            c = Int[c1, c2]=#
            graphdata = compute_graph_data_polygons(rasterdata, flags, pt1, pt2)
            # pairwise_resistance = single_ground_all_pair_resistances(a, c, cfg; orig_pts =[points_rc[3][x], points_rc[3][y]],
            #                                                                        nodemap = nodemap,
            #                                                                        polymap = Polymap(newpoly),
            #                                                                        hbmeta = hbmeta)
            pairwise_resistance = single_ground_all_pairs(graphdata, flags, cfg)
            resistances[i,j] = resistances[j,i] = pairwise_resistance[2,3]
        end
    end
    for i = 1:size(pts, 1)
        resistances[i,i] = 0
    end
    resistances
end

function compute_graph_data_polygons(rasterdata, flags, pt1, pt2)

    # Data
    gmap = rasterdata.cellmap
    polymap = rasterdata.polymap
    points_rc = rasterdata.points_rc
    hbmeta = rasterdata.hbmeta

    # Flags
    avg_res = flags.avg_res
    four_neighbors = flags.four_neighbors

    # Construct new polymap
    newpoly = create_new_polymap(gmap, polymap, points_rc, pt1, pt2)
    nodemap = construct_node_map(gmap, newpoly)

    # Construct graph
    a = construct_graph(gmap, nodemap, avg_res, four_neighbors)
    G = laplacian(a)

    # Find connected components
    cc = connected_components(SimpleWeightedGraph(a))

    # Construct points vector
    x,y = 0,0
    x = find(x -> x == pt1, points_rc[3])[1]
    y = find(x -> x == pt2, points_rc[3])[1]
    c1 = nodemap[points_rc[1][x], points_rc[2][x]]
    c2 = nodemap[points_rc[1][y], points_rc[2][y]]
    points = Int[c1, c2]

    # Exclude pairs array
    exclude_pairs = Tuple{Int,Int}[]
    
    GraphData(G, cc, points, [pt1, pt2], 
            exclude_pairs, nodemap, newpoly, hbmeta)
end

#=function compute_graph_data(rasterdata, flags)
    
    points_rc = rasterdata.points_rc

    pt_file_contains_polygons = 
            length(points_rc[1]) != length(unique(points_rc[3]))

    if !pt_file_contains_polygons
        graphdata = _pairwise_no_polygons(rasterdata, flags)
    end

    graphdata
end=#

function compute_graph_data_no_polygons(data, flags)

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

function create_new_polymap(gmap, polymap, points_rc, 
                pt1 = 0, pt2 = 0, point_map = Matrix{Int64}(0,0))
    
    f(x) = (points_rc[1][x], points_rc[2][x])

    if !isempty(point_map)
        # Combine polymap and pointmap
        newpoly = deepcopy(polymap)
        point_file_no_polygons = length(points_rc[3]) == 
                        length(unique(points_rc[3]))
        if isempty(polymap)
            newpoly = point_map
        elseif point_file_no_polygons
            k = maximum(polymap)
            for i in find(point_map)
                if polymap[i] == 0
                    newpoly[i] = point_map[i] + k
                end
            end
        else
            k = max(maximum(polymap), maximum(point_map))
            for i in find(point_map)
                v1 = point_map[i]
                v2 = newpoly[i]
                if v2 == 0
                    newpoly[i] = k + v1
                    continue
                end
                if v1 != v2
                    ind = find(x -> x == v2, newpoly)
                    newpoly[ind] = v1
                end
            end
        end
        return newpoly
    end

    if isempty(polymap)
        newpoly = zeros(Int, size(gmap)...)
        id1 = find(x -> x == pt1, points_rc[3])
        id2 = find(x -> x == pt2, points_rc[3])
        map(x -> newpoly[f(x)...] = pt1, id1)
        map(x -> newpoly[f(x)...] = pt2, id2)
        return newpoly
    else
        newpoly = deepcopy(polymap)
        k = maximum(polymap)
        for p in (pt1, pt2)
            # find the locations of the point
            idx = find(x -> x == p, points_rc[3])

            if length(idx) == 1
                continue
            end
            allzero = mapreduce(x -> polymap[f(x)...] == 0, &, idx)
            if allzero
                map(x -> newpoly[f(x)...] = k + 1, idx)
                k += 1
            else
                nz = filter(x -> polymap[f(x)...]!= 0, idx)
                if length(nz) == 1
                    map(x -> newpoly[f(x)...] = polymap[overlap[1]], idx)
                else
                    coords = map(x -> f(x), nz)
                    vals = map(x -> polymap[x...], coords)
                    overlap = findin(polymap, vals)
                    newpoly[overlap] = k + 1
                    k += 1
                end
            end
        end
        return newpoly
    end
end