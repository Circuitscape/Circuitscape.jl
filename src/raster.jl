immutable RasterData
    cellmap::Array{Float64,2}
    polymap::Array{Float64,2}
    source_map::Array{Float64,2}
    ground_map::Array{Float64,2}
    points_rc::Tuple{Vector{Int},Vector{Int},Vector{Float64}}
end

function compute_raster(cfg::Inifile)

    # Read inputs
    # gmap, polymap, points_rc = load_maps(cfg)
    rdata = load_maps(cfg)
    gmap = rdata.cellmap
    polymap = rdata.polymap
    points_rc = rdata.points_rc
    c = count(x -> x > 0, gmap)
    info("Resistance/Conductance map has $c nodes")
    four_neighbors = get(cfg, "Connection scheme for raster habitat data",
                                "connect_four_neighbors_only") == "True"
    average_resistances = get(cfg, "Connection scheme for raster habitat data",
                                "connect_using_avg_resistances") == "True"

    scenario = get(cfg, "Circuitscape mode", "scenario")
    if scenario == "pairwise"
        resistances = pairwise_module(gmap, polymap, points_rc, four_neighbors, 
                                        average_resistances)
    elseif scenario == "advanced"
        nodemap = construct_node_map(gmap, polymap)
        a,g = construct_graph(gmap, nodemap, average_resistances, four_neighbors)
        cc = connected_components(g)
        debug("There are $(size(a, 1)) points and $(length(cc)) connected components")
        voltages = advanced(cfg, a, g, rdata.source_map, rdata.ground_map, cc, nodemap = nodemap)
    else
        voltages = onetoall(cfg, gmap, polymap, points_rc)
    end
    #gmap, polymap, points_rc
end

function load_maps(cfg::Inifile)

    # Read raster map
    info("Reading Maps")
    habitat_file = get(cfg, "Habitat raster or graph", "habitat_file")
    is_res = get(cfg, "Habitat raster or graph", "habitat_map_is_resistances") == "True"

    cellmap, habitatmeta = read_cell_map(habitat_file, is_res)

    # Read polygon map
    use_polygons = get(cfg, "Short circuit regions (aka polygons)", "use_polygons") == "True"
    polymap_file = get(cfg, "Short circuit regions (aka polygons)", "polygon_file")
    polymap = use_polygons ? read_polymap(polymap_file, habitatmeta) : Array{Float64,2}()

    scenario = get(cfg, "Circuitscape mode", "scenario")
    point_file = get(cfg, "Options for pairwise and one-to-all and all-to-one modes",
                            "point_file")

    # Default source and ground maps
    source_map = Array{Float64,2}()
    ground_map = Array{Float64,2}()

    points_rc = (Vector{Int}(), Vector{Int}(), Vector{Float64}())
    if scenario == "advanced"
        source_file = get(cfg, "Options for advanced mode", "source_file")
        ground_file = get(cfg, "Options for advanced mode", "ground_file")
        is_res = get(cfg, "Options for advanced mode", "ground_file_is_resistances") == "True"
        source_map, ground_map = read_source_and_ground_maps(source_file, ground_file, habitatmeta, is_res)
    else
        points_rc = read_point_map(point_file, habitatmeta)
    end

    RasterData(cellmap, polymap, source_map, ground_map, points_rc)
end

function pairwise_module(gmap, polymap, points_rc, four_neighbors, average_resistances)

    point_file_contains_polygons = length(points_rc[1]) != length(unique(points_rc[3]))

    if !point_file_contains_polygons
        nodemap = construct_node_map(gmap, polymap)
        a, g = construct_graph(gmap, nodemap, average_resistances, four_neighbors)

        c = zeros(Int, length(points_rc[3]))
        for (i,v) in enumerate(zip(points_rc[1], points_rc[2]))
            c[i] = nodemap[v...]
        end

        resistances = single_ground_all_pair_resistances(a, g, c)
        return resistances
    else
        # get unique list of points
        # for every point pair do
            # construct new polymap
            # construct new nodemap
            # construct new graph
            # solve for two points
        # end

        pts = unique(points_rc[3])
        resistances = -1 * ones(length(pts), length(pts))

        for i = 1:size(pts, 1)
            pt1 = pts[i]
            for j = i+1:size(pts, 1)
                pt2 = pts[j]
                newpoly = create_new_polymap(gmap, polymap, points_rc, pt1, pt2)

                nodemap = construct_node_map(gmap, newpoly)
                a, g = construct_graph(gmap, nodemap, average_resistances, four_neighbors)
                x,y = 0,0
                x = find(x -> x == pt1, points_rc[3])[1]
                y = find(x -> x == pt2, points_rc[3])[1]
                c = Int[nodemap[points_rc[1][x], points_rc[2][x]], nodemap[points_rc[1][y], points_rc[2][y]]]
                pairwise_resistance = single_ground_all_pair_resistances(a, g, c)
                resistances[i,j] = resistances[j,i] = pairwise_resistance[1,2]
            end
        end
        for i = 1:size(pts, 1)
            resistances[i,i] = 0
        end
        return resistances
    end
    return nothing
end

function construct_graph(gmap, nodemap, average_resistances, four_neighbors)
    f1 = average_resistances ? res_avg : cond_avg
    f2 = average_resistances ? weirder_avg : weird_avg
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
    g = Graph(a)
    a, g
end

function create_new_polymap(gmap, polymap, points_rc, pt1, pt2)

    f(x) = (points_rc[1][x], points_rc[2][x])

    if isempty(polymap)
        newpoly = zeros(size(gmap)...)
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

res_avg(x, y) = 1 / ((1/x + 1/y) / 2)
cond_avg(x, y) = (x + y) / 2
weird_avg(x,y) = (x + y) / (2*√2)
weirder_avg(x, y) = 1 / (√2 * (1/x + 1/y) / 2)

function construct_node_map(gmap, polymap)
    nodemap = zeros(size(gmap)) 
    if isempty(polymap)
         ind = find(x -> x > 0, gmap)
         nodemap[ind] = 1:length(ind)
    else
        d = Dict{Int, Vector{Int}}()
        for i in unique(polymap)
            d[i] = find(x -> x == i, polymap)
        end
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
                            if polymap[first(d[key])] != 0 && nodemap[first(d[key])] == 0
                                nodemap[i] = k
                                nodemap[first(d[key])] = k
                                k += 1
                            else
                                nodemap[i] = nodemap[first(d[key])]
                            end
                        end
                    end
                end
            end
        end
    end
    for i in eachindex(polymap)
        if polymap[i] != 0 && nodemap[i] == 0
            val = polymap[i] 
            ind = findfirst(x -> x == val, polymap)
            nodemap[i] = nodemap[ind]
        end
    end

    nodemap 
end

function onetoall(cfg, gmap, polymap, points_rc)

    # Construct point map
    point_map = zeros(size(gmap))
    f(i, x) = points_rc[i][x]
    for x = 1:size(points_rc[1], 1)
        point_map[f(1,x), f(2,x)] = f(3, x)
    end

    points_unique = unique(points_rc[3])

    # Combine polymap and pointmap
    newpoly = deepcopy(polymap)
    if isempty(polymap)
        newpoly = point_map    
    else
        k = maximum(polymap)
        for i in find(point_map)
            if polymap[i] == 0
                newpoly[i] = point_map[i] + k
            end
        end
    end

    nodemap = construct_node_map(gmap, newpoly)

    four_neighbors = get(cfg, "Connection scheme for raster habitat data",
                                "connect_four_neighbors_only") == "True"
    average_resistances = get(cfg, "Connection scheme for raster habitat data",
                                "connect_using_avg_resistances") == "True"
    one_to_all = get(cfg, "Circuitscape mode",
                                "scenario") == "one-to-all"

    a, g = construct_graph(gmap, nodemap, average_resistances, four_neighbors)
    cc = connected_components(g)
    debug("There are $(size(a, 1)) points and $(length(cc)) connected components")

    source_map = Array{Float64,2}()
    ground_map = Array{Float64,2}()
    sources = zeros(size(point_map))
    z = deepcopy(sources)
    res = zeros(size(points_unique, 1))
    num_points_to_solve = size(points_unique, 1)
    for i = 1:num_points_to_solve
        debug("Solving point $i of $num_points_to_solve")
        copy!(sources, z)
        n = points_unique[i]
        if one_to_all
            source_map = map(x -> x == n ? 1 : 0, point_map)
            ground_map = map(x -> x == n ? 0 : x, point_map)
            map!(x -> x > 0 ? Inf : x, ground_map)
        else
            source_map = map(x -> x == n ? 0 : 1, point_map)
            ground_map = map(x -> x == n ? Inf : 0, point_map)
        end

        check_node = nodemap[points_rc[1][i], points_rc[2][i]]
        v = advanced(cfg, a, g, source_map, ground_map, cc; nodemap = nodemap, policy = :rmvgnd, 
                        check_node = check_node)
        res[i] = v[1]
    end
    res
end
