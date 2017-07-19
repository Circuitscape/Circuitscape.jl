immutable RasterData
    cellmap::Matrix{Float64}
    polymap::Matrix{Float64}
    source_map::Matrix{Float64}
    ground_map::Matrix{Float64}
    points_rc::Tuple{Vector{Int},Vector{Int},Vector{Float64}}
    strengths::Matrix{Float64}
    included_pairs::IncludeExcludePairs
end

function compute_raster(cfg)

    # Get all variables
    four_neighbors = cfg["connect_four_neighbors_only"] == "True"
    average_resistances = cfg["connect_using_avg_resistances"] == "True"
    scenario = cfg["scenario"]

    # Read inputs
    rdata, hbmeta = load_maps(cfg)
    gmap = rdata.cellmap
    polymap = rdata.polymap
    points_rc = rdata.points_rc
    included_pairs = rdata.included_pairs
    c = count(x -> x > 0, gmap)
    info("Resistance/Conductance map has $c nodes")

    if scenario == "pairwise"
        resistances = pairwise_module(gmap, polymap, points_rc, four_neighbors, 
                                        average_resistances, included_pairs, cfg, hbmeta)
    elseif scenario == "advanced"
        nodemap = construct_node_map(gmap, polymap)
        a,g = construct_graph(gmap, nodemap, average_resistances, four_neighbors)
        cc = connected_components(g)
        debug("There are $(size(a, 1)) points and $(length(cc)) connected components")
        voltages = advanced(cfg, a, g, rdata.source_map, rdata.ground_map, cc, 
                                                                    nodemap = nodemap, hbmeta = hbmeta, polymap = rdata.polymap)
        return voltages
    else
        voltages = onetoall(cfg, gmap, polymap, points_rc; 
                                    included_pairs = rdata.included_pairs,
                                    strengths = rdata.strengths, hbmeta = hbmeta)
        return voltages
    end
end

function load_maps(cfg)

    # Read all user variables
    habitat_file = cfg["habitat_file"]
    is_res = cfg["habitat_map_is_resistances"] == "True"
    use_polygons = cfg["use_polygons"] == "True"
    polymap_file = cfg["polygon_file"]
    scenario = cfg["scenario"]
    point_file = cfg["point_file"]
    mask_file = cfg["mask_file"]
    use_mask = cfg["use_mask"] == "True"
    use_variable_source_strengths = cfg["use_variable_source_strengths"] == "True"
    variable_source_file = cfg["variable_source_file"]
    use_included_pairs = cfg["use_included_pairs"] == "True"
    included_pairs_file = cfg["included_pairs_file"]
    source_file = cfg["source_file"]
    ground_file = cfg["ground_file"]
    is_gres = cfg["ground_file_is_resistances"] == "True"

    info("Reading Maps")

    # Read raster map
    cellmap, habitatmeta = read_cell_map(habitat_file, is_res)

    # Read polygon map
    polymap = use_polygons ? read_polymap(polymap_file, habitatmeta) : Array{Float64,2}()

    if use_mask
        mask = read_polymap(mask_file, habitatmeta)
        map!(x -> x > 0 ? 1 : 0, mask, mask)
        cellmap = cellmap .* mask
        if sum(cellmap) == 0
            throw(ErrorException("Mask file masks everything!"))
        end
    end

    # Default source and ground maps
    source_map = Array{Float64,2}()
    ground_map = Array{Float64,2}()
    points_rc = (Vector{Int}(), Vector{Int}(), Vector{Float64}())
    strengths = Array{Float64,2}()
    included_pairs = IncludeExcludePairs()

    if use_included_pairs
        included_pairs = read_included_pairs(included_pairs_file)
    end

    if use_variable_source_strengths
        strengths = read_point_strengths(variable_source_file, false)
    end

    if scenario == "advanced"
        source_map, ground_map = read_source_and_ground_maps(source_file, ground_file, habitatmeta, is_gres)
    else
        points_rc = read_point_map(point_file, habitatmeta)
    end

    RasterData(cellmap, polymap, source_map, ground_map, points_rc, strengths, included_pairs), habitatmeta
end

function pairwise_module(gmap, polymap, points_rc, four_neighbors, average_resistances, included_pairs, cfg, hbmeta)

    point_file_contains_polygons = length(points_rc[1]) != length(unique(points_rc[3]))
    mode = included_pairs.mode == :include ? 0 : 1

    if !point_file_contains_polygons
        nodemap = construct_node_map(gmap, polymap)
        a, g = construct_graph(gmap, nodemap, average_resistances, four_neighbors)
        exclude_pairs_array = Tuple{Int,Int}[]
        mat = included_pairs.include_pairs

        if !isempty(included_pairs)
            prune_points!(points_rc, included_pairs.point_ids)
            for j = 1:size(mat, 2)
                for i = 1:size(mat, 1)
                    if mat[i,j] == mode
                        push!(exclude_pairs_array, (i,j))
                    end
                end
            end
        end
        c = zeros(Int, length(points_rc[3]))
        for (i,v) in enumerate(zip(points_rc[1], points_rc[2]))
            c[i] = nodemap[v...]
        end

        resistances = single_ground_all_pair_resistances(a, g, c, cfg; 
                                        exclude = exclude_pairs_array, 
                                        nodemap = nodemap, 
                                        orig_pts = points_rc[3],
                                        polymap = polymap,
                                        hbmeta = hbmeta)
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
                newpoly = create_new_polymap(gmap, polymap, points_rc, pt1 = pt1, pt2 = pt2)

                nodemap = construct_node_map(gmap, newpoly)
                a, g = construct_graph(gmap, nodemap, average_resistances, four_neighbors)
                x,y = 0,0
                x = find(x -> x == pt1, points_rc[3])[1]
                y = find(x -> x == pt2, points_rc[3])[1]
                c1 = nodemap[points_rc[1][x], points_rc[2][x]]
                c2 = nodemap[points_rc[1][y], points_rc[2][y]]
                c = Int[c1, c2]
                pairwise_resistance = single_ground_all_pair_resistances(a, g, c, cfg; orig_pts =[points_rc[3][x], points_rc[3][y]],
                                                                                        nodemap = nodemap,
                                                                                        polymap = newpoly,
                                                                                        hbmeta = hbmeta)
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

function create_new_polymap(gmap, polymap, points_rc; pt1 = 0, pt2 = 0, point_map = Array{Float64,2}())

    f(x) = (points_rc[1][x], points_rc[2][x])
    
    if !isempty(point_map)

       # Combine polymap and pointmap
       newpoly = deepcopy(polymap)
       #= if isempty(polymap)
            newpoly = point_map    
        else
            k = maximum(polymap)
            for i in find(point_map)
                if polymap[i] == 0
                    newpoly[i] = point_map[i] + k
                end
            end
        end =#
        point_file_no_polygons = length(points_rc[3]) == length(unique(points_rc[3]))
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
         ind::Vector{Int64} = find(x -> x > 0, gmap)
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
            val::Float64 = polymap[i] 
            index::Int64 = findfirst(x -> x == val, polymap)
            nodemap[i] = nodemap[index]
        end
    end

    nodemap 
end

function onetoall(cfg, gmap, polymap, points_rc; included_pairs = IncludeExcludePairs(), strengths = Matrix{Float64}(), hbmeta = RasterMeta())

    use_variable_strengths = !isempty(strengths)
    use_included_pairs = !isempty(included_pairs)
    mode = included_pairs.mode == :include ? 0 : 1
    if use_included_pairs
        points_unique = included_pairs.point_ids
        prune_points!(points_rc, included_pairs.point_ids)
        if use_variable_strengths
            prune_strengths!(strengths, included_pairs.point_ids)
        end
    end

    # Construct point map
    point_map = zeros(size(gmap))
    f(i, x) = points_rc[i][x]
    for x = 1:size(points_rc[1], 1)
        point_map[f(1,x), f(2,x)] = f(3, x)
    end

    points_unique = unique(points_rc[3])

    newpoly = create_new_polymap(gmap, polymap, points_rc, point_map = point_map)

    nodemap = construct_node_map(gmap, newpoly)

    four_neighbors = cfg["connect_four_neighbors_only"] == "True"
    average_resistances = cfg["connect_using_avg_resistances"] == "True"
    one_to_all = cfg["scenario"] == "one-to-all"

    a, g = construct_graph(gmap, nodemap, average_resistances, four_neighbors)
    cc = connected_components(g)
    debug("There are $(size(a, 1)) points and $(length(cc)) connected components")

    source_map = Array{Float64,2}()
    ground_map = Array{Float64,2}()
    sources = zeros(size(point_map))
    z = deepcopy(sources)
    
    point_ids = included_pairs.point_ids
    res = zeros(size(points_unique, 1))
    num_points_to_solve = size(points_unique, 1)
    original_point_map = copy(point_map)
    unique_point_map = zeros(gmap)
    for i in points_unique
        ind = findfirst(x -> x == i, points_rc[3])
        unique_point_map[f(1,ind), f(2,ind)] = f(3,ind)
    end

    for i = 1:num_points_to_solve
        copy!(point_map, original_point_map)
        str = use_variable_strengths ? strengths[i,2] : 1
        debug("Solving point $i of $num_points_to_solve")
        copy!(sources, z)
        n = points_unique[i]
        if use_included_pairs
            for j = 1:num_points_to_solve
                if i != j && included_pairs.include_pairs[i,j] == mode
                    exclude = point_ids[j]
                    map!(x -> x == exclude ? 0 : x, point_map, point_map)
                end
            end
            polymap = create_new_polymap(gmap, polymap, points_rc, point_map = point_map)
            nodemap = construct_node_map(gmap, polymap)
            a, g = construct_graph(gmap, nodemap, average_resistances, four_neighbors)
        end
        if one_to_all
            #source_map = map(x -> x == n ? str : 0, point_map)
            source_map = map(x -> x == n ? str : 0, unique_point_map)
            ground_map = map(x -> x == n ? 0 : x, point_map)
            map!(x -> x > 0 ? Inf : x, ground_map, ground_map)
        else
            source_map = map(x -> x != 0 ? x : 0, point_map)
            map!(x -> x == n ? 0 : x, source_map, source_map)
            map!(x -> x != 0 ? 1 : x, source_map, source_map)
            ground_map = map(x -> x == n ? Inf : 0, point_map)
        end

        check_node = nodemap[points_rc[1][i], points_rc[2][i]]
        
        if one_to_all
            v = advanced(cfg, a, g, source_map, ground_map, cc; nodemap = nodemap, policy = :rmvgnd, 
                            check_node = check_node, src = n, polymap = newpoly, hbmeta = hbmeta)
        else
            v = advanced(cfg, a, g, source_map, ground_map, cc; nodemap = nodemap, policy = :rmvsrc, 
                            check_node = check_node, src = n, polymap = newpoly, hbmeta = hbmeta)
        end
        res[i] = v[1]
    end
    res
end

Base.isempty(t::IncludeExcludePairs) = t.mode == :undef

function prune_points!(points_rc, point_ids)
    rmv = Int64[]
    for (i,p) in enumerate(points_rc[3])
        if p in point_ids
            continue
        else
            #for it in 1:3 deleteat!(points_rc[it], i) end
            push!(rmv, i)
        end
    end
    for i in 1:3 deleteat!(points_rc[i], rmv) end
end

function prune_strengths!(strengths, point_ids)
    pts = strengths[:,1]
    l = length(pts)
    rmv = Int[]
    for (i,p) in enumerate(pts)
        if !(p in point_ids)
           push!(rmv, i) 
       end
    end
    rng = collect(1:l)
    deleteat!(rng, rmv)
    strengths[rng,:]
end
