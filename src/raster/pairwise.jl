function raster_pairwise(T, V, cfg)::Matrix{T}

    # Get input
    rasterdata = @timeit CSTIMER[] "load raster data" load_raster_data(T, V, cfg)

    pt_file_contains_polygons = length(rasterdata.points_rc[1]) !=
                                length(unique(rasterdata.points_rc[3]))

    if pt_file_contains_polygons
        _pt_file_polygons_path(rasterdata, cfg)
    else
        _pt_file_no_polygons_path(rasterdata, cfg)
    end
end

function _pt_file_no_polygons_path(rasterdata::RasterData{T,V}, cfg)::Matrix{T} where {T,V}

    graphdata = @timeit CSTIMER[] "construct graph" compute_graph_data_no_polygons(rasterdata, cfg)
    r = @timeit CSTIMER[] "solve pairwise resistances" single_ground_all_pairs(graphdata, cfg)

    if cfg.write_cur_maps || cfg.write_cum_cur_map_only
        @timeit CSTIMER[] "write cumulative current maps" write_cum_maps(graphdata.cum, rasterdata.cellmap, cfg, rasterdata.hbmeta)
    end

    r
end


function _pt_file_polygons_path(rasterdata::RasterData{T,V}, cfg)::Matrix{T} where {T,V}

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
	included_pairs = rasterdata.included_pairs
	exclude_pairs = isempty(included_pairs) ? Vector{Tuple{V,V}}() : 
						generate_exclude_pairs(points_rc, included_pairs)

    # Cumulative maps
    cum = initialize_cum_maps(gmap, cfg.write_max_cur_maps)

    pts = unique(points_rc[3])
    resistances = -1 * ones(length(pts), length(pts))

    pairs = pairs_to_solve(pts, exclude_pairs)
    n = length(pairs)
    @info("Total number of pair solves = $n")

    for (k, (i, j)) in enumerate(pairs)
        pt1, pt2 = pts[i], pts[j]
        @info("Solving pair $k of $n")
        graphdata = compute_graph_data_polygons(rasterdata, pt1, pt2, cum, cfg)
        pairwise_resistance = single_ground_all_pairs(graphdata, cfg, false)
        resistances[i,j] = resistances[j,i] = pairwise_resistance[2,3]
    end
    for i = 1:size(pts, 1)
        resistances[i,i] = 0
    end
    P = [0, pts...]
    r = hcat(P, vcat(pts', resistances))

    if cfg.write_cur_maps || cfg.write_cum_cur_map_only
        write_cum_maps(cum, gmap, cfg, rasterdata.hbmeta)
    end

    # save resistances
    save_resistances(r, cfg)

    r
end

"""
    pairs_to_solve(pts, exclude_pairs)

Index pairs `(i, j)` with `i < j` into `pts` whose node IDs are not listed in
`exclude_pairs` (in either order). This is exactly the set of pairs the
polygon path solves, so counting it keeps the reported number of pair solves
consistent with the include/exclude file (issue #341).
"""
function pairs_to_solve(pts, exclude_pairs)
    excluded = Set(exclude_pairs)
    [(i, j) for i in 1:length(pts) for j in i+1:length(pts)
        if (pts[i], pts[j]) ∉ excluded && (pts[j], pts[i]) ∉ excluded]
end


function compute_graph_data_polygons(rasterdata::RasterData{T,V},
                            pt1, pt2, cum, cfg)::GraphProblem{T,V} where {T,V}


    # Data
    gmap = rasterdata.cellmap
    polymap = rasterdata.polymap
    points_rc = rasterdata.points_rc
    hbmeta = rasterdata.hbmeta

    # Options
    avg_res = cfg.connect_using_avg_resistances
    four_neighbors = cfg.connect_four_neighbors_only

    # Construct new polymap
    newpoly = create_new_polymap(gmap, polymap, points_rc, pt1, pt2)
    nodemap = construct_node_map(gmap, newpoly)

    # Construct graph
    a = construct_graph(gmap, nodemap, avg_res, four_neighbors)
    G = laplacian!(a)

    # Find connected components
    cc = connected_components(a)

    # Construct points vector
    x,y = 0,0
    # x = find(x -> x == pt1, points_rc[3])[1]
    # y = find(x -> x == pt2, points_rc[3])[1]
    x = something(findfirst(isequal(pt1), points_rc[3]), 0)
    y = something(findfirst(isequal(pt2), points_rc[3]), 0)
    c1 = nodemap[points_rc[1][x], points_rc[2][x]]
    c2 = nodemap[points_rc[1][y], points_rc[2][y]]
    points = V[c1, c2]

    # Exclude pairs array
    exclude_pairs = Tuple{V,V}[]

    solver = get_solver(cfg)
    
    GraphProblem(G, cc, points, [pt1, pt2], 
            exclude_pairs, nodemap, newpoly, hbmeta, gmap, cum, solver)
end

function compute_graph_data_no_polygons(data::RasterData{T,V}, cfg)::GraphProblem{T,V} where {T,V}


    # Data
    cellmap = data.cellmap
    polymap = data.polymap
    points_rc = data.points_rc
    included_pairs = data.included_pairs
    hbmeta = data.hbmeta

    # Options
    avg_res = cfg.connect_using_avg_resistances
    four_neighbors = cfg.connect_four_neighbors_only
    write_max_cur_maps = cfg.write_max_cur_maps

    # Nodemap and graph construction
    nodemap = construct_node_map(cellmap, polymap)
    G = construct_graph(cellmap, nodemap, avg_res, four_neighbors)
    G = laplacian!(G)

    # Connected Components
    cc = connected_components(G)

    # Generate exclude pairs array
    if !isempty(included_pairs)
        exclude_pairs = generate_exclude_pairs(points_rc, included_pairs)
    else
        exclude_pairs = Tuple{V,V}[]
    end

    points = zeros(V, length(points_rc[3]))
    for (i,v) in enumerate(zip(points_rc[1], points_rc[2]))
        points[i] = nodemap[v...]
    end

    # Cumulative current maps
    cum = initialize_cum_maps(cellmap, write_max_cur_maps)


    solver = get_solver(cfg)

    GraphProblem(G, cc, points, points_rc[3], 
                exclude_pairs, nodemap, polymap, 
                hbmeta, cellmap, cum, solver)

end

function generate_exclude_pairs(points_rc, included_pairs::IncludeExcludePairs{V}) where V

    # Include mode prunes focal nodes not mentioned in the include file;
    # exclude mode keeps every node and only drops the listed pairs.
    if included_pairs.mode == :include
        prune_points!(points_rc, included_pairs.point_ids)
    end

    exclude_pairs_from(included_pairs)
end

function construct_node_map(gmap, polymap::Matrix{V}) where V

    nodemap = zeros(V, size(gmap))
    ind = gmap .> 0
    nodemap[ind] = 1:sum(ind)

    if isempty(polymap)
        return nodemap
    end

    idx = gmap .> 0 

    polymap_pruned = zeros(V, size(gmap))
    polymap_pruned[idx] = polymap[idx]


    polynums = unique(polymap)
    for i = 1:size(polynums, 1)
        polynum = polynums[i]
        if polynums[i] != 0
            idx1 = findall(x -> x == polynum, polymap_pruned)
            idx2 = findall(x -> x == polynum, polymap)
            if length(idx1) > 0
                nodemap[idx2] .= nodemap[idx1[1]]
            end
        end
    end
    relabel!(nodemap, V(1))

    nodemap
end

function relabel!(nodemap::Matrix{V}, offset = V(0)) where V
    oldlabels = nodemap[findall(x->x!=0,nodemap)]
    newlabels = zeros(V, size(oldlabels))
    s = sort(oldlabels)
    perm = sortperm(oldlabels)
    prepend!(s, s[1] - 1)
    f = findall(x->x!=0,diff(s))
    newlabels[f] .= 1
    newlabels = cumsum(newlabels)
    newlabels[perm] = copy(newlabels)
    nodemap[findall(x->x!=0,nodemap)] = newlabels .- V(1) .+ offset
end

function construct_graph(gmap, nodemap::Matrix{S}, avg_res, four_neighbors) where S
    f1 = avg_res ? res_avg : cond_avg
    f2 = avg_res ? weirder_avg : weird_avg
    I = Vector{S}()
    J = Vector{S}()
    V = Vector{eltype(gmap)}()
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

function create_new_polymap(gmap, polymap::Matrix{V}, points_rc,
                pt1 = 0, pt2 = 0, point_map = Matrix{V}(undef,0,0)) where V

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
            for i in findall(x->x!=0,point_map)
                if polymap[i] == 0
                    newpoly[i] = point_map[i] + k
                end
            end
        else
            k = max(maximum(polymap), maximum(point_map))
            for i in findall(x->x!=0,point_map)
                v1 = point_map[i]
                v2 = newpoly[i]
                if v2 == 0
                    newpoly[i] = k + v1
                    continue
                end
                if v1 != v2
                    ind = findall(x -> x == v2, newpoly)
                    newpoly[ind] .= v1
                end
            end
        end
        return newpoly
    end

    if isempty(polymap)
        newpoly = zeros(V, size(gmap)...)
        id1 = findall(x -> x == pt1, points_rc[3])
        id2 = findall(x -> x == pt2, points_rc[3])
        map(x -> newpoly[f(x)...] = pt1, id1)
        map(x -> newpoly[f(x)...] = pt2, id2)
        return newpoly
    else
        newpoly = deepcopy(polymap)
        k = maximum(polymap)
        for p in (pt1, pt2)
            # find the locations of the point
            idx = findall(x -> x == p, points_rc[3])

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
                    overlap = findall(in(vals), polymap)
                    newpoly[overlap] .= k + 1
                    k += 1
                end
            end
        end
        return newpoly
    end
end
