struct RasterFlags
    is_raster::Bool
    is_pairwise::Bool
    is_advanced::Bool
    is_onetoall::Bool
    is_alltoone::Bool
    grnd_file_is_res::Bool
    policy::Symbol
    four_neighbors::Bool
    avg_res::Bool
    solver::String
    outputflags::OutputFlags
end

function raster_pairwise(T, V, cfg)::Matrix{T}

    # Get input
    rasterdata = load_raster_data(T, V, cfg)

    # Get compute flags
    flags = get_raster_flags(cfg)

    pt_file_contains_polygons = length(rasterdata.points_rc[1]) !=
                                length(unique(rasterdata.points_rc[3]))

    if pt_file_contains_polygons
        _pt_file_polygons_path(rasterdata, flags, cfg)
    else
        _pt_file_no_polygons_path(rasterdata, flags, cfg)
    end
end

function get_raster_flags(cfg)

    # Computation flags
    is_raster = true
    is_pairwise = cfg["scenario"] in PAIRWISE
    is_advanced = cfg["scenario"] in ADVANCED
    is_onetoall = cfg["scenario"] in ONETOALL
    is_alltoone = cfg["scenario"] in ALLTOONE
    four_neighbors = cfg["connect_four_neighbors_only"] in TRUELIST
    avg_res = cfg["connect_using_avg_resistances"] in TRUELIST
    solver = cfg["solver"]
    ground_file_is_resistances =
        cfg["ground_file_is_resistances"] in TRUELIST
    policy = Symbol(cfg["remove_src_or_gnd"])

    # Output Flags
    o = get_output_flags(cfg)

    RasterFlags(is_raster, is_pairwise, is_advanced,
                is_onetoall, is_alltoone,
                ground_file_is_resistances, policy,
                four_neighbors, avg_res, solver, o)
end

function _pt_file_no_polygons_path(rasterdata::RasData{T,V},
                    flags, cfg)::Matrix{T} where {T,V}

    graphdata = compute_graph_data_no_polygons(rasterdata, flags)
    r = single_ground_all_pairs(graphdata, flags, cfg)

    write_cum_maps(graphdata.cum, rasterdata.cellmap, cfg, rasterdata.hbmeta,
                    flags.outputflags.write_max_cur_maps,
                    flags.outputflags.write_cum_cur_map_only)

    r
end

function _pt_file_polygons_path(rasterdata::RasData{T,V},
                        flags, cfg)::Matrix{T} where {T,V}

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

    # Cumulative maps
    cum = initialize_cum_maps(gmap, flags.outputflags.write_max_cur_maps)

    pts = unique(points_rc[3])
    resistances = -1 * ones(length(pts), length(pts))

    n = calc_num_pairs(pts)
    csinfo("Total number of pair solves = $n")

    k = 1
    for i = 1:size(pts, 1)
        pt1 = pts[i]
        for j = i+1:size(pts, 1)
            pt2 = pts[j]
            csinfo("Solving pair $k of $n")
            k += 1
            graphdata = compute_graph_data_polygons(rasterdata, flags, pt1, pt2, cum)
            pairwise_resistance = single_ground_all_pairs(graphdata, flags, cfg, false)
            resistances[i,j] = resistances[j,i] = pairwise_resistance[2,3]
        end
    end
    for i = 1:size(pts, 1)
        resistances[i,i] = 0
    end
    P = [0, pts...]
    r = hcat(P, vcat(pts', resistances))

    write_cum_maps(cum, gmap, cfg, rasterdata.hbmeta,
                   flags.outputflags.write_max_cur_maps,
                   flags.outputflags.write_cum_cur_map_only)

    # save resistances
    save_resistances(r, cfg)

    r
end

function calc_num_pairs(pts)
    n = 0
    for i = 1:size(pts, 1)
        for j = i+1:size(pts, 1)
            n += 1
        end
    end
    n
end

function compute_graph_data_polygons(rasterdata::RasData{T,V},
                            flags, pt1, pt2, cum)::GraphData{T,V} where {T,V}

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
    # x = find(x -> x == pt1, points_rc[3])[1]
    # y = find(x -> x == pt2, points_rc[3])[1]
    x = something(findfirst(isequal(pt1), points_rc[3]), 0)
    y = something(findfirst(isequal(pt2), points_rc[3]), 0)
    c1 = nodemap[points_rc[1][x], points_rc[2][x]]
    c2 = nodemap[points_rc[1][y], points_rc[2][y]]
    points = V[c1, c2]

    # Exclude pairs array
    exclude_pairs = Tuple{V,V}[]

    GraphData(G, cc, points, [pt1, pt2],
            exclude_pairs, nodemap, newpoly, hbmeta, gmap, cum)
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

function compute_graph_data_no_polygons(data::RasData{T,V},
                    flags)::GraphData{T,V} where {T,V}

    # Data
    cellmap = data.cellmap
    polymap = data.polymap
    points_rc = data.points_rc
    included_pairs = data.included_pairs
    hbmeta = data.hbmeta

    # Flags
    avg_res = flags.avg_res
    four_neighbors = flags.four_neighbors
    write_max_cur_maps = flags.outputflags.write_max_cur_maps

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
        exclude_pairs = Tuple{V,V}[]
    end

    points = zeros(V, length(points_rc[3]))
    for (i,v) in enumerate(zip(points_rc[1], points_rc[2]))
        points[i] = nodemap[v...]
    end

    # Cumulative current maps
    cum = initialize_cum_maps(cellmap, write_max_cur_maps)

    GraphData(G, cc, points, points_rc[3],
                exclude_pairs, nodemap, polymap,
                hbmeta, cellmap, cum)
end
Base.isempty(t::IncludeExcludePairs) = t.mode == :undef

function generate_exclude_pairs(points_rc, included_pairs::IncludeExcludePairs{V}) where V

    exclude_pairs_array = Tuple{V,V}[]
    mat = included_pairs.include_pairs
    mode = included_pairs.mode == :include ? 0 : 1

    prune_points!(points_rc, included_pairs.point_ids)
    for j = 1:size(mat, 2)
        for i = 1:size(mat, 1)
            if mat[i,j] == mode && mat[j,i] == mode
                push!(exclude_pairs_array, (i,j))
            end
        end
    end

    exclude_pairs_array
end

function construct_node_map(gmap, polymap::Matrix{V}) where V

    nodemap = zeros(V, size(gmap))
    ind = gmap .> 0
    nodemap[ind] = 1:sum(ind)

    if isempty(polymap)
        return nodemap
    end

    #=d = Dict{Int, Vector{Int}}()
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
            index::INT = findfirst(x -> x == val, polymap)
            nodemap[i] = nodemap[index]
        end
    end=#

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
