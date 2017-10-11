#=immutable RasterData
    cellmap::Matrix{Float64}
    polymap::Matrix{Float64}
    source_map::Matrix{Float64}
    ground_map::Matrix{Float64}
    points_rc::Tuple{Vector{Int},Vector{Int},Vector{Float64}}
    strengths::Matrix{Float64}
    included_pairs::IncludeExcludePairs
end=#

function compute{S}(obj::Raster{S}, cfg)

    flags = inputflags(obj, cfg)
    data, hbmeta = grab_input(obj, flags)
    compflags = computeflags(obj, cfg, data.points_rc)
    compute(obj, data, compflags, hbmeta, cfg)
end
compute(::Raster{Pairwise}, data, compflags, hbmeta, cfg) = pairwise_module(data,
                                                    compflags, hbmeta, cfg)
function compute(::Raster{Advanced}, rdata, compflags, hbmeta, cfg)

    gmap = rdata.cellmap
    polymap = rdata.polymap
    points_rc = rdata.points_rc
    included_pairs = rdata.included_pairs
    avg_res = compflags.avg_res
    four_neighbors = compflags.four_neighbors
    nodemap = construct_node_map(gmap, polymap)
    a = construct_graph(gmap, nodemap, avg_res, four_neighbors)
    voltages = advanced(cfg, a,rdata.source_map, rdata.ground_map,
                        nodemap = nodemap, hbmeta = hbmeta, polymap = rdata.polymap)
    return voltages
end
function compute{S<:Union{OneToAll,AllToOne}}(::Raster{S}, rdata, compflags, hbmeta, cfg)
    gmap = rdata.cellmap
    polymap = rdata.polymap
    points_rc = rdata.points_rc
    included_pairs = rdata.included_pairs
    onetoall(cfg, gmap, polymap, points_rc;
                            included_pairs = rdata.included_pairs,
                            strengths = rdata.strengths, hbmeta = hbmeta)
end

function pairwise_module(rdata, compflags, hbmeta, cfg)
    ptpoly = compflags.ptpoly
    _pairwise(ptpoly, rdata, compflags, hbmeta, cfg)
end

function _pairwise(::PointFileNoPolygons, rdata, compflags, hbmeta, cfg)

    gmap = rdata.cellmap
    polymap = rdata.polymap
    points_rc = rdata.points_rc
    avg_res = compflags.avg_res
    included_pairs = rdata.included_pairs
    four_neighbors = compflags.four_neighbors

    mode = included_pairs.mode == :include ? 0 : 1
    nodemap = construct_node_map(gmap, polymap)
    a = construct_graph(gmap, nodemap, avg_res, four_neighbors)
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

    resistances = single_ground_all_pair_resistances(a, c, cfg;
                                    exclude = exclude_pairs_array,
                                    nodemap = nodemap,
                                    orig_pts = points_rc[3],
                                    polymap = polymap,
                                    hbmeta = hbmeta)
end
function _pairwise(::PointFileContainsPolygons, rdata, compflags, hbmeta, cfg)

    # get unique list of points
    # for every point pair do
        # construct new polymap
        # construct new nodemap
        # construct new graph
        # solve for two points
    # end
    gmap = rdata.cellmap
    polymap = rdata.polymap
    points_rc = rdata.points_rc
    avg_res = compflags.avg_res
    four_neighbors = compflags.four_neighbors

    pts = unique(points_rc[3])
    resistances = -1 * ones(length(pts), length(pts))

    for i = 1:size(pts, 1)
        pt1 = pts[i]
        for j = i+1:size(pts, 1)
            pt2 = pts[j]
            newpoly = create_new_polymap(gmap, polymap, points_rc, pt1 = pt1, pt2 = pt2)
            nodemap = construct_node_map(gmap, newpoly)
            a = construct_graph(gmap, nodemap, avg_res, four_neighbors)
            x,y = 0,0
            x = find(x -> x == pt1, points_rc[3])[1]
            y = find(x -> x == pt2, points_rc[3])[1]
            c1 = nodemap[points_rc[1][x], points_rc[2][x]]
            c2 = nodemap[points_rc[1][y], points_rc[2][y]]
            c = Int[c1, c2]
            pairwise_resistance = single_ground_all_pair_resistances(a, c, cfg; orig_pts =[points_rc[3][x], points_rc[3][y]],
                                                                                    nodemap = nodemap,
                                                                                    polymap = Polymap(newpoly),
                                                                                    hbmeta = hbmeta)
            resistances[i,j] = resistances[j,i] = pairwise_resistance[1,2]
        end
    end
    for i = 1:size(pts, 1)
        resistances[i,i] = 0
    end
    resistances
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

function create_new_polymap(gmap, p, points_rc; pt1 = 0, pt2 = 0, point_map = Matrix{Float64}(0,0))

    f(x) = (points_rc[1][x], points_rc[2][x])
    polymap = typeof(p) == NoPoly ? Matrix{Float64}(0,0) : p.polymap
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
Base.isempty(p::Polymap) = isempty(p.polymap)

res_avg(x, y) = 1 / ((1/x + 1/y) / 2)
cond_avg(x, y) = (x + y) / 2
weird_avg(x,y) = (x + y) / (2*√2)
weirder_avg(x, y) = 1 / (√2 * (1/x + 1/y) / 2)

function construct_node_map(gmap, ::NoPoly)
    nodemap = zeros(Int, size(gmap))
    ind::Vector{Int64} = find(x -> x > 0, gmap)
    nodemap[ind] = 1:length(ind)
    nodemap
end

construct_node_map(gmap, p::Polymap) = construct_node_map(gmap, p.polymap)
function construct_node_map{T}(gmap, polymap::Matrix{T})

    nodemap = zeros(size(gmap))
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
    avg_res = cfg["connect_using_avg_resistances"] == "True"
    one_to_all = cfg["scenario"] == "one-to-all"

    a = construct_graph(gmap, nodemap, avg_res, four_neighbors)
    cc = connected_components(SimpleWeightedGraph(a))
    debug("There are $(size(a, 1)) points and $(length(cc)) connected components")

    source_map = Matrix{Float64}(0, 0)
    ground_map = Matrix{Float64}(0, 0)
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
            polymap = create_new_polymap(gmap, Polymap(polymap), points_rc, point_map = point_map)
            nodemap = construct_node_map(gmap, polymap)
            a = construct_graph(gmap, nodemap, avg_res, four_neighbors)
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
            v = advanced(cfg, a, source_map, ground_map; nodemap = nodemap, policy = :rmvgnd,
                            check_node = check_node, src = n, polymap = Polymap(newpoly), hbmeta = hbmeta)
        else
            v = advanced(cfg, a, source_map, ground_map; nodemap = nodemap, policy = :rmvsrc,
                            check_node = check_node, src = n, polymap = Polymap(newpoly), hbmeta = hbmeta)
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

function update_voltmatrix!(voltmatrix, local_nodemap, voltages, hbmeta, c, r, j, cc)

    for i = 2:size(c, 1)
        ind = findfirst(cc, c[i])
        if ind != 0
            voltageAtPoint = voltages[ind]
            voltageAtPoint = 1 - (voltageAtPoint/r)
            voltmatrix[i,j] = voltageAtPoint
        end
    end
end

function update_shortcut_resistances!(anchor, shortcut, resistances, voltmatrix, c, cc)
    check = map(x -> x in cc, c)
    l = size(resistances, 1)
    for pointx = 1:l
        if check[pointx]
            R1x = resistances[anchor, pointx]
            if R1x != -1
                shortcut[pointx, anchor] = shortcut[anchor, pointx] = R1x
                for point2 = pointx:l
                    if check[point2]
                        R12 = resistances[anchor, point2]
                        if R12 != -1
                            if R1x != -777
                                shortcut[anchor, point2] = shortcut[point2, anchor] = R12
                                Vx = voltmatrix[pointx, point2]
                                R2x = 2*R12*Vx + R1x - R12
                                if shortcut[point2, pointx] != -777
                                    shortcut[point2, pointx] = shortcut[pointx, point2] = R2x
                                end
                            else
                                shortcut[pointx, :] = shortcut[:, pointx] = -777
                            end
                        end
                    end
                end
            end
        end
    end
end
