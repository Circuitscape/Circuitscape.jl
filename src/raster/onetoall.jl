function raster_one_to_all(T, V, cfg)::Matrix{T}

    # Load the data
    rasterdata = load_raster_data(T, V, cfg)

    # Get flags
    flags = get_raster_flags(cfg)

    # Send to main kernel
    onetoall_kernel(rasterdata, flags, cfg)
end

function onetoall_kernel(data::RasData{T,V}, flags, cfg)::Matrix{T} where {T,V}

    # Data
    strengths = data.strengths
    included_pairs = data.included_pairs
    points_rc = data.points_rc
    gmap = data.cellmap
    polymap = data.polymap
    hbmeta = data.hbmeta
    source_map = data.source_map

    # Flags
    use_variable_strengths = !isempty(strengths)
    use_included_pairs = !isempty(included_pairs)
    mode = included_pairs.mode == :include ? 0 : 1
    one_to_all = flags.is_onetoall 
    avg_res = flags.avg_res
    four_neighbors = flags.four_neighbors   

    if use_included_pairs
        points_unique = included_pairs.point_ids
        prune_points!(points_rc, included_pairs.point_ids)
        if use_variable_strengths
            prune_strengths!(strengths, included_pairs.point_ids)
        end
    end

    # Construct point map
    point_map = zeros(V, size(gmap))
    f(i, x) = points_rc[i][x]
    for x = 1:size(points_rc[1], 1)
        point_map[f(1,x), f(2,x)] = f(3, x)
    end

    points_unique = unique(points_rc[3])

    newpoly = create_new_polymap(gmap, polymap, points_rc, 0, 0, point_map)

    nodemap = construct_node_map(gmap, newpoly)

    a = construct_graph(gmap, nodemap, avg_res, four_neighbors)
    cc = connected_components(SimpleWeightedGraph(a))
    G = laplacian(a)
    csinfo("There are $(size(a, 1)) points and $(length(cc)) connected components")

    # source_map = Matrix{eltype(a)}(0, 0)
    # ground_map = Matrix{eltype(a)}(0, 0)
    s = zeros(eltype(a), size(point_map))
    z = deepcopy(s)
    cum = initialize_cum_maps(gmap, flags.outputflags.write_max_cur_maps)

    point_ids = included_pairs.point_ids
    res = zeros(eltype(a), size(points_unique, 1)) |> SharedArray
    num_points_to_solve = size(points_unique, 1)
    original_point_map = copy(point_map)
    unique_point_map = zeros(V, size(gmap))

    for i in points_unique
        ind = findfirst(x -> x == i, points_rc[3])
        unique_point_map[f(1,ind), f(2,ind)] = f(3,ind)
    end

    # @distributed for i = 1:num_points_to_solve
    function f(i)
        # copyto!(point_map, original_point_map)
        point_map = copy(original_point_map)
        str = use_variable_strengths ? strengths[i,2] : 1
        csinfo("Solving point $i of $num_points_to_solve")
        # copyto!(s, z)
        s = copy(z)
        n = points_unique[i]
        if use_included_pairs
            for j = 1:size(point_ids,1)
                if i != j && included_pairs.include_pairs[i,j] == mode
                    exclude = point_ids[j]
                    map!(x -> x == exclude ? 0 : x, point_map, point_map)
                end
            end
            # polymap = create_new_polymap(gmap, Polymap(polymap), points_rc, point_map = point_map)
            newpoly = create_new_polymap(gmap, polymap, points_rc, 0, 0, point_map)            
            nodemap = construct_node_map(gmap, polymap)
            a = construct_graph(gmap, nodemap, avg_res, four_neighbors)
        end
        # T = eltype(a)
        if one_to_all
            #source_map = map(x -> x == n ? str : 0, point_map)
            source_map = map(x -> x == n ? T(str) : T(0), unique_point_map)
            ground_map = map(x -> x == n ? T(0) : T(x), point_map)
            map!(x -> x > 0 ? Inf : x, ground_map, ground_map)
        else
            source_map = map(x -> x != 0 ? T(x) : T(0), point_map)
            map!(x -> x == n ? 0 : x, source_map, source_map)
            map!(x -> x != 0 ? 1 : x, source_map, source_map)
            ground_map = map(x -> x == n ? Inf : T(0), point_map)
        end

        check_node = nodemap[points_rc[1][i], points_rc[2][i]]
        
        policy = one_to_all ? :rmvgnd : :rmvsrc
        sources, grounds, finite_grounds = 
                    _get_sources_and_grounds(source_map, ground_map, 
                            flags, G, nodemap, policy)
      
        advanced_data = AdvancedData(G, cc, nodemap, newpoly, hbmeta,
                        sources, grounds, source_map, finite_grounds, 
                        check_node, n, gmap)

        
        if one_to_all
            # v = advanced(cfg, a, source_map, ground_map; nodemap = nodemap, policy = :rmvgnd,
            #                check_node = check_node, src = n, polymap = Polymap(newpoly), hbmeta = hbmeta)
            v, curr = advanced_kernel(advanced_data, flags, cfg)
        else
            # v = advanced(cfg, a, source_map, ground_map; nodemap = nodemap, policy = :rmvsrc,
            #                check_node = check_node, src = n, polymap = Polymap(newpoly), hbmeta = hbmeta)
            v, curr = advanced_kernel(advanced_data, flags, cfg)
        end
        res[i] = v[1]

        cum.cum_curr[mycsid()] .+= curr
        flags.outputflags.write_max_cur_maps && (cum.max_curr[mycsid()] .= max.(cum.max_curr[mycsid()], curr))
    end

    pmap(x -> f(x), 1:num_points_to_solve)

    if flags.outputflags.write_cur_maps
        write_cum_maps(cum, gmap, cfg, hbmeta, 
                       flags.outputflags.write_max_cur_maps, 
                       flags.outputflags.write_cum_cur_map_only)
    end

    hcat(points_unique, res)
end

function prune_points!(points_rc, point_ids::Vector{V}) where V
    rmv = V[]
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

function prune_strengths!(strengths, point_ids::Vector{V}) where V
    pts = strengths[:,1]
    l = length(pts)
    rmv = V[]
    for (i,p) in enumerate(pts)
        if !(p in point_ids)
           push!(rmv, i)
       end
    end
    rng = collect(1:l)
    deleteat!(rng, rmv)
    strengths[rng,:]
end
