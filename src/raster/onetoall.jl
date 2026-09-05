"""
    run_onetoall(data::RasterData, cfg)

One-to-all and all-to-one modes. Every focal point becomes a node of one
shared graph; then for each focal node in turn the source and ground
vectors are built (one-to-all: it is the source and the others are direct
grounds; all-to-one: it is the ground and the others are sources) and the
advanced kernel solves its component. Returns the focal ids with the
effective resistance of each.
"""
function run_onetoall(data::RasterData{T,V}, cfg)::Matrix{T} where {T,V}

    # Data
    strengths = data.strengths
    included_pairs = data.included_pairs
    points_rc = data.points_rc
    gmap = data.cellmap
    polymap = data.polymap
    hbmeta = data.hbmeta

    # Options
    use_variable_strengths = !isempty(strengths)
    use_included_pairs = !isempty(included_pairs)
    mode = included_pairs.mode == :include ? 0 : 1
    one_to_all = is_onetoall(cfg)

    if use_included_pairs
        prune_points!(points_rc, included_pairs.point_ids)
        if use_variable_strengths
            strengths = prune_strengths(strengths, included_pairs.point_ids)
        end
    end

    # Construct point map
    point_map = zeros(V, size(gmap))
    f(i, x) = points_rc[i][x]
    for x = 1:size(points_rc[1], 1)
        point_map[f(1,x), f(2,x)] = f(3, x)
    end

    points_unique = unique(points_rc[3])

    # The shared graph: focal points merged into the polygons
    newpoly = create_new_polymap(gmap, polymap, points_rc, 0, 0, point_map)
    G, cc, geometry = build_graph(data, cfg, newpoly)
    nodemap = geometry.nodemap
    @info("There are $(size(G, 1)) points and $(length(cc)) connected components")

    cum = initialize_cum_maps(gmap, cfg.write_max_cur_maps)

    point_ids = included_pairs.point_ids
    res = zeros(T, size(points_unique, 1))
    num_points_to_solve = size(points_unique, 1)
    original_point_map = copy(point_map)
    unique_point_map = zeros(V, size(gmap))
    strength_map = use_variable_strengths ? zeros(T, size(gmap)) : zeros(T, 0, 0)

    for i in points_unique
        ind = findfirst(x -> x == i, points_rc[3])
        unique_point_map[f(1,ind), f(2,ind)] = f(3,ind)
    end

    # With an include/exclude file the set of focal points that act as
    # sources/grounds changes from one focal point to the next. When every focal
    # point is a single cell that never changes the graph: an inactive point is
    # still an ordinary node, so the graph built above is valid for every solve
    # and only the source/ground vectors differ. Only when a focal point spans
    # several cells (a focal region) does the topology depend on which points
    # are active, because a region collapses into one node only while it is
    # active. That is the one case that needs a rebuild per focal point. (Same
    # test as `has_focal_regions`.)
    point_file_no_polygons = length(points_rc[1]) == length(points_unique)

    # Solve focal point `i`; returns its current map, or nothing if skipped
    function f(i)
        let point_map = copy(original_point_map), strength_map = copy(strength_map), newpoly = newpoly, nodemap = nodemap, G = G, cc = cc
        str = use_variable_strengths ? strengths[i,2] : 1
        @info("Solving point $i of $num_points_to_solve")
        n = points_unique[i]
        if use_included_pairs
            # `i` indexes the focal points present in the raster, but the
            # include matrix is indexed by position in `point_ids`. These differ
            # whenever the file names an ID the point raster does not contain,
            # so look the row up by node ID rather than reusing `i` (issue #341).
            row = findfirst(isequal(n), point_ids)
            for j = 1:size(point_ids,1)
                if j != row && included_pairs.include_pairs[row,j] == mode
                    exclude = point_ids[j]
                    map!(x -> x == exclude ? 0 : x, point_map, point_map)
                end
            end
            if !point_file_no_polygons
                # Focal regions: excluded regions must not be merged into a
                # single node for this solve, so rebuild the graph from the
                # pruned point map.
                newpoly = create_new_polymap(gmap, polymap, points_rc, 0, 0, point_map)
                G, cc, point_geometry = build_graph(data, cfg, newpoly)
                nodemap = point_geometry.nodemap
            end
        end
        if use_variable_strengths
            _tmp = [point_map[points_rc[1][x], points_rc[2][x]] for x = 1:size(points_rc[1], 1)]
            idx = findall(x -> x == 0, _tmp)
            _strengths = deepcopy(strengths)
            _strengths[idx, 2] .= 1
            for x = 1:size(points_rc[1], 1)
                strength_map[points_rc[1][x], points_rc[2][x]] = _strengths[x,2]
            end
        end
        if sum(point_map) == n
            res[i] = -1
            return nothing
        end
        if one_to_all
            source_map = map(x -> x == n ? T(str) : T(0), unique_point_map)
            ground_map = map(x -> x == n ? T(0) : T(x), point_map)
            map!(x -> x > 0 ? Inf : x, ground_map, ground_map)
        else
            if use_variable_strengths
                source_map = map((x,y) -> x == n ? T(0) : T(y), unique_point_map, strength_map)
            else
                source_map = map(x -> x != 0 ? T(1) : T(0), unique_point_map)
                source_map = map((x,y) -> x == n ? T(0) : y, point_map, source_map)
            end
            ground_map = map(x -> x == n ? Inf : T(0), point_map)
        end

        # Only the component holding this focal node is solved. `i` indexes
        # `points_unique`, not the rows of `points_rc`; look the cell up by id
        # so focal regions get the right node.
        ind = findfirst(isequal(n), points_rc[3])
        check_node = nodemap[points_rc[1][ind], points_rc[2][ind]]

        point_geometry = RasterGeometry(nodemap, newpoly, hbmeta, gmap)
        policy = one_to_all ? :rmvgnd : :rmvsrc
        sources, grounds, finite_grounds =
                    sources_and_grounds(point_geometry, source_map, ground_map, G, cfg, policy)

        prob = AdvancedProblem(G, cc, point_geometry, sources, grounds, finite_grounds, get_solver(cfg))
        voltages, curr, solver_called = advanced_kernel(prob, cfg; check_node, name = "_$(V(n))")

        res[i] = onetoall_resistance(voltages, solver_called, one_to_all, n, str,
                                     unique_point_map, nodemap)
        return curr
        end
    end

    is_parallel = cfg.parallelize
    if is_parallel
        results = fetch.(map(x -> Threads.@spawn(f(x)), 1:num_points_to_solve))
    else
        results = map(f, 1:num_points_to_solve)
    end

    # Reduce: accumulate current maps on main thread
    for curr in results
        curr === nothing && continue
        cum.cum_curr .+= curr
        cfg.write_max_cur_maps && (cum.max_curr .= max.(cum.max_curr, curr))
    end

    write_cum_maps(cum, hbmeta, cfg)

    hcat(points_unique, res)
end

# The effective resistance of focal node `n` from one solve: in one-to-all its
# voltage over the injected strength, -1 if nothing was solved or the voltage
# is zero; all-to-one only writes maps and reports 0 for a solved node.
function onetoall_resistance(voltages::Vector{T}, solver_called, one_to_all, n, str,
                             unique_point_map, nodemap) where T
    solver_called || return T(-1)
    one_to_all || return T(0)
    cell = findfirst(isequal(n), unique_point_map)
    node = nodemap[cell]
    val = node == 0 ? T(0) : voltages[node] / T(str)
    val ≈ 0 ? T(-1) : val
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

function prune_strengths(strengths, point_ids::Vector{V}) where V
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
