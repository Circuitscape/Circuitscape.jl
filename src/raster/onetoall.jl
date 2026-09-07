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
    included_pairs = data.included_pairs
    points_rc = data.points_rc
    gmap = data.cellmap
    polymap = data.polymap
    hbmeta = data.hbmeta

    # Options
    use_variable_strengths = !isempty(data.strengths)
    use_included_pairs = !isempty(included_pairs)
    mode = included_pairs.mode == :include ? 0 : 1
    one_to_all = is_onetoall(cfg)

    use_included_pairs && prune_points!(points_rc, included_pairs.point_ids)

    # Variable source strengths keyed by focal id, as Python's
    # `get_strengths_rc`: the row order of the file is irrelevant, an id the
    # file does not list gets strength 1, and ids the point file does not
    # contain are ignored (so nothing has to be pruned against the include
    # file). The first row wins when an id is repeated.
    strength_of = Dict{V,T}()
    if use_variable_strengths
        for r in axes(data.strengths, 1)
            get!(strength_of, V(data.strengths[r,1]), data.strengths[r,2])
        end
    end

    # Construct point map
    point_map = zeros(V, size(gmap))
    _pt(i, x) = points_rc[i][x]
    for x = 1:size(points_rc[1], 1)
        point_map[_pt(1,x), _pt(2,x)] = _pt(3, x)
    end

    points_unique = unique(points_rc[3])

    # The shared graph: focal points merged into the polygons
    newpoly = create_new_polymap(gmap, polymap, points_rc, 0, 0, point_map)
    G, cc, geometry = build_graph(data, cfg, newpoly)
    nodemap = geometry.nodemap
    @info("There are $(size(G, 1)) points and $(length(cc)) connected components")

    cum = initialize_cum_maps(V, gmap, cfg.write_max_cur_maps)

    res = zeros(T, size(points_unique, 1))
    num_points_to_solve = size(points_unique, 1)
    original_point_map = copy(point_map)

    # One representative cell per focal id (its first cell in `points_rc`),
    # aligned with `points_unique`; `unique_point_map` marks those cells.
    unique_cells = map(points_unique) do id
        ind = findfirst(isequal(id), points_rc[3])
        (_pt(1,ind), _pt(2,ind))
    end
    unique_point_map = zeros(V, size(gmap))
    for (k, (r, c)) in enumerate(unique_cells)
        unique_point_map[r, c] = points_unique[k]
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

    # Solve focal point `i`; returns its current map, or nothing if skipped.
    # A top-level function rather than a closure: the per-point body rebuilds
    # the graph for focal regions, and a closure that reassigns captured
    # variables boxes them, which made the whole body dynamically typed.
    solve_point(i) = solve_onetoall_point(i, data, cfg, G, cc, nodemap, newpoly,
        original_point_map, unique_point_map, unique_cells, points_unique,
        strength_of, res, mode, use_variable_strengths, use_included_pairs,
        point_file_no_polygons, one_to_all)

    is_parallel = cfg.parallelize
    if is_parallel
        results = fetch.(map(x -> Threads.@spawn(solve_point(x)), 1:num_points_to_solve))
    else
        results = map(solve_point, 1:num_points_to_solve)
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

# One focal point of `run_onetoall`: build its source and ground maps (and,
# for focal regions with an include/exclude file, its own graph), solve its
# component with the advanced kernel, record its effective resistance in
# `res[i]` and return its current map, or nothing if the point was skipped.
# `point_map` is copied because the include/exclude handling edits it per
# point; every other argument is read only, so the function is safe to run
# on several threads at once.
function solve_onetoall_point(i, data::RasterData{T,V}, cfg, G, cc, nodemap, newpoly,
                              original_point_map, unique_point_map, unique_cells,
                              points_unique, strength_of, res, mode,
                              use_variable_strengths, use_included_pairs,
                              point_file_no_polygons, one_to_all) where {T,V}
    included_pairs = data.included_pairs
    point_ids = included_pairs.point_ids
    points_rc = data.points_rc
    gmap = data.cellmap
    polymap = data.polymap
    hbmeta = data.hbmeta
    num_points_to_solve = length(points_unique)

    point_map = copy(original_point_map)
    @info("Solving point $i of $num_points_to_solve")
    n = points_unique[i]
    str = use_variable_strengths ? get(strength_of, n, one(T)) : one(T)
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
    # Nothing to solve when every other focal point is excluded for this one
    # (Python: `unique_point_map.sum() == src`). Test the cells rather than
    # sum them so a lone multi-cell focal region is recognised as well.
    if all(x -> x == 0 || x == n, point_map)
        res[i] = -1
        return nothing
    end
    if one_to_all
        source_map = map(x -> x == n ? str : zero(T), unique_point_map)
        ground_map = map(x -> x == n ? T(0) : T(x), point_map)
        map!(x -> x > 0 ? Inf : x, ground_map, ground_map)
    else
        if use_variable_strengths
            # Python's `get_strength_map`: one entry per focal id at its
            # representative cell. A point the include file drops for this
            # solve keeps a unit strength (Python zeroes its id before the
            # lookup, which then defaults to 1), so with variable strengths
            # an excluded point still acts as a unit source. This reproduces
            # Python's own behaviour and allToOneVerify12 depends on it; the
            # unit-strength branch below does exclude such points.
            strength_map = zeros(T, size(gmap))
            for (k, (r, c)) in enumerate(unique_cells)
                strength_map[r, c] = point_map[r, c] == 0 ? one(T) :
                                     get(strength_of, points_unique[k], one(T))
            end
            source_map = map((x,y) -> x == n ? T(0) : T(y), unique_point_map, strength_map)
        else
            # Every active focal point other than `n` is a unit source; the
            # points the include file excludes are zero in `point_map`.
            source_map = map((x,y) -> (x != 0 && y != 0 && x != n) ? one(T) : zero(T),
                             unique_point_map, point_map)
        end
        ground_map = map(x -> x == n ? Inf : T(0), point_map)
    end

    # Only the component holding this focal node is solved. `i` indexes
    # `points_unique`, not the rows of `points_rc`; use its representative
    # cell so focal regions get the right node.
    check_node = nodemap[unique_cells[i]...]

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

# The effective resistance of focal node `n` from one solve: in one-to-all its
# voltage over the injected strength, -1 if nothing was solved; all-to-one
# only writes maps and reports 0 for a solved node.
function onetoall_resistance(voltages::Vector{T}, solver_called, one_to_all, n, str,
                             unique_point_map, nodemap) where T
    solver_called || return T(-1)
    one_to_all || return T(0)
    cell = findfirst(isequal(n), unique_point_map)
    node = nodemap[cell]
    node == 0 ? T(0) : voltages[node] / T(str)
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
