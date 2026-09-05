# Graph construction for a raster landscape: one node per habitat cell, the
# cells of one short-circuit polygon merged into a node, edges to the 4 or 8
# neighbours. The output geometry is the nodemap this produces.

"""
    build_graph(data::RasterData, cfg, polymap = data.polymap) -> (G, cc, geometry)

Laplacian, connected components and geometry of the landscape graph. A
`polymap` other than the loaded one lets callers merge extra cells into
nodes, as the focal-region paths do.
"""
function build_graph(data::RasterData{T,V}, cfg, polymap = data.polymap) where {T,V}
    cellmap = data.cellmap
    nodemap = construct_node_map(cellmap, polymap)
    A = construct_graph(cellmap, nodemap, cfg.connect_using_avg_resistances,
                        cfg.connect_four_neighbors_only)
    cc = connected_components(A)
    G = laplacian!(A)
    G, cc, RasterGeometry(nodemap, polymap, data.hbmeta, cellmap)
end

# Focal nodes are the nodes under the focal cells, paired with their ids.
function focal_nodes(data::RasterData{T,V}, geometry::RasterGeometry) where {T,V}
    points_rc = data.points_rc
    included_pairs = data.included_pairs
    exclude_pairs = isempty(included_pairs) ? Tuple{V,V}[] :
                        generate_exclude_pairs(points_rc, included_pairs)
    nodemap = geometry.nodemap
    points = V[nodemap[i, j] for (i, j) in zip(points_rc[1], points_rc[2])]
    points, points_rc[3], exclude_pairs
end

initialize_cum(data::RasterData, cfg, num_nodes) =
    initialize_cum_maps(data.cellmap, cfg.write_max_cur_maps)

# A point file that names one id in several cells describes focal regions.
has_focal_regions(data::RasterData) =
    length(data.points_rc[1]) != length(unique(data.points_rc[3]))

"""
    pairwise_regions(rasterdata, cfg)

Pairwise mode over focal regions: every pair gets its own graph, with the
two regions being solved collapsed into nodes and the others left as
ordinary habitat.
"""
function pairwise_regions(rasterdata::RasterData{T,V}, cfg)::Matrix{T} where {T,V}

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
        graphdata = region_pair_problem(rasterdata, pt1, pt2, cum, cfg)
        pairwise_resistance = single_ground_all_pairs(graphdata, cfg, false)
        resistances[i,j] = resistances[j,i] = pairwise_resistance[2,3]
    end
    for i = 1:size(pts, 1)
        resistances[i,i] = 0
    end
    P = [0, pts...]
    r = hcat(P, vcat(pts', resistances))

    if cfg.write_cur_maps || cfg.write_cum_cur_map_only
        write_cum_maps(cum, rasterdata.hbmeta, cfg)
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


function region_pair_problem(rasterdata::RasterData{T,V},
                            pt1, pt2, cum, cfg)::GraphProblem{T,V} where {T,V}

    points_rc = rasterdata.points_rc

    # Construct new polymap and the graph over it
    newpoly = create_new_polymap(rasterdata.cellmap, rasterdata.polymap, points_rc, pt1, pt2)
    G, cc, geometry = build_graph(rasterdata, cfg, newpoly)

    # Construct points vector
    nodemap = geometry.nodemap
    x = something(findfirst(isequal(pt1), points_rc[3]), 0)
    y = something(findfirst(isequal(pt2), points_rc[3]), 0)
    c1 = nodemap[points_rc[1][x], points_rc[2][x]]
    c2 = nodemap[points_rc[1][y], points_rc[2][y]]
    points = V[c1, c2]

    # Exclude pairs array
    exclude_pairs = Tuple{V,V}[]

    GraphProblem(G, cc, points, [pt1, pt2], exclude_pairs, geometry, cum, get_solver(cfg))
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

"""
    construct_graph(gmap, nodemap, avg_res, four_neighbors)

Build the symmetric weighted adjacency matrix of the raster graph. Every
non-zero cell of `nodemap` is a node; each node is connected to its right and
lower neighbours (and, unless `four_neighbors`, to its two right-hand diagonal
neighbours). Edge weights are the conductance between the two cells, averaged
either as resistances (`avg_res`) or as conductances.

Edges are counted in a first pass so that the COO buffers can be allocated
exactly, and both orientations of every edge are emitted in the second pass so
that `sparse` yields the symmetric matrix directly. This avoids the repeated
`push!` reallocations and the `a + a'` step (a second full matrix plus its
transpose) of the previous implementation, cutting peak memory at load time,
which is where the largest grids used to run out of memory.
"""
function construct_graph(gmap, nodemap::Matrix{S}, avg_res, four_neighbors) where S
    f1 = avg_res ? res_avg : cond_avg
    f2 = avg_res ? weirder_avg : weird_avg
    nrows, ncols = size(gmap)
    size(nodemap) == (nrows, ncols) ||
        throw(DimensionMismatch("gmap is $(size(gmap)) but nodemap is $(size(nodemap))"))

    nedges = count_graph_edges(nodemap, four_neighbors)

    I = Vector{S}(undef, 2nedges)
    J = Vector{S}(undef, 2nedges)
    V = Vector{eltype(gmap)}(undef, 2nedges)

    k = 0
    @inbounds for j = 1:ncols
        for i = 1:nrows
            n = nodemap[i,j]
            n == 0 && continue
            g = gmap[i,j]

            # Horizontal neighbour
            if j != ncols && nodemap[i,j+1] != 0
                k = _emit_edge!(I, J, V, k, n, nodemap[i,j+1], f1(g, gmap[i,j+1]))
            end

            # Vertical neighbour
            if i != nrows && nodemap[i+1,j] != 0
                k = _emit_edge!(I, J, V, k, n, nodemap[i+1,j], f1(g, gmap[i+1,j]))
            end

            if !four_neighbors
                # Diagonal neighbours
                if i != nrows && j != ncols && nodemap[i+1,j+1] != 0
                    k = _emit_edge!(I, J, V, k, n, nodemap[i+1,j+1], f2(g, gmap[i+1,j+1]))
                end

                if i != 1 && j != ncols && nodemap[i-1,j+1] != 0
                    k = _emit_edge!(I, J, V, k, n, nodemap[i-1,j+1], f2(g, gmap[i-1,j+1]))
                end
            end
        end
    end
    @assert k == 2nedges

    m = maximum(nodemap)

    # `sparse!` is the driver behind `sparse`; it lets us hand it the COO
    # buffers J and V as the storage for the result's rowval and nzval, so no
    # fresh CSC is allocated. colptr gets its own (m + 1)-vector: recycling I
    # for it would shrink I's length but keep its 2*nedges buffer alive.
    # (`sparse!` sums duplicate (i, j) pairs, exactly as `a + a'` did when a
    # polygon maps adjacent cells to the same node.)
    csrrowptr = Vector{S}(undef, m + 1)
    csrcolval = Vector{S}(undef, length(I))
    csrnzval = Vector{eltype(gmap)}(undef, length(I))
    klasttouch = Vector{S}(undef, m)
    colptr = Vector{S}(undef, m + 1)
    SparseArrays.sparse!(I, J, V, m, m, +, klasttouch,
                         csrrowptr, csrcolval, csrnzval, colptr, J, V)
end

# Write edge (a, b, v) in both orientations at positions k+1 and k+2.
@inline function _emit_edge!(I, J, V, k, a, b, v)
    @inbounds begin
        I[k+1] = a; J[k+1] = b; V[k+1] = v
        I[k+2] = b; J[k+2] = a; V[k+2] = v
    end
    k + 2
end

"""
    count_graph_edges(nodemap, four_neighbors)

Number of undirected edges `construct_graph` will emit for `nodemap`, using the
same neighbour rules, so that its output buffers can be sized exactly.
"""
function count_graph_edges(nodemap, four_neighbors)
    nrows, ncols = size(nodemap)
    n = 0
    @inbounds for j = 1:ncols
        for i = 1:nrows
            nodemap[i,j] == 0 && continue
            if j != ncols && nodemap[i,j+1] != 0
                n += 1
            end
            if i != nrows && nodemap[i+1,j] != 0
                n += 1
            end
            if !four_neighbors
                if i != nrows && j != ncols && nodemap[i+1,j+1] != 0
                    n += 1
                end
                if i != 1 && j != ncols && nodemap[i-1,j+1] != 0
                    n += 1
                end
            end
        end
    end
    n
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
