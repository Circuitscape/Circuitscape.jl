# Grid-to-graph construction for a raster landscape: one node per habitat
# cell, the cells of one short-circuit polygon merged into a node, edges to
# the 4 or 8 neighbours weighted by the average conductance of the two cells.
# Focal points and regions are merged into nodes by editing the polymap.

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
non-zero cell of `nodemap` is a node; each cell with positive conductance is
connected to its right and lower neighbours (and, unless `four_neighbors`, to
its two right-hand diagonal neighbours) that also have positive conductance,
so a nodata cell of a polygon belongs to the polygon's node but has no edges. Edge weights are the conductance between the two cells, averaged
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

    nedges = count_graph_edges(gmap, nodemap, four_neighbors)

    I = Vector{S}(undef, 2nedges)
    J = Vector{S}(undef, 2nedges)
    V = Vector{eltype(gmap)}(undef, 2nedges)

    k = 0
    @inbounds for j = 1:ncols
        for i = 1:nrows
            _linked(gmap, nodemap, i, j) || continue
            n = nodemap[i,j]
            g = gmap[i,j]

            # Horizontal neighbour
            if j != ncols && _linked(gmap, nodemap, i, j+1)
                k = _emit_edge!(I, J, V, k, n, nodemap[i,j+1], f1(g, gmap[i,j+1]))
            end

            # Vertical neighbour
            if i != nrows && _linked(gmap, nodemap, i+1, j)
                k = _emit_edge!(I, J, V, k, n, nodemap[i+1,j], f1(g, gmap[i+1,j]))
            end

            if !four_neighbors
                # Diagonal neighbours
                if i != nrows && j != ncols && _linked(gmap, nodemap, i+1, j+1)
                    k = _emit_edge!(I, J, V, k, n, nodemap[i+1,j+1], f2(g, gmap[i+1,j+1]))
                end

                if i != 1 && j != ncols && _linked(gmap, nodemap, i-1, j+1)
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
    count_graph_edges(gmap, nodemap, four_neighbors)

Number of undirected edges `construct_graph` will emit for `gmap` and
`nodemap`, using the same neighbour rules (both cells must be nodes with
nonzero conductance), so that its output buffers can be sized exactly.
"""
function count_graph_edges(gmap, nodemap, four_neighbors)
    nrows, ncols = size(nodemap)
    n = 0
    @inbounds for j = 1:ncols
        for i = 1:nrows
            _linked(gmap, nodemap, i, j) || continue
            if j != ncols && _linked(gmap, nodemap, i, j+1)
                n += 1
            end
            if i != nrows && _linked(gmap, nodemap, i+1, j)
                n += 1
            end
            if !four_neighbors
                if i != nrows && j != ncols && _linked(gmap, nodemap, i+1, j+1)
                    n += 1
                end
                if i != 1 && j != ncols && _linked(gmap, nodemap, i-1, j+1)
                    n += 1
                end
            end
        end
    end
    n
end

# Whether cell (i, j) takes part in edges: it must be a node and have positive
# conductance. A cell of a polygon has the polygon's node id even when it is
# nodata (conductance 0); like Circuitscape 4's _neighbors_* (which link two
# cells only when both conductances are nonzero), such a cell is not linked to
# its neighbours, so a polygon straddling nodata does not reach across it.
# `> 0` rather than `!= 0` matches `construct_node_map`, and keeps negative or
# -9999 cells passed in by API callers (e.g. Omniscape) out of the graph.
@inline _linked(gmap, nodemap, i, j) =
    @inbounds nodemap[i,j] != 0 && gmap[i,j] > 0

res_avg(x, y) = 1 / ((1/x + 1/y) / 2)
cond_avg(x, y) = (x + y) / 2
weird_avg(x,y) = (x + y) / (2*√2)
weirder_avg(x, y) = 1 / (√2 * (1/x + 1/y) / 2)

"""
    merge_regions_with_polygons!(newpoly, polymap, point_map)

Write into `newpoly` the polygon map for one-to-all and all-to-one runs with
focal regions and polygons, as Circuitscape 4's `get_poly_map_temp` builds it:
for each point id in ascending order, the polygons it overlaps on the running
map and all of its own cells get a fresh id above every polygon number. Two
regions that share a polygon therefore collapse into one node.

Relabelling the whole map once per point id costs O(points × cells), which is
seconds per call on large grids with many regions (and the call is repeated
for every focal node when an include/exclude file is used). Instead, each
point id's fresh label absorbs the current root of every polygon label it
overlaps (a union-find over labels; a parent is always a later, larger label,
so the forest has no cycles), and one final pass writes every cell's root.
"""
function merge_regions_with_polygons!(newpoly::AbstractMatrix{V}, polymap,
                                      point_map) where V
    k0 = V(maximum(polymap))
    # (point id, polygon label) for every cell of a region lying on a polygon
    overlaps = Tuple{V,V}[]
    ids = V[]
    @inbounds for i in eachindex(point_map, polymap)
        p = point_map[i]
        p == 0 && continue
        push!(ids, p)
        polymap[i] != 0 && push!(overlaps, (V(p), V(polymap[i])))
    end
    unique!(sort!(ids))
    unique!(sort!(overlaps))
    label = Dict{V,V}(p => k0 + V(n) for (n, p) in enumerate(ids))

    parent = Dict{V,V}()
    function root(x)
        while true
            y = get(parent, x, x)
            y == x && return x
            x = y
        end
    end
    for (p, l) in overlaps   # ascending point id, as in get_poly_map_temp
        r = root(l)
        k = label[p]
        r != k && (parent[r] = k)
    end

    @inbounds for i in eachindex(newpoly, point_map, polymap)
        p = point_map[i]
        if p != 0
            newpoly[i] = root(label[p])
        elseif polymap[i] != 0
            newpoly[i] = isempty(parent) ? polymap[i] : root(V(polymap[i]))
        else
            newpoly[i] = 0
        end
    end
    newpoly
end

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
            merge_regions_with_polygons!(newpoly, polymap, point_map)
        end
        return newpoly
    end

    if isempty(polymap)
        newpoly = zeros(V, size(gmap)...)
        id1 = findall(x -> x == pt1, points_rc[3])
        id2 = findall(x -> x == pt2, points_rc[3])
        # Plain loops rather than `map` over closures: `newpoly` is assigned
        # on several branches, and a captured variable that is reassigned is
        # boxed, which made the return type `Any` for every caller.
        for x in id1; newpoly[f(x)...] = pt1; end
        for x in id2; newpoly[f(x)...] = pt2; end
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
            # A focal region becomes one node together with every polygon
            # it overlaps, as in Circuitscape 4's get_overlap_polymap: the
            # overlapped polygons and the region's own cells all get a new
            # id. Detecting the overlap on `newpoly` rather than `polymap`
            # means the second region merges with the first if they share
            # a polygon, again as in the original.
            k += 1
            # A loop, not a generator: a generator over `newpoly` captures
            # it, and since `newpoly` is assigned on several branches it is
            # then boxed and the function's return type becomes `Any`.
            vals = V[]
            for x in idx
                v = newpoly[f(x)...]
                v != 0 && !(v in vals) && push!(vals, v)
            end
            if !isempty(vals)
                for i in eachindex(newpoly)
                    newpoly[i] in vals && (newpoly[i] = k)
                end
            end
            for x in idx; newpoly[f(x)...] = k; end
        end
        return newpoly
    end
end
