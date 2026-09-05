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
