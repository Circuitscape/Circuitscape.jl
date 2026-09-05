# Graph construction for a network: the edge list becomes a symmetric
# adjacency matrix, and the output geometry is the identity on node ids.

function build_graph(data::NetworkData{T,V}, cfg) where {T,V}

    i,j,v = data.coords

    idx = findfirst(x -> x < 1, i)
    idx !== nothing && throw("Indices no good")
    idx = findfirst(x -> x < 1, j)
    idx !== nothing && throw("Indices no good")

    m = max(maximum(i), maximum(j))
    A = sparse(i,j,v,m,m)
    A = A + A'

    cc = connected_components(A)
    c = size(A,1)
    @info("Graph has $c nodes and $(length(cc)) connected components")

    G = @timeit CSTIMER[] "construct graph laplacian" laplacian!(A)

    G, cc, NetworkGeometry(collect(V, 1:m))
end

# Focal nodes are the user's node ids themselves. In include mode focal
# points absent from the file are dropped, mirroring the raster path
# (issue #341).
function focal_nodes(data::NetworkData{T,V}, ::NetworkGeometry) where {T,V}
    included_pairs = data.included_pairs
    fp = data.fp
    if isempty(included_pairs)
        exclude_pairs = Tuple{V,V}[]
    else
        included_pairs.mode == :include &&
            (fp = filter(x -> x in included_pairs.point_ids, fp))
        exclude_pairs = exclude_pairs_from(included_pairs)
    end
    fp, fp, exclude_pairs
end

initialize_cum(data::NetworkData, cfg, num_nodes) = initialize_cum_vectors(data.coords, num_nodes)

has_focal_regions(::NetworkData) = false
