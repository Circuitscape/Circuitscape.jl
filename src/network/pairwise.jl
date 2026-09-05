"""
Primary driver for network pairwise.
"""
function network_pairwise(T, V, cfg)::Matrix{T}

    # Get input
    networkdata = get_network_data(T, V, cfg)

    # Compute graph data
    graphdata = compute_graph_data(networkdata, cfg)

    # Send to main kernel
    ret = single_ground_all_pairs(graphdata, cfg)

	# Write cum maps
	if cfg.write_cur_maps
		cum_node_curr = graphdata.cum.cum_node_curr
		cum_branch_curr = graphdata.cum.cum_branch_curr
		cum_node_curr = hcat(1:length(cum_node_curr), cum_node_curr)
		coords = graphdata.cum.coords
		cum_branch_curr = hcat(getindex.(coords, 1), getindex.(coords, 2), cum_branch_curr)
		write_currents(cum_node_curr, cum_branch_curr, "_cum", cfg)
	end

	ret
end

function compute_graph_data(data::NetworkData{T,V}, cfg)::GraphProblem{T,V} where {T,V}


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

    # Include/exclude pairs. In include mode focal points absent from the file
    # are dropped, mirroring the raster path (issue #341).
    included_pairs = data.included_pairs
    fp = data.fp
    if isempty(included_pairs)
        exclude_pairs = Tuple{V,V}[]
    else
        included_pairs.mode == :include &&
            (fp = filter(x -> x in included_pairs.point_ids, fp))
        exclude_pairs = exclude_pairs_from(included_pairs)
    end
    solver = get_solver(cfg)

    nodemap = Matrix{V}(undef,0,0)
    polymap = Matrix{V}(undef,0,0)
    hbmeta = RasterMeta()
    cellmap = Matrix{T}(undef,0,0)

	cum = initialize_cum_vectors(data.coords, size(G,1))

    GraphProblem(G, cc, fp, fp, 
                exclude_pairs, nodemap, polymap, hbmeta, cellmap, cum, solver)

end

