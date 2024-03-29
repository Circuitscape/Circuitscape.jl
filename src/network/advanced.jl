function network_advanced(T, V, cfg)::Matrix{T}

    # Get data
    data = get_network_data(T, V, cfg)

    # Get flags
    flags = get_network_flags(cfg)

    # Compute advanced data
    advanced_data = compute_advanced_data(data, flags, cfg)

    # Send to main kernel
    v , _ = advanced_kernel(advanced_data, flags, cfg)

    v
end


function compute_advanced_data(data::NetworkData{T,V}, 
                            flags, cfg)::AdvancedProblem{T,V} where {T,V}

    i,j,v = data.coords

    idx = findfirst(x -> x < 1, i)
    idx != nothing && throw("Indices no good")
    idx = findfirst(x -> x < 1, j)
    idx != nothing && throw("Indices no good")

    m = max(maximum(i), maximum(j))
    A = sparse(i,j,v,m,m)
    A = A + A'

    cc = connected_components(SimpleWeightedGraph(A))
	c = size(A,1)
	csinfo("Graph has $c nodes and $(length(cc)) connected components", cfg["suppress_messages"] in TRUELIST)

    t = @elapsed G = laplacian(A)
    csinfo("Time taken to construct graph laplacian = $t", cfg["suppress_messages"] in TRUELIST)

    nodemap, polymap = Matrix{V}(undef,0,0), Matrix{V}(undef,0,0)
    cellmap = Matrix{T}(undef,0,0)

    solver = get_solver(cfg)

    sources, grounds, finite_grounds = 
                get_sources_and_grounds(data, flags, G, nodemap)

    source_map = Matrix{eltype(A)}(undef,0,0)
    AdvancedProblem(G, cc, nodemap, polymap, RasterMeta(), 
                 sources, grounds, source_map, finite_grounds, V(-1), V(0), cellmap, solver)

end
