function network_advanced(T, V, cfg)::Matrix{T}

    # Get data
    data = get_network_data(T, V, cfg)

    # Get flags
    flags = get_network_flags(cfg)

    # Compute advanced data
    advanced_data = compute_advanced_data(data, flags)

    # Send to main kernel
    v , _ = advanced_kernel(advanced_data, flags, cfg)

    v
end

function compute_advanced_data(data::NetworkData{T,V}, 
                            flags)::AdvancedData{T,V} where {T,V}

    i,j,v = data.coords
    
    idx = findfirst(x -> x < 1, i)
    idx != nothing && throw("Indices no good")
    idx = findfirst(x -> x < 1, j)
    idx != nothing && throw("Indices no good")
    
    m = max(i[end], j[end])
    A = sparse(i,j,v,m,m)
    A = A + A'

    cc = connected_components(SimpleWeightedGraph(A))
    
    t = @elapsed G = laplacian(A)
    csinfo("Time taken to construct graph laplacian = $t")

    nodemap, polymap = Matrix{V}(undef,0,0), Matrix{V}(undef,0,0)
    cellmap = Matrix{T}(undef,0,0)

    sources, grounds, finite_grounds = 
                get_sources_and_grounds(data, flags, G, nodemap)

    source_map = Matrix{eltype(A)}(undef,0,0)
    AdvancedData(G, cc, nodemap, polymap, RasterMeta(), 
                 sources, grounds, source_map, finite_grounds, V(-1), V(0), cellmap)
end
