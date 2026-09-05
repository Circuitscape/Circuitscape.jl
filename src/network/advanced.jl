function network_advanced(T, V, cfg)::Matrix{T}

    # Get data
    data = get_network_data(T, V, cfg)

    # Compute advanced data
    advanced_data = compute_advanced_data(data, cfg)

    # Send to main kernel
    v , _ = advanced_kernel(advanced_data, cfg)

    v
end


function compute_advanced_data(data::NetworkData{T,V}, cfg)::AdvancedProblem{T,V} where {T,V}

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

    geometry = NetworkGeometry(collect(V, 1:size(G,1)))

    solver = get_solver(cfg)

    sources, grounds, finite_grounds =
                get_sources_and_grounds(data, cfg, G, geometry)

    source_map = Matrix{eltype(A)}(undef,0,0)
    AdvancedProblem(G, cc, geometry,
                 sources, grounds, source_map, finite_grounds, V(-1), V(0), solver)

end
