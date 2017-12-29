function model_problem(size, el = Float64)
    
    # Cell map is the uniform 1 across the grid
    cellmap = ones(el, size, size)
    
    # Nodemap is just 1:endof(cellmap)
    nodemap = reshape(1:endof(cellmap), size(cellmap)...)

    # Construct graph with 4 neighbors
    G = construct_graph(cellmap, nodemap, true, false)

    # Return laplacian
    laplacian(G)
end
    