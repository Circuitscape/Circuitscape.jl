export model_problem

function model_problem(T, s)
    
    # Cell map is the uniform 1 across the grid
    cellmap = ones(T, s, s)
    
    # Nodemap is just 1:endof(cellmap)
    nodemap = reshape(1:endof(cellmap), s, s) |> collect

    # Construct graph with 4 neighbors
    G = construct_graph(cellmap, nodemap, true, true)

    # Return laplacian
    laplacian(G)
end
model_problem(s::Integer) = model_problem(Float64, s)
    