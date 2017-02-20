function read_graph(gpath::String)
    i,j,v = load_graph(gpath)
    idx = findfirst(x -> x < 1, i)
    idx != 0 && throw("Indices no good")
    idx = findfirst(x -> x < 1, j)
    idx != 0 && throw("Indices no good")
    m = max(i[end], j[end])
    A = sparse(i,j,v,m,m)
    A + A'
end

function load_graph(gpath::String)
    g = readdlm(gpath)
    i = Int.(g[:,1]) + 1
    j = Int.(g[:,2]) + 1
    v = Float64.(g[:,3])
    i,j,v
end
