function read_graph(a::Inifile, gpath::String)
    i,j,v = load_graph(gpath)
    idx = findfirst(x -> x < 1, i)
    idx != 0 && throw("Indices no good")
    idx = findfirst(x -> x < 1, j)
    idx != 0 && throw("Indices no good")
    is_res = get(a, "Habitat raster or graph", "habitat_map_is_resistances")
    if is_res == "True"
        v = 1./v
    end
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

read_focal_points(path::String) = Int.(vec(readcsv(path)) + 1)

function read_point_strengths(path::String)
    a = readdlm(path)
    a[:,1] = a[:,1] + 1
    a
end
