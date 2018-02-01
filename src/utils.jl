export  model_problem,
        test_problem

"""
Calculate laplacian of the adjacency matrix of a graph
"""
function laplacian(G)
    n = size(G, 1)
    s = Vector{eltype(G)}(n)
    for i = 1:n
        s[i] = sum_off_diag(G, i)
        for j in nzrange(G, i)
            if i == G.rowval[j]
                G.nzval[j] = 0
            else
                G.nzval[j] = -G.nzval[j]
            end
        end
    end
    G + spdiagm(s)
end

function sum_off_diag(G, i)
     sum = zero(eltype(G))
     for j in nzrange(G, i)
         if G.rowval[j] != i
             sum += G.nzval[j]
         end
     end
     sum
 end

 """
 Construct nodemap specific to a connected component
 """
function construct_local_node_map(nodemap, component, polymap)
    local_nodemap = zeros(eltype(nodemap), size(nodemap))
    idx = findin(nodemap, component)
    local_nodemap[idx] = nodemap[idx]
    _construct_local_nodemap(local_nodemap, polymap, idx)
end

function _construct_local_nodemap(local_nodemap, polymap, idx)
    if isempty(polymap)
        i = find(local_nodemap)
        local_nodemap[idx] = 1:length(i)
        return local_nodemap
    else
        local_polymap = zeros(eltype(local_nodemap), size(local_nodemap))
        local_polymap[idx] = polymap[idx]
        return construct_node_map(local_nodemap, local_polymap)
    end
end

"""
Define model circuitscape problem - helps in tests
"""
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

# Testing utilities 

function test_problem(str)
    base_path = joinpath(Pkg.dir("Circuitscape"), "test", "input")
    str2 = replace(str, ".ini", "")
    if contains(str, "sgVerify")
        config_path = joinpath(base_path, "raster", "pairwise", replace(str2, "sgVerify", ""), str)
    elseif contains(str, "mgVerify")
        config_path = joinpath(base_path, "raster", "advanced", replace(str2, "mgVerify", ""), str)
    elseif contains(str, "oneToAll")
        config_path = joinpath(base_path, "raster", "one_to_all", replace(str2, "oneToAllVerify", ""), str)
    elseif contains(str, "allToOne")
        config_path = joinpath(base_path, "raster", "all_to_one", replace(str2, "allToOneVerify", ""), str)
    elseif contains(str, "sgNetworkVerify")
        config_path = joinpath(base_path, "network", replace(str2, "sgNetworkVerify", ""), str)
    elseif contains(str, "mgNetworkVerify")
        config_path = joinpath(base_path, "network", replace(str2, "mgNetworkVerify", ""), str)
    end
end    
