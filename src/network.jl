using LightGraphs
using Logging
using IterativeSolvers
Logging.configure(level = DEBUG)

function thing(network_file, current_file)

    # Read network file
    a = readdlm(network_file, ' ')

    # Create sparse matrix
    i = Int.(a[:,1] + 1)
    j = Int.(a[:,2] + 1)
    val = float.(a[:,3] + 1)
    m = max(i[end], j[end])
    A = sparse(i, j, val, m, m)
    A = A + A'

    # Create graph
    g = Graph(A)

    # Read currents
    c = Int.(vec(readcsv(current_file)) + 1)

    voltages = single_ground_all_pair_resistances(g, c)

    A,c, voltages
end

function single_ground_all_pair_resistances{T}(g::Graph, c::Vector{T})
    numpoints = size(c, 1)
    debug("There are $numpoints points")
    cc = connected_components(g)
    debug("There are $(length(cc)) connected components")

    volt = Vector{Float64}(size(g, 1))
    total = Int(numpoints * (numpoints-1) / 2)
    
    p = 1
    for i = 1:numpoints
        for j = i+1:numpoints
            pt1 = c[i]
            pt2 = c[j]
            info("Solving for pair $p of $total")
            debug("pt1 = $pt1, pt2 = $pt2")
            p +=1
            curr = zeros(size(g, 1))
            curr[pt1] = -1
            curr[pt2] = 1
            cond = laplacian_matrix(g)
            volt = cg(cond, curr)
        end
    end
end
