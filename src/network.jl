using Logging
using LightGraphs
using IterativeSolvers
Logging.configure(level = DEBUG)

function thing(network_file, current_file)

    # Read network file
    a = readdlm(network_file, ' ')

    # Create sparse matrix
    i = Int.(a[:,1] + 1)
    j = Int.(a[:,2] + 1)
    val = float.(a[:,3])
    m = max(i[end], j[end])
    A = sparse(i, j, val, m, m)
    A = A + A'

    # Create graph
    g = Graph(A)

    # Read currents
    c = Int.(vec(readcsv(current_file)) + 1)

    voltages = single_ground_all_pair_resistances(A, g, c)

    A,c, voltages
end

function single_ground_all_pair_resistances{T}(a::SparseMatrixCSC, g::Graph, c::Vector{T})
    numpoints = size(c, 1)
    debug("There are $numpoints points")
    cc = connected_components(g)
    debug("There are $(length(cc)) connected components")

    cond = laplacian(a)

    volt = Vector{Float64}(size(g, 1))
    total = Int(numpoints * (numpoints-1) / 2)
    
    p = 1
    for i = 1:numpoints
        pt1 = c[i]
        d = cond[pt1, pt1]
        cond[pt1, pt1] = 0
        for j = i+1:numpoints
            pt2 = c[j]
            info("Solving for pair $p of $total")
            debug("pt1 = $pt1, pt2 = $pt2")
            p +=1
            curr = zeros(size(g, 1))
            curr[pt1] = -1
            curr[pt2] = 1
            println('\n')
            volt = cg(cond, curr)
        end
        cond[pt1,pt1] = d
    end
end

function laplacian(G::SparseMatrixCSC)
    G = G - spdiagm(diag(G))
    G = -G + spdiagm(vec(sum(G, 1)))
end
