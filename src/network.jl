using Logging
using LightGraphs
using IterativeSolvers
Logging.configure(level = DEBUG)

function network(network_file, current_file)

    # Read network file
    A = read_graph(network_file)

    # Create graph
    g = Graph(A)

    # Read currents
    c = Int.(vec(readcsv(current_file)) + 1)

    resistances = single_ground_all_pair_resistances(A, g, c)

    resistances
end

function single_ground_all_pair_resistances{T}(a::SparseMatrixCSC, g::Graph, c::Vector{T})
    numpoints = size(c, 1)
    debug("There are $numpoints points")
    cc = connected_components(g)
    debug("There are $(length(cc)) connected components")
    resistances = -1 * ones(numpoints, numpoints) 

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
            println("Conductance matrix")
            Base.print_matrix(STDOUT, cond)
            println('\n')
            @show curr
            volt = gmres(cond, curr)
            @show volt
            postprocess(volt[1], c, pt1, pt2, resistances)
        end
        cond[pt1,pt1] = d
    end
    for i = 1:size(resistances,1)
        resistances[i,i] = 0
    end
    resistances
end

function laplacian(G::SparseMatrixCSC)
    G = G - spdiagm(diag(G))
    G = -G + spdiagm(vec(sum(G, 1)))
end

function postprocess(volt, cond, p1, p2, resistances)
    fname = "/tmp/voltages_$(p1)_$(p2).txt"

    open(fname, "a") do f
        for i in eachindex(cond)
            write(f, string(p1), '\t', string(cond[i]), '\t', string(volt[cond[i]] - volt[p1]))
            write(f, '\n')
        end
    end

    r = resistances[p1, p2] = resistances[p2, p1] = volt[p2] - volt[p1]
end
