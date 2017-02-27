function single_ground_all_pair_resistances{T}(a::SparseMatrixCSC, g::Graph, c::Vector{T})
    numpoints = size(c, 1)
    cc = connected_components(g)
    debug("There are $numpoints nodes and $(length(cc)) connected components")
    resistances = -1 * ones(numpoints, numpoints) 

    cond = laplacian(a)

    volt = Vector{Float64}(size(g, 1))
    total = Int(numpoints * (numpoints-1) / 2)
    
    p = 0 
    for i = 1:numpoints
        rcc = rightcc(cc, c[i])
        cond_pruned = cond[rcc, rcc]
        pt1 = ingraph(rcc, c[i])
        d = cond_pruned[pt1, pt1]
        cond_pruned[pt1, pt1] = 0
        M = aspreconditioner(SmoothedAggregationSolver(cond_pruned))
        for j = i+1:numpoints
            pt2 = ingraph(rcc, c[j])
            if pt2 == 0
                continue
            end
            #info("Solving for pair $p of $total")
            debug("pt1 = $pt1, pt2 = $pt2")
            p +=1
            curr = zeros(size(cond_pruned, 1))
            curr[pt1] = -1
            curr[pt2] = 1
            println('\n')
            volt = cg(cond_pruned, curr, M; tol = 1e-6, maxiter = 100000)
            postprocess(volt[1], c, i, j, resistances, pt1, pt2)
        end
        cond_pruned[pt1,pt1] = d
    end
    debug("solved $p equations")
    for i = 1:size(resistances,1)
        resistances[i,i] = 0
    end
    resistances
end

@inline function rightcc{T}(cc::Vector{Vector{T}}, c::T)
    for i in eachindex(cc)
        if c in cc[i]
            return cc[i]
        end
    end
end

@inline function ingraph{T}(cc::Vector{T}, c::T)
    findfirst(cc, c)
end

function laplacian(G::SparseMatrixCSC)
    G = G - spdiagm(diag(G))
    G = -G + spdiagm(vec(sum(G, 1)))
end

function postprocess(volt, cond, i, j, resistances, pt1, pt2)
    #fname = "/tmp/voltages_$(p1)_$(p2).txt"

    #=open(fname, "a") do f
        for i in eachindex(cond)
            write(f, string(p1), '\t', string(cond[i]), '\t', string(volt[cond[i]] - volt[p1]))
            write(f, '\n')
        end
    end=#

    r = resistances[i, j] = resistances[j, i] = volt[pt2] - volt[pt1]
end

function compute_network(a::Inifile)
    network_file = get(a, "Habitat raster or graph", "habitat_file")
    point_file = get(a, "Options for pairwise and one-to-all and all-to-one modes",
                        "point_file")
    A = read_graph(a, network_file)
    g = Graph(A)
    scenario = get(a, "Circuitscape mode", "scenario")
    if scenario == "pairwise"
        fp = read_focal_points(point_file)
        resistances = single_ground_all_pair_resistances(A, g, fp)
        return resistances
    elseif scenario == "advanced"
        # advanced logic
    end
end
