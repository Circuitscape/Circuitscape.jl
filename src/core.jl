struct GraphData{T,V}
    G::SparseMatrixCSC{T,V}
    cc::Vector{Vector{V}}
    points::Vector{V}
    user_points::Vector{V}
    exclude_pairs::Vector{Tuple{V,V}}
    nodemap::Matrix{V}
    polymap::Matrix{V}
end

"""
Core kernel of Circuitscape - used to solve several pairs

Input:
* data::GraphData
"""
function single_ground_all_pairs(data, flags, cfg)

    # Data
    a = data.G
    cc = data.cc
    points = data.points
    exclude = data.exclude_pairs
    nodemap = data.nodemap
    polymap = data.polymap
    orig_pts = data.user_points
    hbmeta = RasterMeta()

    # Flags
    outputflags = flags.outputflags
    is_raster = flags.is_raster
    write_volt_maps = outputflags.write_volt_maps
    write_cur_maps = outputflags.write_cur_maps

    #pairs = get_num_pairs(cc, points, exclude_pairs)
    #info("Total number of pair solves = $(length(pairs))")

    # Get number of focal points
    numpoints = size(points, 1)
    
    info("Graph has $(size(a,1)) nodes, $numpoints focal points and $(length(cc)) connected components")
    
    # Initialize pairwise resistance
    resistances = -1 * ones(eltype(a), numpoints, numpoints)
    voltmatrix = zeros(eltype(a), size(resistances))
    
    # Take laplacian of matrix
    # t = @elapsed lap = laplacian(a)
    # info("Time to construct laplacian = $t seconds")
    
    # Get a vector of connected components
    comps = getindex.([a], cc, cc)
    
    get_shortcut_resistances = false
    if is_raster && !write_volt_maps && !write_cur_maps
        get_shortcut_resistances = true
        shortcut = -1 * ones(size(resistances))
        covered = Dict{Int, Bool}()
        for i = 1:size(cc, 1)
            covered[i] = false
        end
    end
    
    for (cid, comp) in enumerate(cc)
    
        # Subset of points relevant to CC
        csub = filter(x -> x in comp, points) |> unique
        #idx = findin(c, csub)
    
        if isempty(csub)
            continue
        end

        # Conductance matrix corresponding to CC
        matrix = comps[cid]
    
        # Construct preconditioner *once* for every CC
        t1 = @elapsed P = aspreconditioner(smoothed_aggregation(matrix))
        info("Time taken to construct preconditioner = $t1 seconds")

        # Get local nodemap for CC - useful for output writing
        t2 = @elapsed local_nodemap = construct_local_node_map(nodemap, comp, polymap)
        info("Time taken to construct local nodemap = $t2 seconds")
    
        function f(i)

            # Generate return type
            ret = Vector{Tuple{Int,Int,Float64}}()

            pi = csub[i]
            comp_i = findfirst(comp, pi)
            I = find(x -> x == pi, points)
            smash_repeats!(ret, I)

            # Preprocess matrix
            d = matrix[comp_i, comp_i]

            # Iteration space through all possible pairs
            rng = i+1:size(csub, 1)

            # Loop through all possible pairs
            for j in rng

                pj = csub[j]
                comp_j = findfirst(comp, pj)
                J = find(x -> x == pj, points)

                # Forget excluded pairs
                ex = false
                for c_i in I, c_j in J
                    if (c_i, c_j) in exclude
                        ex = true
                        break
                    end
                end
                ex && continue

                if pi == pj
                    continue
                end

                # Initialize currents
                current = zeros(size(matrix, 1))
                current[comp_i] = -1
                current[comp_j] = 1

                # Solve system
                info("Solving points $pi and $pj")
                t2 = @elapsed v = solve_linear_system(cfg, matrix, current, P)
                info("Time taken to solve linear system = $t2 seconds")
                v .= v .- v[comp_i]

                # Calculate resistance
                r = v[comp_j] - v[comp_i]

                # Return resistance value
                for c_i in I, c_j in J
                    push!(ret, (c_i, c_j, r))
                    postprocess(v, points, c_i, c_j, r, comp_i, comp_j, matrix, comp, cfg, voltmatrix,
                                                get_shortcut_resistances;
                                                local_nodemap = local_nodemap,
                                                orig_pts = orig_pts,
                                                polymap = polymap,
                                                hbmeta = hbmeta)
                end
            end

        matrix[comp_i, comp_i] = d

        ret
        end

        X = pmap(x ->f(x), 1:size(csub,1))

        # Set all resistances
        for x in X
            for i = 1:size(x, 1)
                resistances[x[i][1], x[i][2]] = x[i][3]
                resistances[x[i][2], x[i][1]] = x[i][3]
            end
        end

    end

    for i = 1:size(resistances,1)
        resistances[i,i] = 0
    end

    # Pad it with the user points
    vcat(vcat(0,orig_pts)', hcat(orig_pts, resistances))

end

"""
Returns all possible pairs to solve.

Input:
* ccs::Vector{Vector{Int}} - vector of connected components
* fp::Vector{Int} - vector of focal points 
* exclude_pairs::Vector{Tuple{Int,Int}} - vector of point pairs (tuples) to exclude 

Output: 
* n - total number of pairs
"""
function get_num_pairs(ccs, fp, exclude_pairs)

    num = 0

    for (i,cc) in enumerate(ccs)
        sub_fp = findin(x -> x in cc, fp)
        l = endof(sub_fp)
        for ii = 1:l
            pt1 = sub_fp[ii]
            for jj = ii+1:l
                pt2 = sub_fp[jj]
                if (pt1, pt2) in exclude_pairs
                    continue
                else
                    num += 1
                end
            end
        end
    end
    num
end

function smash_repeats!(ret, I)
    for i = 1:size(I,1)
        for j = i+1:size(I,1)
            push!(ret, (I[i], I[j], 0))
        end
    end
end

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

function solve_linear_system(v, G, curr, M)
    cg(G, curr, Pl = M, tol = 1e-6, maxiter = 100000)
end
#solve_linear_system(cfg, G, curr, M) = solve_linear_system!(cfg, zeros(size(curr)), G, curr, M)

function postprocess(volt, cond, i, j, r, pt1, pt2, cond_pruned, cc, cfg, voltmatrix,
                                                                get_shortcut_resistances;
                                                                local_nodemap = Matrix{Float64}(),
                                                                orig_pts = Vector{Int}(),
                                                                polymap = NoPoly(),
                                                                hbmeta = hbmeta)

    if get_shortcut_resistances
        update_voltmatrix!(voltmatrix, local_nodemap, volt, hbmeta, cond, r, j, cc)
        return nothing
    end

    name = "_$(cond[i])_$(cond[j])"

    if cfg["data_type"] == "raster"
        name = "_$(Int(orig_pts[i]))_$(Int(orig_pts[j]))"
    end

    if cfg["write_volt_maps"] == "True"
        t = @elapsed write_volt_maps(name, volt, cc, cfg, hbmeta = hbmeta, nodemap = local_nodemap)
        info("Time taken to write voltage maps = $t seconds")
    end

    if cfg["write_cur_maps"] == "True"
        t = @elapsed write_cur_maps(cond_pruned, volt, [-9999.], cc, name, cfg;
        nodemap = local_nodemap,
        hbmeta = hbmeta)
        info("Time taken to write current maps = $t seconds")
    end
    nothing
end