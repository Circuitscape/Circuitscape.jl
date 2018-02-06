struct GraphData{T,V}
    G::SparseMatrixCSC{T,V}
    cc::Vector{Vector{V}}
    points::Vector{V}
    user_points::Vector{V}
    exclude_pairs::Vector{Tuple{V,V}}
    nodemap::Matrix{V}
    polymap::Matrix{V}
end

struct ComponentData{T,V}
    cc::Vector{V}
    matrix::SparseMatrixCSC{T,V}
    local_nodemap::Matrix{V}
    hbmeta::RasterMeta
end

struct Output{T,V}
    points::Vector{V}
    voltages::Vector{T}
    orig_pts::Tuple{V,V}
    comp_idx::Tuple{V,V}
    resistance::T
    col::V
end

struct Shortcut{T}
    get_shortcut_resistances::Bool
    voltmatrix::Matrix{T}
end

"""
Core kernel of Circuitscape - used to solve several pairs

Input:
* data::GraphData
"""
function single_ground_all_pairs(data, flags, cfg)

    if flags.solver in AMG
        info("Solver used: AMG accelerated by CG")
        amg_solver_path(data, flags, cfg)
    else
        info("Solver used: CHOLMOD")
        cholmod_solver_path(data, flags, cfg)
    end
end

function amg_solver_path(data, flags, cfg)

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

    # Get number of focal points
    numpoints = size(points, 1)
    
    info("Graph has $(size(a,1)) nodes, $numpoints focal points and $(length(cc)) connected components")

    num = get_num_pairs(cc, points, exclude)
    info("Total number of pair solves = $num")
    
    # Initialize pairwise resistance
    resistances = -1 * ones(eltype(a), numpoints, numpoints)
    voltmatrix = zeros(eltype(a), size(resistances))
    
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
    shortcut = Shortcut(get_shortcut_resistances, voltmatrix)    
  
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

        component_data = ComponentData(comp, matrix, local_nodemap, hbmeta)        

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
                current = zeros(eltype(a), size(matrix, 1))
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
                    output = Output(points, v, (orig_pts[c_i], orig_pts[c_j]),
                                        (comp_i, comp_j), r, c_j)
                    postprocess(output, component_data, flags, shortcut, cfg)
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

function cholmod_solver_path(data, flags, cfg)

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

    # Get number of focal points
    numpoints = size(points, 1)
    
    info("Graph has $(size(a,1)) nodes, $numpoints focal points and $(length(cc)) connected components")

    num = get_num_pairs(cc, points, exclude)
    info("Total number of pair solves = $num")
    
    # Initialize pairwise resistance
    resistances = -1 * ones(eltype(a), numpoints, numpoints)
    voltmatrix = zeros(eltype(a), size(resistances))
    
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
    shortcut = Shortcut(get_shortcut_resistances, voltmatrix)
  
    for (cid, comp) in enumerate(cc)
    
        # Subset of points relevant to CC
        csub = filter(x -> x in comp, points) |> unique
        #idx = findin(c, csub)
    
        if isempty(csub)
            continue
        end

        # Conductance matrix corresponding to CC
        matrix = comps[cid]

        # Check if positive definite (laplacians may be semi definite too)
        posdef = isposdef(matrix)

        if posdef
            t = @elapsed factor = cholfact(matrix)
        else
            t = @elapsed factor = cholfact(matrix + speye(size(matrix,1))/10^6)
        end
        info("Time taken to construct cholesky factor = $t")

        # Get local nodemap for CC - useful for output writing
        t2 = @elapsed local_nodemap = construct_local_node_map(nodemap, comp, polymap)
        info("Time taken to construct local nodemap = $t2 seconds")

        component_data = ComponentData(comp, matrix, local_nodemap, hbmeta)

        ret = Vector{Tuple{Int,Int,Float64}}()
        
        # Batched backsubstitution
        for i = 1:size(csub,1)

            pi = csub[i]
            comp_i = findfirst(comp, pi)
            I = find(x -> x == pi, points)
            smash_repeats!(ret, I)

            # Generate return type
            
            pi = csub[i]
            comp_i = findfirst(comp, pi)
            I = find(x -> x == pi, points)
            smash_repeats!(ret, I)

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
                current = zeros(eltype(a), size(matrix, 1))
                current[comp_i] = -1
                current[comp_j] = 1

                # Solve system
                info("Solving points $pi and $pj")
                t2 = @elapsed v = factor \ current
                info("Time taken to solve linear system = $t2 seconds")
                v .= v .- v[comp_i]

                # Calculate resistance
                r = v[comp_j] - v[comp_i]

                # Return resistance value
                for c_i in I, c_j in J
                    push!(ret, (c_i, c_j, r))
                    output = Output(points, v, (orig_pts[c_i], orig_pts[c_j]),
                                        (comp_i, comp_j), r, c_j)
                    postprocess(output, component_data, flags, shortcut, cfg)
                end
            end
        end

        # Set all resistances
        for i = 1:size(ret, 1)
            resistances[ret[i][1], ret[i][2]] = ret[i][3]
            resistances[ret[i][2], ret[i][1]] = ret[i][3]
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
        sub_fp = filter(x -> x in cc, fp) |> unique
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

solve_linear_system(cfg, 
            G::SparseMatrixCSC{T,V}, 
            curr::Vector{T}, M) where {T,V} = 
            cg(G, curr, Pl = M, tol = T(1e-6), maxiter = 100_000)

function postprocess(output, component_data, flags, shortcut, cfg)

    
    voltages = output.voltages
    matrix = component_data.matrix
    local_nodemap = component_data.local_nodemap
    hbmeta = component_data.hbmeta
    orig_pts = output.orig_pts

    # Shortcut flags and data
    get_shortcut_resistances = shortcut.get_shortcut_resistances

    if get_shortcut_resistances
        # update_voltmatrix!(voltmatrix, local_nodemap, voltages, hbmeta, points, r, j, cc)
        update_voltmatrix!(shortcut, output, component_data)
        return nothing
    end

    # name = "_$(points[i])_$(points[j])"

    #if flags.is_raster
        # name = "_$(Int(orig_pts[i]))_$(Int(orig_pts[j]))"
    name = "_$(orig_pts[1])_$(orig_pts[2])"
    #end

    if flags.outputflags.write_volt_maps
        t = @elapsed write_volt_maps(name, output, component_data, flags, cfg)
        info("Time taken to write voltage maps = $t seconds")
    end

    if flags.outputflags.write_cur_maps
        t = @elapsed write_cur_maps(name, output, component_data, 
                                    [-9999.], flags, cfg)
        info("Time taken to write current maps = $t seconds")
    end
    nothing
end