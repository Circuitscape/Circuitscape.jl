struct Cumulative{T}
    cum_curr::Vector{SharedMatrix{T}}
    max_curr::Vector{SharedMatrix{T}}
    cum_branch_curr::Vector{SharedVector{T}}
    cum_node_curr::Vector{SharedVector{T}}
end

struct GraphData{T,V}
    G::SparseMatrixCSC{T,V}
    cc::Vector{Vector{V}}
    points::Vector{V}
    user_points::Vector{V}
    exclude_pairs::Vector{Tuple{V,V}}
    nodemap::Matrix{V}
    polymap::Matrix{V}
    hbmeta::RasterMeta
    cellmap::Matrix{T}
    cum::Cumulative{T}
end

struct ComponentData{T,V}
    cc::Vector{V}
    matrix::SparseMatrixCSC{T,V}
    local_nodemap::Matrix{V}
    hbmeta::RasterMeta
    cellmap::Matrix{T}
end

struct Output{T,V}
    points::Vector{V}
    voltages::Vector{T}
    orig_pts::Tuple{V,V}
    comp_idx::Tuple{V,V}
    resistance::T
    col::V
    cum::Cumulative{T}
end

struct Shortcut{T}
    get_shortcut_resistances::Bool
    voltmatrix::Matrix{T}
    shortcut_res::Matrix{T}
end

"""
Core kernel of Circuitscape - used to solve several pairs

Input:
* data::GraphData
"""
function single_ground_all_pairs(data::GraphData{T,V}, flags, cfg, log = true) where {T,V}

    if flags.solver in AMG
        csinfo("Solver used: AMG accelerated by CG")
        amg_solver_path(data, flags, cfg, log)
    else
        csinfo("Solver used: CHOLMOD")
        bs = parse(Int, cfg["cholmod_batch_size"])
        _cholmod_solver_path(data, flags, cfg, log, bs)
    end
end

function amg_solver_path(data::GraphData{T,V}, flags, cfg, log)::Matrix{T} where {T,V}

    # Data
    a = data.G
    cc = data.cc
    points = data.points
    exclude = data.exclude_pairs
    nodemap = data.nodemap
    polymap = data.polymap
    orig_pts = data.user_points
    hbmeta = data.hbmeta
    cellmap = data.cellmap

    # Flags
    outputflags = flags.outputflags
    is_raster = flags.is_raster
    write_volt_maps = outputflags.write_volt_maps
    write_cur_maps = outputflags.write_cur_maps
    write_cum_cur_map_only = outputflags.write_cum_cur_map_only
    write_max_cur_maps = outputflags.write_max_cur_maps

    # Get number of focal points
    numpoints = size(points, 1)

    # Cumulative currents
    cum = data.cum
    
    csinfo("Graph has $(size(a,1)) nodes, $numpoints focal points and $(length(cc)) connected components")

    num, d = get_num_pairs(cc, points, exclude)
    log && csinfo("Total number of pair solves = $num")
    
    # Initialize pairwise resistance
    resistances = -1 * ones(T, numpoints, numpoints)::Matrix{T}
    voltmatrix = zeros(T, size(resistances))::Matrix{T}
    shortcut_res = deepcopy(resistances)::Matrix{T}
    
    # Get a vector of connected components
    comps = getindex.([a], cc, cc)
    
    get_shortcut_resistances = false
    if is_raster && !write_volt_maps && !write_cur_maps && 
            !write_cum_cur_map_only && !write_max_cur_maps && 
            isempty(exclude)
        get_shortcut_resistances = true
        csinfo("Triggering resistance calculation shortcut")
        num, d = get_num_pairs_shortcut(cc, points, exclude)
        csinfo("Total number of pair solves has been reduced to $num ")
    end
    shortcut = Shortcut(get_shortcut_resistances, voltmatrix, shortcut_res)    
  
    for (cid, comp) in enumerate(cc)
    
        # Subset of points relevant to CC
        csub = filter(x -> x in comp, points) |> unique
        #idx = findin(c, csub)
    
        if isempty(csub)
            continue
        end

        # Conductance matrix corresponding to CC
        matrix = comps[cid]        

        # Regularization step
        matrix.nzval .+= eps(eltype(matrix)) * norm(matrix.nzval)

        # Construct preconditioner *once* for every CC
        t1 = @elapsed P = aspreconditioner(smoothed_aggregation(matrix))
        csinfo("Time taken to construct preconditioner = $t1 seconds")

        # Get local nodemap for CC - useful for output writing
        t2 = @elapsed local_nodemap = construct_local_node_map(nodemap, comp, polymap)
        csinfo("Time taken to construct local nodemap = $t2 seconds")

        component_data = ComponentData(comp, matrix, local_nodemap, hbmeta, cellmap)        

        function f(i)

            # Generate return type
            ret = Vector{Tuple{V,V,T}}()

            pi = csub[i]
            comp_i = something(findfirst(isequal(pi),comp), 0)
            comp_i = V(comp_i)
            I = findall(x -> x == pi, points)
            smash_repeats!(ret, I)

            # Preprocess matrix
            # d = matrix[comp_i, comp_i]

            # Iteration space through all possible pairs
            rng = i+1:size(csub, 1)
            if nprocs() > 1 
                for j in rng
                    pj = csub[j]
                    csinfo("Scheduling pair $(d[(pi,pj)]) of $num to be solved")
                end
            end

            # Loop through all possible pairs
            for j in rng

                pj = csub[j]
                comp_j = something(findfirst(isequal(pj), comp), 0)
                comp_j = V(comp_j)
                J = findall(x -> x == pj, points)

                if pi == pj
                    continue
                end

                # Forget excluded pairs
                ex = false
                for c_i in I
                    for c_j in J
                        if (c_i, c_j) in exclude
                            continue
                        end

                        # Initialize currents
                        current = zeros(T, size(matrix, 1))
                        current[comp_i] = -1
                        current[comp_j] = 1

                        # Solve system
                        # csinfo("Solving points $pi and $pj")
                        log && csinfo("Solving pair $(d[(pi,pj)]) of $num")
                        t2 = @elapsed v = solve_linear_system(cfg, matrix, current, P)
                        csinfo("Time taken to solve linear system = $t2 seconds")
                        v .= v .- v[comp_i]

                        # Calculate resistance
                        r = v[comp_j] - v[comp_i]

                        # Return resistance value
                        push!(ret, (c_i, c_j, r))
                        if get_shortcut_resistances
                            resistances[c_i, c_j] = r
                            resistances[c_j, c_i] = r
                        end
                        output = Output(points, v, (orig_pts[c_i], orig_pts[c_j]),
                                        (comp_i, comp_j), r, V(c_j), cum)
                        postprocess(output, component_data, flags, shortcut, cfg)
                    end
                end
            end

        # matrix[comp_i, comp_i] = d

        ret
        end

        if get_shortcut_resistances        
            idx = something(findfirst(isequal(csub[1]), points), 0)
            f(1)
            update_shortcut_resistances!(idx, shortcut, resistances, points, comp)
        else
            X = pmap(x ->f(x), 1:size(csub,1))

            # Set all resistances
            for x in X
                for i = 1:size(x, 1)
                    resistances[x[i][1], x[i][2]] = x[i][3]
                    resistances[x[i][2], x[i][1]] = x[i][3]
                end
            end
        end

    end

    if get_shortcut_resistances
        resistances = shortcut.shortcut_res
    end

    for i = 1:size(resistances,1)
        resistances[i,i] = 0
    end

    # Pad it with the user points
    r = vcat(vcat(0,orig_pts)', hcat(orig_pts, resistances))

    # Save resistances
    save_resistances(r, cfg)

    r
end

struct CholmodNode{T}
    cc_idx::Tuple{T,T}
    points_idx::Tuple{T,T}
end

function _cholmod_solver_path(data::GraphData{T,V}, flags, 
                                  cfg, log, batch_size = 1000) where {T,V}
    
    # Data
    a = data.G
    cc = data.cc
    points = data.points
    exclude = data.exclude_pairs
    nodemap = data.nodemap
    polymap = data.polymap
    orig_pts = data.user_points
    hbmeta = data.hbmeta
    cellmap = data.cellmap

    # Flags
    outputflags = flags.outputflags
    is_raster = flags.is_raster
    write_volt_maps = outputflags.write_volt_maps
    write_cur_maps = outputflags.write_cur_maps
    write_cum_cur_map_only = outputflags.write_cum_cur_map_only
    write_max_cur_maps = outputflags.write_max_cur_maps

    # Cumulative current map
    cum = data.cum

    # CHOLMOD solver mode works only in double precision
    if eltype(a) == Float32
        cswarn("Converting single precision matrix to double")
        a = Float64.(a)
    end

    # Get number of focal points
    numpoints = size(points, 1)
    
    csinfo("Graph has $(size(a,1)) nodes, $numpoints focal points and $(length(cc)) connected components")

    num, d = get_num_pairs(cc, points, exclude)
    log && csinfo("Total number of pair solves = $num")
    
    # Initialize pairwise resistance
    resistances = -1 * ones(eltype(a), numpoints, numpoints)
    voltmatrix = zeros(eltype(a), size(resistances))
    shortcut_res = -1 * ones(eltype(a), size(resistances))
    
    # Get a vector of connected components
    comps = getindex.([a], cc, cc)
    
    get_shortcut_resistances = false
    if is_raster && !write_volt_maps && !write_cur_maps && 
            !write_cum_cur_map_only  && !write_max_cur_maps &&
            isempty(exclude)
        get_shortcut_resistances = true
        csinfo("Triggering resistance calculation shortcut")
        num, d = get_num_pairs_shortcut(cc, points, exclude)
        csinfo("Total number of pair solves has been reduced to $num ")
    end
    shortcut = Shortcut(get_shortcut_resistances, voltmatrix, shortcut_res)
    
    for (cid, comp) in enumerate(cc)
    
        # Subset of points relevant to CC
        csub = filter(x -> x in comp, points) |> unique
        #idx = findin(c, csub)
    
        if isempty(csub)
            continue
        end

        # Conductance matrix corresponding to CC
        matrix = comps[cid]

        t = @elapsed factor = construct_cholesky_factor(matrix)
        csinfo("Time taken to construct cholesky factor = $t")

        # Get local nodemap for CC - useful for output writing
        t2 = @elapsed local_nodemap = construct_local_node_map(nodemap, comp, polymap)
        csinfo("Time taken to construct local nodemap = $t2 seconds")

        component_data = ComponentData(comp, matrix, local_nodemap, hbmeta, cellmap)

        ret = Vector{Tuple{V,V,Float64}}()

        cholmod_batch = CholmodNode[]
        
        # Batched backsubstitution
        function g(i)

            pi = csub[i]
            comp_i = V(something(findfirst(isequal(pi),comp),0))
            I = findall(x -> x == pi, points)
            # smash_repeats!(ret, I)
            smash_repeats!(resistances, I)

            # Iteration space through all possible pairs
            rng = i+1:size(csub, 1)

            # Loop through all possible pairs
            for j in rng

                pj = csub[j]
                comp_j = V(something(findfirst(isequal(pj), comp),0))
                J = findall(x -> x == pj, points)

                if pi == pj
                    continue
                end

                # Forget excluded pairs
                for c_i in I, c_j in J
                    if (c_i, c_j) in exclude
                        continue
                    else
                        push!(cholmod_batch, 
                          CholmodNode((comp_i, comp_j), (V(c_i), V(c_j))))
                    end
                end
            end
        end

        function f(i, rng, lhs)
            v = rng[i]
            output = Output(points, lhs[:,i], 
                (orig_pts[cholmod_batch[v].points_idx[1]], 
                orig_pts[cholmod_batch[v].points_idx[2]]), 
                cholmod_batch[v].cc_idx, 
                lhs[cholmod_batch[v].cc_idx[2], i] - lhs[cholmod_batch[v].cc_idx[1], i],
                V(cholmod_batch[v].points_idx[2]), cum)
            postprocess(output, component_data, flags, shortcut, cfg)
        end
        if get_shortcut_resistances
            idx = something(findfirst(isequal(csub[1]), points),0)
            g(1)
        else 
            g.(1:size(csub, 1))
        end

        # l = min(num, length(cholmod_batch))
        l = length(cholmod_batch)

        for st in 1:batch_size:l

            rng = st + batch_size <= l ?
                            (st:(st+batch_size-1)) : (st:l)

            csinfo("Solving points $(rng.start) to $(rng.stop)")

            rhs = zeros(eltype(matrix), size(matrix, 1), length(rng))

            for (i,v) in enumerate(rng)
                node = cholmod_batch[v]
                rhs[node.cc_idx[1], i] = -1
                rhs[node.cc_idx[2], i] = 1 
            end
            lhs = factor \ rhs

            # Normalisation step
            for (i,val) in enumerate(rng)
                n = cholmod_batch[val].cc_idx[1]
                v = lhs[n,i]
                for j = 1:size(matrix, 1)
                    lhs[j,i] = lhs[j,i] - v
                end
            end


            pmap(x -> f(x, rng, lhs), 1:length(rng))
            for (i,v) in enumerate(rng)
                coords = cholmod_batch[v].points_idx
                r = lhs[cholmod_batch[v].cc_idx[2], i] - 
                            lhs[cholmod_batch[v].cc_idx[1], i]
                resistances[coords...] = r
                resistances[reverse(coords)...] = r
            end
        end

        if get_shortcut_resistances
            update_shortcut_resistances!(idx, shortcut, resistances, points, comp)
        end
    end

    if get_shortcut_resistances
        resistances = shortcut.shortcut_res
    end

    for i = 1:size(resistances,1)
        resistances[i,i] = 0
    end

    # Pad it with the user points
    r = vcat(vcat(0,orig_pts)', hcat(orig_pts, resistances))

    # Save resistances
    save_resistances(r, cfg)

    r
end

function construct_cholesky_factor(matrix)
    cholesky(matrix + sparse(10eps()*I,size(matrix)...))
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
function get_num_pairs(ccs, fp::Vector{V}, exclude_pairs) where V

    num = 0
    d = Dict{Tuple{V,V}, V}()

    for (i,cc) in enumerate(ccs)
        sub_fp = filter(x -> x in cc, fp) |> unique
        l = lastindex(sub_fp)
        for ii = 1:l
            pt1 = sub_fp[ii]
            for jj = ii+1:l
                pt2 = sub_fp[jj]
                if (pt1, pt2) in exclude_pairs
                    continue
                else
                    num += 1
                    d[(pt1, pt2)] = num
                end
            end
        end
    end
    num, d
end

function get_num_pairs_shortcut(ccs, fp::Vector{V}, exclude_pairs) where V

    num = 0
    d = Dict{Tuple{V,V}, V}()

    for (i,cc) in enumerate(ccs)
        sub_fp = filter(x -> x in cc, fp) |> unique
        l = lastindex(sub_fp)
        l == 0 && continue
        for ii = 1:1
            pt1 = sub_fp[ii]
            for jj = ii+1:l
                pt2 = sub_fp[jj]
                if (pt1, pt2) in exclude_pairs
                    continue
                else
                    num += 1
                    d[(pt1, pt2)] = num
                end
            end
        end
    end
    num, d
end
function smash_repeats!(ret, I)
    for i = 1:size(I,1)
        for j = i+1:size(I,1)
            push!(ret, (I[i], I[j], 0))
        end
    end
end

function smash_repeats!(resistances::Matrix{T}, I) where T
    for i = 1:size(I,1)
        for j = i+1:size(I,1)
            resistances[I[i], I[j]] = 0
            resistances[I[j], I[i]] = 0
        end
    end
end

"""
Calculate laplacian of the adjacency matrix of a graph
"""
function laplacian(G::SparseMatrixCSC{T,V}) where {T,V}
    n = size(G, 1)
    s = Vector{eltype(G)}(undef,n)
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
    r = V(1):V(n)
    S = sparse(r, r, s)
    G + S
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

function solve_linear_system(cfg, 
            G::SparseMatrixCSC{T,V}, 
            curr::Vector{T}, M)::Vector{T} where {T,V} 
    cg(G, curr, Pl = M, tol = T(1e-6), maxiter = 100_000)
end

function postprocess(output, component_data, flags, shortcut, cfg)

    
    voltages = output.voltages
    matrix = component_data.matrix
    local_nodemap = component_data.local_nodemap
    hbmeta = component_data.hbmeta
    orig_pts = output.orig_pts

    # Shortcut flags and data
    get_shortcut_resistances = shortcut.get_shortcut_resistances

    if get_shortcut_resistances
        update_voltmatrix!(shortcut, output, component_data)
        return nothing
    end


    name = "_$(orig_pts[1])_$(orig_pts[2])"

    if flags.outputflags.write_volt_maps
        t = @elapsed write_volt_maps(name, output, component_data, flags, cfg)
        csinfo("Time taken to write voltage maps = $t seconds")
    end

    if flags.outputflags.write_cur_maps
        t = @elapsed write_cur_maps(name, output, component_data, 
                                    [-9999.], flags, cfg)
        csinfo("Time taken to write current maps = $t seconds")
    end
    nothing
end

function update_voltmatrix!(shortcut, output, component_data)

    # Data
    voltmatrix = shortcut.voltmatrix
    c = output.points
    cc = component_data.cc
    voltages = output.voltages
    r = output.resistance
    j = output.col

    for i = 2:size(c, 1)
        ind = something(findfirst(isequal(c[i]), cc),0)
        if ind != 0
            voltageAtPoint = voltages[ind]
            voltageAtPoint = 1 - (voltageAtPoint/r)
            voltmatrix[i,j] = voltageAtPoint
        end
    end
end

                            
function update_shortcut_resistances!(anchor, sc, resistances, points, comp)

    # Data
    voltmatrix = sc.voltmatrix
    shortcut = sc.shortcut_res

    check = map(x -> x in comp, points)
    l = size(resistances, 1)
    for pointx = 1:l
        if check[pointx]
            R1x = resistances[anchor, pointx]
            if R1x != -1
                shortcut[pointx, anchor] = shortcut[anchor, pointx] = R1x
                for point2 = pointx:l
                    if check[point2]
                        R12 = resistances[anchor, point2]
                        if R12 != -1
                            if R1x != -777
                                shortcut[anchor, point2] = shortcut[point2, anchor] = R12
                                Vx = voltmatrix[pointx, point2]
                                R2x = 2*R12*Vx + R1x - R12
                                if shortcut[point2, pointx] != -777
                                    shortcut[point2, pointx] = shortcut[pointx, point2] = R2x
                                end
                            else
                                shortcut[pointx, :] = shortcut[:, pointx] = -777
                            end
                        end
                    end
                end
            end
        end
    end
end
