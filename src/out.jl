function compute_3col{T}(resistances::Matrix{T}, fp)
    l = length(fp)
    r3col = zeros(T, div(l * (l-1), 2), 3)
    k = 1
    for i = 1:l
        for j = i+1:l
            r3col[k,1] = fp[i]
            r3col[k,2] = fp[j]
            r3col[k,3] = resistances[j,i]
            k += 1
        end
    end
    r3col
end

function write_cur_maps(cond, voltages, finitegrounds, cc, name)

    node_currents, branch_currents = _create_current_maps(cond, voltages, finitegrounds)
    branch_currents_array = _convert_to_3col(branch_currents, cc)
    node_currents_array = _append_name_to_node_currents(node_currents, cc)
   
    write_currents(node_currents_array, branch_currents_array, name)
end

function write_currents(node_curr_arr, branch_curr_arr, name)
   pref = split(cfg["output_file"], '.')[1]
   writedlm("$(pref)_$(name).txt", node_curr_arr, '\t')
end

function _append_name_to_node_currents(node_currents, cc)
    @show node_currents
    @show cc

    [cc node_currents]
end
   

function _convert_to_3col(branch_currents, cc)

    l = length(branch_currents.nzval)
    graph = zeros(l, 3)
    graph[:,3] = branch_currents.nzval
    idx = find(branch_currents)
    I,J,V = findnz(branch_currents)

    for i = 1:size(I, 1)
        graph[i,1] = cc[I[i]]
        graph[i,2] = cc[J[i]]
    end

    graph
end

function _create_current_maps(G, voltages, finitegrounds)

    node_currents = get_node_currents(G, voltages, finitegrounds)
    branch_currents = _get_branch_currents(G, voltages, true)
    branch_currents = abs(branch_currents)

    node_currents, branch_currents
end

function get_node_currents(G, voltages, finitegrounds)

    node_currents_pos = _get_node_currents_posneg(G, voltages, finitegrounds, true)
    node_currents_neg = _get_node_currents_posneg(G, voltages, finitegrounds, false)
    node_currents = map((x,y) -> x > y ? x : y, node_currents_pos, node_currents_neg)
    
end

function _get_node_currents_posneg(G, voltages, finitegrounds, pos)

    branch_currents = _get_branch_currents(G, voltages, pos)
    branch_currents = branch_currents - branch_currents'
    I,J,V = findnz(branch_currents)
    mask = V .> 0
    n = size(G, 1)
    branch_currents = sparse(I[mask], J[mask], V[mask], n, n)
    vec(sum(branch_currents, 1))
end

function _get_branch_currents(G, voltages, pos)

	branch_currents = _get_branch_currents_posneg(G, voltages, pos)
	I,J,V = findnz(G)
    n = size(G, 1)
	mask = I .< J
	branch_currents = sparse(I[mask], J[mask], branch_currents, n, n)

	branch_currents  
end

function _get_branch_currents_posneg{T}(G, v::Vector{T}, pos)
    I,J,V = findnz(G)
    mask = I .< J
    vdiff = zeros(T, sum(mask))
    if pos
        vdiff = v[I[mask]] - v[J[mask]]
        #for (i,v) in enumerate(find(mask))
        #    vdiff[i] = v[I[v]] - v[J[v]]
        #end
    else
        vdiff = v[J[mask]] - v[I[mask]]
        #for (i,v) in enumerate(find(mask))
        #    vdiff[i] = v[J[v]] - v[I[v]]
        #end
    end
    map!(x -> x < 0 ? -x : 0, V)

    branch_currents = vdiff .* V[mask]
end
