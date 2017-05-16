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

function write_cur_maps(G, voltages, finitegrounds, cc, name, cfg; nodemap = Matrix{Float64}(),
                                                                    hbmeta = RasterMeta())

    node_currents, branch_currents = _create_current_maps(G, voltages, finitegrounds, cfg, nodemap = nodemap, hbmeta = hbmeta)

    if cfg["data_type"] == "network"
        branch_currents_array = _convert_to_3col(branch_currents, cc)
        node_currents_array = _append_name_to_node_currents(node_currents, cc)
        write_currents(node_currents_array, branch_currents_array, name, cfg)
    else
       current_map = node_currents
       write_aagrid(current_map, name, cfg, hbmeta)
   end
end

function write_currents(node_curr_arr, branch_curr_arr, name, cfg)
   pref = split(cfg["output_file"], '.')[1]
   writedlm("$(pref)_node_currents_$(name).txt", node_curr_arr, '\t')
   writedlm("$(pref)_branch_currents_$(name).txt", branch_curr_arr, '\t')
end

_append_name_to_node_currents(node_currents, cc) = [cc node_currents]
   

function _convert_to_3col(branch_currents, cc)

    l = length(branch_currents.nzval)
    graph = zeros(l, 3)
    graph[:,3] = branch_currents.nzval

    # Inspired by show method for sparse matrices in Julia Base
    b = branch_currents

    k = 1
    for i = 1:size(b, 1)
        for j in nzrange(b, i)
            row = b.rowval[j]
            graph[k,1] = cc[row]
            graph[k,2] = cc[i]
            k += 1
        end
    end

    graph
end

function _create_current_maps(G, voltages, finitegrounds, cfg; nodemap = Matrix{Float64}(), hbmeta = RasterMeta())

    node_currents = get_node_currents(G, voltages, finitegrounds)

    if cfg["data_type"] == "network"

        branch_currents = _get_branch_currents(G, voltages, true)
        branch_currents = abs(branch_currents)
        return node_currents, branch_currents

    else

        idx = find(nodemap)
        current_map = zeros(hbmeta.nrows, hbmeta.ncols)
        current_map[idx] = node_currents[Int.(nodemap[idx])]
        return current_map, spzeros(0,0)

    end

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
    dropzeros!(branch_currents)

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
    maxcur = maximum(branch_currents)
    map!(x -> abs(x / maxcur) < 1e-8 ? 0 : x, branch_currents)
    branch_currents
end

function write_aagrid(map, name, cfg, hbmeta; voltage = false, cum = false, max = false)
    pref = split(cfg["output_file"], '.')[1]

    str = "curmap"
    if cum
        str = "cum_$(str)"
    elseif max 
        str = "max_$(str)"
    elseif voltage
        str = "voltmap"
    end

    writedlm("$(pref)_$(str)_$(name).asc", round(map, 8), ' ')
end

function write_volt_maps(name, voltages, cc, nodemap, cfg, hbmeta)

    if cfg["data_type"] == "network"
        write_voltages(cfg["output_file"], name, voltages, cc)
    else
        vm = _create_voltage_map(voltages, nodemap, hbmeta)
        write_aagrid(vm, name, cfg, hbmeta, voltage = true)
    end
end

function write_voltages{T}(output, name, voltages::Vector{T}, cc)

    volt_arr = zeros(T, size(voltages, 1), 2)
    volt_arr[:,1] = cc
    volt_arr[:,2] = voltages
    
    pref = split(output, '.')[1]
    writedlm("$(pref)_voltages_$(name).txt", volt_arr)

end

function _create_voltage_map{T}(voltages::Vector{T}, nodemap, hbmeta)
    voltmap = zeros(T, hbmeta.nrows, hbmeta.ncols)
    idx = find(nodemap)
    voltmap[idx] = voltages[Int.(nodemap[idx])]
    voltmap
end
