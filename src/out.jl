struct OutputFlags
    write_volt_maps::Bool
    write_cur_maps::Bool
    write_cum_cur_maps_only::Bool
    write_max_cur_maps::Bool
    set_null_currents_to_nodata::Bool
    set_null_voltages_to_nodata::Bool
    compress_grids::Bool
    log_transform_maps::Bool
end

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

    
function write_cur_maps(name, output, component_data, finitegrounds, flags, cfg)
    
    # Get desired data
    G = component_data.matrix
    voltages = flags.is_advanced ? output : output.voltages
    cc = component_data.cc
    nodemap = component_data.local_nodemap
    hbmeta = component_data.hbmeta

    node_currents, branch_currents = _create_current_maps(G, voltages, finitegrounds, cfg, nodemap = nodemap, hbmeta = hbmeta)

    if !flags.is_raster
        branch_currents_array = _convert_to_3col(branch_currents, cc)
        node_currents_array = _append_name_to_node_currents(node_currents, cc)
        write_currents(node_currents_array, branch_currents_array, name, cfg)
    else
       current_map = node_currents
       # log_transform = cfg["log_transform_maps"] == "True"
       log_transform = flags.outputflags.log_transform_maps
       write_aagrid(current_map, name, cfg, hbmeta; log_transform = log_transform)
   end
end

function write_currents(node_curr_arr, branch_curr_arr, name, cfg)
   pref = split(cfg["output_file"], '.')[1]
   writedlm("$(pref)_node_currents$(name).txt", node_curr_arr, '\t')
   writedlm("$(pref)_branch_currents$(name).txt", branch_curr_arr, '\t')
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
        branch_currents = abs.(branch_currents)
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
    branch_currents = sparse( I[mask], J[mask], V[mask], n, n)

	if finitegrounds[1]!= -9999
        finiteground_currents = finitegrounds .* voltages
        if pos
            map!(x -> x < 0 ? -x : 0, finiteground_currents, finiteground_currents)
        else
            map!(x -> x > 0 ? x : 0, finiteground_currents, finiteground_currents)
        end
        n = size(G, 1)
        branch_currents = branch_currents + spdiagm(finiteground_currents, 0, n, n)
    end

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
    map!(x -> x < 0 ? -x : 0, V, V)

    branch_currents = vdiff .* V[mask]
    maxcur = maximum(branch_currents)
    map!(x -> abs(x / maxcur) < 1e-8 ? 0 : x, branch_currents, branch_currents)
    branch_currents
end

function write_aagrid(cmap, name, cfg, hbmeta;
                        voltage = false, cum = false, max = false,
                        log_transform = false)

    pref = split(cfg["output_file"], '.')[1]

    if log_transform
        map!(x -> x > 0 ? log10(x) : float(hbmeta.nodata), cmap, cmap)
    end

    str = "curmap"
    if cum
        str = "cum_$(str)"
    elseif max
        str = "max_$(str)"
    elseif voltage
        str = "voltmap"
    end

    filename = "$(pref)_$(str)$(name).asc"
    f = open(filename, "w")

    write(f, "ncols         $(hbmeta.ncols)\n")
    write(f, "nrows         $(hbmeta.nrows)\n")
    write(f, "xllcorner     $(hbmeta.xllcorner)\n")
    write(f, "yllcorner     $(hbmeta.yllcorner)\n")
    write(f, "cellsize      $(hbmeta.cellsize)\n")
    write(f, "NODATA_value  $(hbmeta.nodata)\n")

    writedlm(f, round.(cmap, 8), ' ')
    close(f)
end

function write_volt_maps(name, output, component_data, flags, cfg)

    voltages = flags.is_advanced ? output : output.voltages

    if !flags.is_raster

        cc = component_data.cc
        write_voltages(cfg["output_file"], name, voltages, cc)

    else

        # Desired data
        # voltages = output.voltages
        cc = component_data.cc
        hbmeta = component_data.hbmeta
        nodemap = component_data.local_nodemap

        vm = _create_voltage_map(voltages, nodemap, hbmeta)
        write_aagrid(vm, name, cfg, hbmeta, voltage = true)
    end
end

function write_voltages{T}(output, name, voltages::Vector{T}, cc)

    volt_arr = zeros(T, size(voltages, 1), 2)
    volt_arr[:,1] = cc
    volt_arr[:,2] = voltages

    pref = split(output, '.')[1]
    writedlm("$(pref)_voltages$(name).txt", volt_arr)

end

function _create_voltage_map{T}(voltages::Vector{T}, nodemap, hbmeta)
    voltmap = zeros(T, hbmeta.nrows, hbmeta.ncols)
    idx = find(nodemap)
    voltmap[idx] = voltages[Int.(nodemap[idx])]
    voltmap
end

alloc_map(hbmeta) = zeros(hbmeta.nrows, hbmeta.ncols)

function accum_voltages!(base, newvolt, nodemap, hbmeta)
    voltmap = _create_voltage_map(newvolt, nodemap, hbmeta)
    for i in eachindex(base)
        base[i] += voltmap[i]
    end
end

function accum_currents!(base, newcurr, cfg, G, voltages, finitegrounds, nodemap, hbmeta)
    node_currents, branch_currents = _create_current_maps(G, voltages, finitegrounds, cfg,
                                                            nodemap = nodemap, hbmeta = hbmeta)

    for i in eachindex(base)
        base[i] += node_currents[i]
    end
end

