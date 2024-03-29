struct OutputFlags
    write_volt_maps::Bool
    write_cur_maps::Bool
    write_cum_cur_map_only::Bool
    write_max_cur_maps::Bool
    set_null_currents_to_nodata::Bool
    set_null_voltages_to_nodata::Bool
    compress_grids::Bool
    log_transform_maps::Bool
end

function compute_3col(resistances::Matrix{T}) where {T}
    fp = deleteat!(resistances[:,1], 1)
    l = length(fp)
    r3col = zeros(T, div(l * (l-1), 2), 3)
    k = 1
    for i = 1:l
        for j = i+1:l
            r3col[k,1] = fp[i]
            r3col[k,2] = fp[j]
            r3col[k,3] = resistances[j+1,i+1]
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
    cellmap = component_data.cellmap

    # Flags
    log_transform = flags.outputflags.log_transform_maps
    set_null_currents_to_nodata =
        flags.outputflags.set_null_currents_to_nodata
    write_max_cur_maps = flags.outputflags.write_max_cur_maps
    write_cum_cur_map_only = flags.outputflags.write_cum_cur_map_only

    node_currents, branch_currents = _create_current_maps(G, voltages, finitegrounds, cfg, nodemap = nodemap, hbmeta = hbmeta)

    if !flags.is_raster

        # Branch currents
        branch_currents_array = _convert_to_3col(branch_currents, cc)

        # Node currents
        node_currents_array = _append_name_to_node_currents(node_currents, cc)

		if flags.is_advanced
			write_currents(node_currents_array, branch_currents_array, name, cfg)
			return nothing
		end

        # TODO: implement cumulative maps for netowrk mode
		cum_branch_curr = output.cum.cum_branch_curr
		cum_node_curr = output.cum.cum_node_curr

        # Accumulate branch currents
        # cum_branch_curr[mycsid()] .+= branch_currents_array[:,3]
		bca = branch_currents_array
		cbc = cum_branch_curr[mycsid()]

		k = output.cum.coords
		@inbounds for i = 1:size(branch_currents_array, 1)
			idx = findfirst(isequal((Int(bca[i,1]), Int(bca[i,2]))), k)
			idx == nothing && (idx = findfirst(isequal((Int(bca[i,2]), Int(bca[i,1]))), k))
			cbc[idx] += bca[i,3]
		end


        # Accumulate node currents
        cnc = cum_node_curr[mycsid()] 
		nca = node_currents_array
		@inbounds for i = 1:size(nca, 1)
			cnc[Int(nca[i,1])] += nca[i,2]
		end

        # !write_cum_cur_map_only &&
        write_currents(node_currents_array, branch_currents_array, name, cfg)
		return nothing
    else

        cmap = node_currents
        cum_curr = output.cum.cum_curr
        max_curr = output.cum.max_curr

        # Process the current map
        process_grid!(cmap, cellmap, hbmeta, log_transform = log_transform,
                            set_null_to_nodata = set_null_currents_to_nodata)

        # Accumulate by default
        cum_curr[mycsid()] .+= cmap

        # Max current if user asks for it
        if write_max_cur_maps
            max_curr[mycsid()] .= max.(max_curr[mycsid()], cmap)
        end

        # Write current maps
        !write_cum_cur_map_only && flags.outputflags.write_cur_maps && 
                        write_grid(cmap, name, cfg, hbmeta)

		return nothing
   end
end

function write_currents(node_curr_arr, branch_curr_arr, name, cfg)
   pref = split(cfg["output_file"], ".out")[1]
   # 1e-6 because we guarantee only 6 digits of precision on solve
   idx = findall(x -> !isapprox(x, 0.0, atol = 1e-6), branch_curr_arr[:,3])
   branch_curr_arr = branch_curr_arr[idx, :]
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

        current_map = zeros(eltype(G), hbmeta.nrows, hbmeta.ncols)
        for j = 1:size(nodemap, 2)
            for i = 1:size(nodemap, 1)
                idx = nodemap[i,j]
                if idx == 0
                    continue
                else
                    current_map[i,j] = node_currents[idx]
                end
            end
        end

        return current_map, spzeros(0,0)
    end
end

function get_node_currents(G, voltages, finitegrounds)

    node_currents_pos = _get_node_currents_posneg(G, voltages, finitegrounds, true)
    node_currents_neg = _get_node_currents_posneg(G, voltages, finitegrounds, false)
    node_currents = map((x,y) -> x > y ? x : y, node_currents_pos, node_currents_neg)

end

function _get_node_currents_posneg(G::SparseMatrixCSC{T,V},
                            voltages, finitegrounds, pos) where {T,V}

    branch_currents = _get_branch_currents(G, voltages, pos)
    branch_currents = branch_currents - branch_currents'
    dropnonzeros!(branch_currents)

	if finitegrounds[1]!= -9999
        finiteground_currents = finitegrounds .* voltages
        if pos
            map!(x -> x < 0 ? -x : 0, finiteground_currents, finiteground_currents)
        else
            map!(x -> x > 0 ? x : 0, finiteground_currents, finiteground_currents)
        end
        n = size(G, 1)
        branch_currents = branch_currents + spdiagm(0 => finiteground_currents)
    end

    s = vec(sum(branch_currents, dims=1))

    s
end

function dropnonzeros!(G)
    for i = 1:size(G, 1)
        for j in nzrange(G, i)
            row = G.rowval[j]
            val = G.nzval[j]
            if val < 0
                G.nzval[j] = 0
            end
        end
    end
    dropzeros!(G)
end

function _get_branch_currents(G::SparseMatrixCSC{T,V}, voltages, pos) where {T,V}

    branch_currents = _get_branch_currents_posneg(G, voltages, pos)

    # Make sparse matrix with branch_currents as right upper triangle
    N = size(G, 1)
    n = size(branch_currents, 1)
    I = zeros(V, n)
    J = zeros(V, n)
    k = 1
    for i = 1:N
        for j in nzrange(G, i)
            row = G.rowval[j]
            if i > row
                I[k] = row
                J[k] = i
                k += 1
            end
        end
    end
    @assert n + 1 == k

    B = sparse(I, J, branch_currents, N, N)

    B
	# branch_currents
end

function _get_branch_currents_posneg(G::SparseMatrixCSC{T,V},
                                v::Vector{T}, pos) where {T,V}

    n = count_upper(G)
    b = zeros(T, n)
    k = 1
    if pos
        for i = 1:size(G, 1)
            for j in nzrange(G, i)
                row = G.rowval[j]
                val = G.nzval[j]
                if i > row
                    b[k] = abs(val) * (v[row] - v[i])
                    k += 1
                end
            end
        end
    else
        for i = 1:size(G, 1)
            for j in nzrange(G, i)
                row = G.rowval[j]
                val = G.nzval[j]
                if i > row
                    b[k] = abs(val) * (v[i] - v[row])
                    k += 1
                end
            end
        end
    end
    # @show n, k
    @assert n + 1 == k
    maxcur = maximum(b)
    # map!(x -> abs(x / maxcur) < 1e-8 ? 0 : x, b, b)
    for i = 1:size(b, 1)
        if abs(b[i]/ maxcur) < 1e-8
            b[i] = 0
        end
    end

    b
end

function count_upper(G)
    n = 0
    for i = 1:size(G, 1)
        for j in nzrange(G,i)
            row = G.rowval[j]
            if i > row
                n += 1
            end
        end
    end
    n
end

function process_grid!(cmap, cellmap, hbmeta; log_transform = false,
                                set_null_to_nodata = false)
    if log_transform
        map!(x -> x > 0 ? log10(x) : float(hbmeta.nodata), cmap, cmap)
    end

    if set_null_to_nodata
        for i in eachindex(cmap)
            if cellmap[i] == 0
                cmap[i] = hbmeta.nodata
            end
        end
    end

end

function write_grid(cmap, name, cfg, hbmeta;
                        voltage = false, cum = false,
                        max = false)


    str = "curmap"
    if cum
        str = "cum_$(str)"
    elseif max
        str = "max_$(str)"
    elseif voltage
        str = "voltmap"
    end

    pref = split(cfg["output_file"], ".out")[1]
    filename = "$(pref)_$(str)$(name)"

    cfg["write_as_tif"] in TRUELIST ? (file_format = "tif") :
            (file_format = "asc")

    write_raster(filename,
                 cmap,
                 hbmeta.wkt,
                 hbmeta.transform,
                 file_format)
end


function write_grid(cmap, name, cfg, hbmeta, cellmap;
                        voltage = false, cum = false, max = false,
                        log_transform = false, set_null_to_nodata = false)

    pref = split(cfg["output_file"], ".out")[1]

    if log_transform
        map!(x -> x > 0 ? log10(x) : float(hbmeta.nodata), cmap, cmap)
    end

    if set_null_to_nodata
        for i in eachindex(cmap)
            if cellmap[i] == 0
                cmap[i] = hbmeta.nodata
            end
        end
    end

    str = "curmap"
    if cum
        str = "cum_$(str)"
    elseif max
        str = "max_$(str)"
    elseif voltage
        str = "voltmap"
    end

    filename = "$(pref)_$(str)$(name)"

    cfg["write_as_tif"] in TRUELIST ? (file_format = "tif") :
            (file_format = "asc")

    write_raster(filename,
                 cmap,
                 hbmeta.wkt,
                 hbmeta.transform,
                 file_format)
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
        set_null_voltages_to_nodata = flags.outputflags.set_null_voltages_to_nodata

        vm = _create_voltage_map(voltages, nodemap, hbmeta)
        write_grid(vm, name, cfg, hbmeta, component_data.cellmap, voltage = true,
                        set_null_to_nodata = set_null_voltages_to_nodata)
    end
end

function write_voltages(output, name, voltages::Vector{T}, cc) where {T}

    volt_arr = zeros(T, size(voltages, 1), 2)
    volt_arr[:,1] = cc
    volt_arr[:,2] = voltages

    pref = split(output, ".out")[1]
    writedlm("$(pref)_voltages$(name).txt", volt_arr)

end

function _create_voltage_map(voltages::Vector{T}, nodemap, hbmeta) where {T}
    voltmap = zeros(T, hbmeta.nrows, hbmeta.ncols)
    for j = 1:size(nodemap, 2)
        for i = 1:size(nodemap, 1)
            idx = nodemap[i,j]
            if idx == 0
                continue
            else
                voltmap[i,j] = voltages[idx]
            end
        end
    end
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

function save_resistances(r, cfg)
    pref = split(cfg["output_file"], ".out")[1]
    filename = "$(pref)_resistances.out"
    filename_3col = "$(pref)_resistances_3columns.out"
    rcol = compute_3col(r)
    open(filename, "w") do f
        writedlm(f, r, ' ')
    end
    open(filename_3col, "w") do f
        writedlm(f, rcol, ' ')
    end
end

function write_cum_maps(cum, cellmap::Matrix{T}, cfg, hbmeta, write_max, write_cum) where T

    if write_cum || cfg["write_cur_maps"] in TRUELIST
        cum_curr = zeros(T, size(cellmap)...)
        for i = 1:nprocs()
            cum_curr .+= cum.cum_curr[i]
        end
        postprocess_cum_curmap!(cum_curr)
        write_grid(cum_curr, "", cfg, hbmeta, cum = true)
    end

    if write_max
        max_curr = fill(T(-9999), size(cellmap)...)
        for i = 1:nprocs()
            max_curr .= max.(cum.max_curr[i], max_curr)
        end
        postprocess_cum_curmap!(max_curr)
        write_grid(max_curr, "", cfg, hbmeta, max = true)
    end

end

# Write a single band raster, either in .tif or .asc format,
# inspired by GeoArrays.write()
function write_raster(fn_prefix::String,
                      array::Matrix{T} where T <: Number,
                      wkt::String,
                      transform,
                      file_format::String)
    # transponse array back to columns by rows
    array_t = permutedims(array, [2, 1])

    width, height = size(array_t)

    # Define extension and driver based in file_format
    file_format == "tif" ? (ext = ".tif"; driver = "GTiff") :
            (ext = ".asc"; driver = "AAIGrid")

    file_format == "tif" ? (options = ["COMPRESS=LZW"]) :
                           (options = [])

    # Append file extention to filename
    fn = string(fn_prefix, ext)

    # Create raster in memory *NEEDED* because no create driver for .asc
    ArchGDAL.create(fn_prefix,
                    driver = ArchGDAL.getdriver("MEM"),
                    width = width,
                    height = height,
                    nbands = 1,
                    dtype = eltype(array_t),
                    options = options) do dataset
        band = ArchGDAL.getband(dataset, 1)
        # Write data to band
        ArchGDAL.write!(band, array_t)

        # Write nodata and projection info
        ArchGDAL.setnodatavalue!(band, -9999.0)
        ArchGDAL.setgeotransform!(dataset, transform)
        ArchGDAL.setproj!(dataset, wkt)

        # Copy memory object to disk (necessary because ArchGDAL.create
        # does not support creation of ASCII rasters)
        ArchGDAL.write(dataset, fn,
                       driver = ArchGDAL.getdriver(driver),
                       options = options)
    end

end
