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


"""
    write_cur_maps(name, output, component_data, finitegrounds, cfg)

Current map of one solve on a connected component, accumulated into the
cumulative maps when `output` is a pairwise `Output` and written per solve
when the options ask for it. Where the currents go is decided by the
component's geometry: a grid for rasters, node and branch lists for networks.
"""
write_cur_maps(name, output, component_data, finitegrounds, cfg) =
    write_cur_maps(component_data.geometry, name, output, component_data, finitegrounds, cfg)

function write_cur_maps(geometry::NetworkGeometry, name, output, component_data, finitegrounds, cfg)

    node_currents, branch_currents = _create_current_maps(component_data.matrix,
                                            _voltages(output), finitegrounds, geometry)

    nodes = geometry.nodes
    branch_currents_array = _convert_to_3col(branch_currents, nodes)
    node_currents_array = _append_name_to_node_currents(node_currents, nodes)

    accumulate_currents!(output, node_currents_array, branch_currents_array)

    write_currents(node_currents_array, branch_currents_array, name, cfg)
    nothing
end

# Advanced mode writes each solve's currents as they are; only the pairwise
# kernel keeps cumulative node and branch currents.
accumulate_currents!(::AbstractVector, node_currents_array, branch_currents_array) = nothing

function accumulate_currents!(output::Output, node_currents_array, branch_currents_array)
    cum_branch_curr = output.cum.cum_branch_curr
    cum_node_curr = output.cum.cum_node_curr

    # Resolve every branch to its slot before taking the lock, so the
    # critical section is just the adds.
    bca = branch_currents_array
    bidx = output.cum.branch_index
    V = keytype(bidx).parameters[1]
    slots = [bidx[(V(bca[i,1]), V(bca[i,2]))] for i in 1:size(bca, 1)]

    lock(output.cum.lock) do
        # Accumulate branch currents
        cbc = cum_branch_curr
        @inbounds for i = 1:size(bca, 1)
            cbc[slots[i]] += bca[i,3]
        end

        # Accumulate node currents
        cnc = cum_node_curr
        nca = node_currents_array
        @inbounds for i = 1:size(nca, 1)
            cnc[Int(nca[i,1])] += nca[i,2]
        end
    end
    nothing
end

function write_cur_maps(geometry::RasterGeometry, name, output::Output, component_data, finitegrounds, cfg)

    nodemap = geometry.nodemap
    hbmeta = geometry.hbmeta

    # Output options
    log_transform = cfg.log_transform_maps
    set_null_currents_to_nodata = cfg.set_null_currents_to_nodata
    write_max_cur_maps = cfg.write_max_cur_maps
    write_cum_cur_map_only = cfg.write_cum_cur_map_only

    cmap, _ = _create_current_maps(component_data.matrix, output.voltages, finitegrounds, geometry)
    cum_curr = output.cum.cum_curr
    max_curr = output.cum.max_curr

    # Issue 342: with set_focal_node_currents_to_zero the two focal nodes
    # being solved carry no current in this pair's map, so the cumulative
    # map shows only current that flows *through* a focal region when
    # other pairs are active (Dickson et al. 2013).
    if cfg.set_focal_node_currents_to_zero
        zero_focal_cells!(cmap, nodemap, output.comp_idx)
    end

    # Process the current map
    process_grid!(cmap, geometry.cellmap, hbmeta, log_transform = log_transform,
                        set_null_to_nodata = set_null_currents_to_nodata)

    # Accumulate by default
    lock(output.cum.lock) do
        cum_curr .+= cmap

        # Max current if user asks for it
        if write_max_cur_maps
            max_curr .= max.(max_curr, cmap)
        end
    end

    # Write current maps
    !write_cum_cur_map_only && cfg.write_cur_maps &&
                    write_grid(cmap, name, cfg, hbmeta)

    nothing
end

# Zero every cell of `cmap` whose local node index is one of `nodes`. Focal
# regions collapse to a single node, so this clears the whole region.
function zero_focal_cells!(cmap, nodemap, nodes)
    @inbounds for k in eachindex(nodemap)
        nodemap[k] in nodes && (cmap[k] = 0)
    end
    cmap
end

function write_currents(node_curr_arr, branch_curr_arr, name, cfg)
    pref = split(cfg.output_file, ".out")[1]
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

# Node and branch currents of one solve. For a network both are returned as
# vectors over the matrix rows / sparse branch matrix; for a raster the node
# currents are spread onto the grid and there are no branch currents.
function _create_current_maps(G, voltages, finitegrounds, ::NetworkGeometry)
    node_currents = get_node_currents(G, voltages, finitegrounds)
    branch_currents = _get_branch_currents(G, voltages, true)
    branch_currents = abs.(branch_currents)
    node_currents, branch_currents
end

function _create_current_maps(G, voltages, finitegrounds, geometry::RasterGeometry)
    node_currents = get_node_currents(G, voltages, finitegrounds)
    nodemap = geometry.nodemap
    hbmeta = geometry.hbmeta

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

    current_map, nothing
end

"""
    get_node_currents(G, voltages, finitegrounds)

Per-node current `max(total inflow, total outflow)`, in one sweep over the
CSC structure of `G` with no sparse temporaries.

For every off-diagonal entry `(row, col)` of the symmetric Laplacian `G` the
branch current `b = |G[row,col]| * (v[row] - v[col])` is the current arriving
at `col` from `row`; `b > 0` is inflow at `col`, `b < 0` is outflow. Every
undirected edge is therefore seen from both ends, once in each column, which
is exactly the column sum of the old `B - B'` construction and keeps the
same (row-sorted) summation order. A finite ground at node `i` carries
`finitegrounds[i] * v[i]`: negative is inflow, positive is outflow.

The historical implementation zeroed branch currents smaller than `1e-8`
relative to the largest *signed* branch current of each direction before
summing; that reference differs between the two directions, so it is
recovered first (a pass that touches only `nzval`/`voltages`) and applied
in the accumulation pass to reproduce the previous results exactly.
`finitegrounds[1] != -9999` is the sentinel meaning "no finite grounds".
"""
function get_node_currents(G::SparseMatrixCSC{T}, voltages::AbstractVector,
                           finitegrounds) where {T}
    n = size(G, 1)
    v = voltages
    rowval = rowvals(G)
    nzval = nonzeros(G)

    # Pass 1: maxima of the signed branch currents over the upper triangle
    # (`maxin`: rows above the diagonal, the old "pos" reference) and the
    # lower triangle (`maxout`, the old "neg" reference). They are -Inf when
    # the graph has no edges, in which case nothing is thresholded.
    maxin = typemin(T)
    maxout = typemin(T)
    @inbounds for col = 1:n
        vc = v[col]
        for k in nzrange(G, col)
            row = rowval[k]
            row == col && continue
            b = abs(nzval[k]) * (v[row] - vc)
            if row < col
                b > maxin && (maxin = b)
            else
                b > maxout && (maxout = b)
            end
        end
    end

    # Pass 2: accumulate inflow and outflow per node and keep the larger.
    node_currents = Vector{T}(undef, n)
    has_finitegrounds = finitegrounds[1] != -9999
    @inbounds for col = 1:n
        vc = v[col]
        inflow = zero(T)
        outflow = zero(T)
        for k in nzrange(G, col)
            row = rowval[k]
            row == col && continue
            b = abs(nzval[k]) * (v[row] - vc)
            if b > 0
                abs(b / maxin) < 1e-8 || (inflow += b)
            elseif b < 0
                abs(b / maxout) < 1e-8 || (outflow -= b)
            end
        end
        if has_finitegrounds
            fg = finitegrounds[col] * vc
            if fg < 0
                inflow -= fg
            elseif fg > 0
                outflow += fg
            end
        end
        node_currents[col] = inflow > outflow ? inflow : outflow
    end

    node_currents
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

function write_grid(cmap, name, cfg, hbmeta, cellmap = nothing;
                        voltage = false, cum = false, max = false,
                        log_transform = false, set_null_to_nodata = false)

    if cellmap !== nothing
        process_grid!(cmap, cellmap, hbmeta, log_transform = log_transform,
                            set_null_to_nodata = set_null_to_nodata)
    end

    str = "curmap"
    if cum
        str = "cum_$(str)"
    elseif max
        str = "max_$(str)"
    elseif voltage
        str = "voltmap"
    end

    pref = split(cfg.output_file, ".out")[1]
    filename = "$(pref)_$(str)$(name)"

    cfg.write_as_tif ? (file_format = "tif") :
            (file_format = "asc")

    write_raster(filename,
                 cmap,
                 hbmeta.wkt,
                 hbmeta.transform,
                 file_format)
end

write_volt_maps(name, output, component_data, cfg) =
    write_volt_maps(component_data.geometry, name, _voltages(output), cfg)

write_volt_maps(geometry::NetworkGeometry, name, voltages, cfg) =
    write_voltages(cfg.output_file, name, voltages, geometry.nodes)

function write_volt_maps(geometry::RasterGeometry, name, voltages, cfg)
    hbmeta = geometry.hbmeta
    vm = _create_voltage_map(voltages, geometry.nodemap, hbmeta)
    write_grid(vm, name, cfg, hbmeta, geometry.cellmap, voltage = true,
                    set_null_to_nodata = cfg.set_null_voltages_to_nodata)
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

# Advanced mode solves every component into one map: each component's
# voltages / currents are spread onto the grid through its local geometry and
# added to `base`. Networks keep whole-graph vectors instead, so nothing to do.
function accum_voltages!(base, voltages, geometry::RasterGeometry)
    voltmap = _create_voltage_map(voltages, geometry.nodemap, geometry.hbmeta)
    for i in eachindex(base)
        base[i] += voltmap[i]
    end
end
accum_voltages!(base, voltages, ::NetworkGeometry) = nothing

function accum_currents!(base, G, voltages, finitegrounds, geometry::RasterGeometry)
    node_currents, _ = _create_current_maps(G, voltages, finitegrounds, geometry)
    for i in eachindex(base)
        base[i] += node_currents[i]
    end
end
accum_currents!(base, G, voltages, finitegrounds, ::NetworkGeometry) = nothing

function save_resistances(r, cfg)
    pref = split(cfg.output_file, ".out")[1]
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

"""
    write_cum_maps(cum, geometry, cfg)

Write the cumulative (and maximum) current maps accumulated over a pairwise
run, in whatever form the geometry calls for.
"""
write_cum_maps(cum, geometry::RasterGeometry, cfg) = write_cum_maps(cum, geometry.hbmeta, cfg)

function write_cum_maps(cum, geometry::NetworkGeometry, cfg)
    cfg.write_cur_maps || return nothing
    cum_node_curr = cum.cum_node_curr
    cum_branch_curr = cum.cum_branch_curr
    cum_node_curr = hcat(1:length(cum_node_curr), cum_node_curr)
    coords = cum.coords
    cum_branch_curr = hcat(getindex.(coords, 1), getindex.(coords, 2), cum_branch_curr)
    write_currents(cum_node_curr, cum_branch_curr, "_cum", cfg)
end

# Grid form, shared with the paths that keep a cumulative map of their own
# (focal regions, one-to-all).
function write_cum_maps(cum, hbmeta::RasterMeta, cfg)
    (cfg.write_cur_maps || cfg.write_cum_cur_map_only) || return nothing

    cum_curr = cum.cum_curr
    postprocess_cum_curmap!(cum_curr)
    write_grid(cum_curr, "", cfg, hbmeta, cum = true)

    if cfg.write_max_cur_maps
        max_curr = cum.max_curr
        postprocess_cum_curmap!(max_curr)
        write_grid(max_curr, "", cfg, hbmeta, max = true)
    end

end

"""
    write_raster(fn_prefix, array, wkt, transform, file_format)

Write a single-band raster. `.asc` is written natively with `writedlm` —
several times faster than the GDAL AAIGrid driver, and thread-safe, so it
takes no lock; `.tif` goes through GDAL under `IO_LOCK`.
"""
function write_raster(fn_prefix::String,
                      array::Matrix{T} where T <: Number,
                      wkt::String,
                      transform,
                      file_format::String)
    if file_format == "tif"
        write_raster_gdal(fn_prefix, array, wkt, transform)
    else
        write_asc(fn_prefix, array, wkt, transform)
    end
end

function write_asc(fn_prefix::String, array::AbstractMatrix, wkt::String, transform)
    nrows, ncols = size(array)
    if length(transform) >= 6
        xll, dx, dy = transform[1], transform[2], -transform[6]
        yll = transform[4] - nrows * dy
    else
        xll, yll, dx, dy = 0.0, 0.0, 1.0, 1.0
    end
    open(fn_prefix * ".asc", "w") do io
        println(io, "ncols         ", ncols)
        println(io, "nrows         ", nrows)
        println(io, "xllcorner     ", xll)
        println(io, "yllcorner     ", yll)
        if dx ≈ dy
            println(io, "cellsize      ", dx)
        else
            println(io, "dx            ", dx)
            println(io, "dy            ", dy)
        end
        println(io, "NODATA_value  -9999")
        writedlm(io, array, ' ')
    end
    # Projection travels in a sidecar, as with the GDAL AAIGrid driver
    isempty(wkt) || write(fn_prefix * ".prj", wkt)
    nothing
end

# Replaced by a more specific method in CircuitscapeArchGDALExt when ArchGDAL
# is loaded.
write_raster_gdal(fn_prefix, array, wkt, transform) =
    error(gdal_missing_message("Writing GeoTIFF (write_as_tif = True)"))
