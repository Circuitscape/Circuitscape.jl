export  accumulate_current_maps,
        calculate_cum_current_map,
        calculate_max_current_map

mutable struct MutablePair{T,V}
	first::T
	second::V
end

 """
 Construct nodemap specific to a connected component
 """
function construct_local_node_map(nodemap, component, polymap)
    local_nodemap = zeros(eltype(nodemap), size(nodemap))
    idx = findall(in(component), nodemap)
    local_nodemap[idx] = nodemap[idx]
    if nodemap == local_nodemap
        return local_nodemap
    end
    _construct_local_nodemap(local_nodemap, polymap, idx)
end

function _construct_local_nodemap(local_nodemap, polymap, idx)
    if isempty(polymap)
        i = findall(x ->x!=0, local_nodemap)
        local_nodemap[i] = 1:length(i)
        return local_nodemap
    else
        local_polymap = zeros(eltype(local_nodemap), size(local_nodemap))
        local_polymap[idx] = polymap[idx]
        return construct_node_map(local_nodemap, local_polymap)
    end
end


function get_output_flags(cfg)

    # Output flags
    write_volt_maps = cfg.write_volt_maps
    write_cur_maps = cfg.write_cur_maps
    write_cum_cur_map_only = cfg.write_cum_cur_map_only
    write_max_cur_maps = cfg.write_max_cur_maps
    set_null_currents_to_nodata = cfg.set_null_currents_to_nodata
    set_null_voltages_to_nodata = cfg.set_null_voltages_to_nodata
    compress_grids = cfg.compress_grids
    log_transform_maps = cfg.log_transform_maps

    o = OutputFlags(write_volt_maps, write_cur_maps,
                    write_cum_cur_map_only, write_max_cur_maps,
                    set_null_currents_to_nodata, set_null_voltages_to_nodata,
                    compress_grids, log_transform_maps)
end

# Helps start new processes from the INI file
function myaddprocs(n)
    addprocs(n)
    @everywhere Core.eval(Main, :(using Circuitscape))
end

# Reads the directory with the current maps
# and accumulates all current maps
function accumulate_current_maps(path, f)
    dir = dirname(path)
    base = basename(path)

    # If base file has a dot
    name = split(base, ".out")[1]

    cmap_list = readdir(dir) |>
                    x -> filter(y -> startswith(y, "$(name)_"), x) |>
                    x -> filter(y -> occursin("_curmap_", y), x)
    isempty(cmap_list) && return

    headers = ""
    first_file = joinpath(dir, cmap_list[1])
    nrow = 0
    ncol = 0

    # Read the headers from the first file
    open(first_file, "r") do f

        # Get num cols
        str = readline(f)
        headers = headers * str * "\n"
        ncol = split(str)[2] |> x -> parse(Int, x)

        # Get num rows
        str = readline(f)
        headers = headers * str * "\n"
        nrow = split(str)[2] |> x -> parse(Int, x)

        # Just append the rest
        for i = 3:6
            headers = headers * readline(f) * "\n"
        end
    end

    accum = zeros(nrow, ncol)
    for file in cmap_list
        csinfo("Accumulating $file", cfg.suppress_messages)
        cmap_path = joinpath(dir, file)
        cmap = readdlm(cmap_path, skipstart = 6)
        f_in_place!(accum, cmap, f)
    end
    for i in eachindex(accum)
        if accum[i] < -9999
            accum[i] = -9999
        end
    end

    name =  if isequal(f, +)
                "cum"
            elseif isequal(f, max)
                "max"
            end

    accum_path = joinpath(dir, name * "_$(name)_curmap.asc")
    csinfo("Writing to $accum_path", cfg.suppress_messages)
    open(accum_path, "w") do f
        write(f, headers)
        writedlm(f, round.(accum, digits=8), ' ')
    end

end

function f_in_place!(accum, cmap, f)
    accum .= f.(accum, cmap)
end

calculate_cum_current_map(path) = accumulate_current_maps(path, +)
calculate_max_current_map(path) = accumulate_current_maps(path, max)

function postprocess_cum_curmap!(accum)
    for i in eachindex(accum)
        if accum[i] < -9999
            accum[i] = -9999
        end
    end
end

mycsid() = myid() - minimum(workers()) + 1

function initialize_cum_maps(cellmap::Matrix{T}, max = false) where T
    cum_curr = Vector{SharedMatrix{T}}(undef,nprocs())
    for i = 1:nprocs()
        cum_curr[i] = SharedArray(zeros(T, size(cellmap)...))
    end
    max_curr = Vector{SharedMatrix{T}}()
    if max
        max_curr = Vector{SharedMatrix{T}}(undef,nprocs())
        for i = 1:nprocs()
            max_curr[i] = SharedArray(fill(T(-9999), size(cellmap)...))
        end
    end
    cum_branch_curr = Vector{SharedVector{T}}()
    cum_node_curr = Vector{SharedVector{T}}()

    Cumulative(cum_curr, max_curr,
			   cum_branch_curr, cum_node_curr, Vector{Tuple{Int,Int}}())
end

function initialize_cum_vectors(coords::Tuple{Vector{V},Vector{V},Vector{T}}, num_nodes::Int64) where {T,V}
    cum_curr = Vector{SharedMatrix{T}}()
    max_curr = Vector{SharedMatrix{T}}()
	cum_branch_curr = Vector{SharedVector{T}}(undef,nprocs())
	cum_node_curr = Vector{SharedVector{T}}(undef,nprocs())
	_i, _j, _v = coords
    for i = 1:nprocs()
        #cum_branch_curr[i] = SharedArray(zeros(T, size(v)...))
		cum_branch_curr[i] = SharedVector(zeros(T,length(_v)))
		cum_node_curr[i] = SharedVector(zeros(T,num_nodes))
    end

	coords = map((x,y) -> (x,y), _i, _j)
    Cumulative(cum_curr, max_curr,
			   cum_branch_curr, cum_node_curr, coords)
end

# Function to calculate current for Omniscape moving window solves
function compute_omniscape_current(
        conductance::Array{T, 2} where T <: Union{Float32, Float64},
        source::Array{T, 2} where T <: Union{Float32, Float64},
        ground::Array{T, 2} where T <: Union{Float32, Float64},
        cs_cfg::Dict{String, String}
    )
    V = Int64
    T = eltype(conductance)

    # get raster data
    cellmap = conductance
    polymap = Matrix{V}(undef, 0, 0)
    source_map = source
    ground_map = ground
    points_rc = (V[], V[], V[])
    strengths = Matrix{T}(undef, 0, 0)

    included_pairs = IncludeExcludePairs(:undef,
                                         V[],
                                         Matrix{V}(undef,0,0))

    # This is just to satisfy type requirements, most of it not used
    hbmeta = RasterMeta(size(cellmap)[2],
                        size(cellmap)[1],
                        0.,
                        0.,
                        1.,
                        -9999.,
                        Array{T, 1}(undef, 1),
                        "")

    rasterdata = RasterData(cellmap,
                            polymap,
                            source_map,
                            ground_map,
                            points_rc,
                            strengths,
                            included_pairs,
                            hbmeta)

    # Generate advanced data
    o = OutputFlags(
        false, false, false, false,
        false, false, false, false
    )

    cfg = CSConfig(cs_cfg)

    flags = RasterFlags(
        true, false, true, false, false, false, Symbol("rmvsrc"),
        cfg.connect_four_neighbors_only, false, 
        _solver_str(cfg.solver), o
    )

    data = compute_advanced_data(rasterdata, flags, cfg)

    G = data.G
    nodemap = data.nodemap
    polymap = data.polymap
    hbmeta = data.hbmeta
    sources = data.sources
    grounds = data.grounds
    finitegrounds = data.finitegrounds
    cc = data.cc
    check_node = data.check_node
    source_map = data.source_map # Need it for one to all mode
    cellmap = data.cellmap

    f_local = Vector{eltype(G)}()
    voltages = Vector{eltype(G)}()
    outcurr = alloc_map(hbmeta)

    for c in cc
        if check_node != -1 && !(check_node in c)
            continue
        end

        # a_local = laplacian(G[c, c])
        a_local = G[c,c]
        s_local = sources[c]
        g_local = grounds[c]

        if sum(s_local) == 0 || sum(g_local) == 0
            continue
        end

        if finitegrounds != [-9999.]
            f_local = finitegrounds[c]
        else
            f_local = finitegrounds
        end

        voltages = multiple_solver(cfg,
                                   data.solver,
                                   a_local,
                                   s_local,
                                   g_local,
                                   f_local)

        local_nodemap = construct_local_node_map(nodemap,
                                                 c,
                                                 polymap)

        accum_currents!(outcurr,
                        voltages,
                        cfg,
                        a_local,
                        voltages,
                        f_local,
                        local_nodemap,
                        hbmeta)
    end

    return outcurr
end
