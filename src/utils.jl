export  accumulate_current_maps,
        calculate_cum_current_map,
        calculate_max_current_map

const IO_LOCK = ReentrantLock()

"""
    connected_components(A::SparseMatrixCSC) -> Vector{Vector{Ti}}

Connected components of the graph whose edges are the nonzero off-diagonal
entries of the symmetric matrix `A`, computed directly on the CSC structure.
Components are ordered by their smallest vertex and vertices ascend within
each — the same contract as `Graphs.connected_components(SimpleGraph(A))`,
which this replaces. That path built an adjacency-list copy of `A` first: at
4M cells, 1.4 GB allocated and 458 MB retained to produce 31 MB of labels.
"""
function connected_components(A::SparseMatrixCSC{Tv,Ti}) where {Tv,Ti}
    n = size(A, 1)
    size(A, 2) == n || throw(ArgumentError("adjacency matrix must be square"))
    rv = rowvals(A)
    nz = nonzeros(A)
    label = zeros(Ti, n)
    stack = Ti[]
    ncomp = 0
    for u in 1:n
        label[u] == 0 || continue
        ncomp += 1
        label[u] = ncomp
        push!(stack, Ti(u))
        while !isempty(stack)
            v = pop!(stack)
            @inbounds for k in nzrange(A, v)
                w = rv[k]
                (nz[k] == 0 || label[w] != 0) && continue
                label[w] = ncomp
                push!(stack, w)
            end
        end
    end
    comps = [Ti[] for _ in 1:ncomp]
    for u in 1:n
        push!(comps[label[u]], Ti(u))
    end
    comps
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
    OutputFlags(cfg.write_volt_maps, cfg.write_cur_maps,
                cfg.write_cum_cur_map_only, cfg.write_max_cur_maps,
                cfg.set_null_currents_to_nodata, cfg.set_null_voltages_to_nodata,
                cfg.compress_grids, cfg.log_transform_maps,
                cfg.set_focal_node_currents_to_zero)
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
        @info("Accumulating $file")
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
    @info("Writing to $accum_path")
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

function initialize_cum_maps(cellmap::Matrix{T}, max = false) where T
    cum_curr = zeros(T, size(cellmap)...)
    max_curr = max ? fill(T(-9999), size(cellmap)...) : zeros(T, 0, 0)
    cum_branch_curr = Vector{T}()
    cum_node_curr = Vector{T}()

    Cumulative(cum_curr, max_curr, cum_branch_curr, cum_node_curr,
               Vector{Tuple{Int,Int}}(), Dict{Tuple{Int,Int},Int}(), ReentrantLock())
end

function initialize_cum_vectors(coords::Tuple{Vector{V},Vector{V},Vector{T}}, num_nodes::Int64) where {T,V}
    cum_curr = zeros(T, 0, 0)
    max_curr = zeros(T, 0, 0)
	_i, _j, _v = coords
	cum_branch_curr = zeros(T, length(_v))
	cum_node_curr = zeros(T, num_nodes)

    pair_coords = map((x,y) -> (x,y), _i, _j)
    # First occurrence wins for a duplicated edge, as the linear scan it
    # replaces did; both orientations map to the same slot.
    branch_index = Dict{Tuple{V,V},Int}()
    sizehint!(branch_index, 2 * length(pair_coords))
    for (idx, (a, b)) in enumerate(pair_coords)
        get!(branch_index, (a, b), idx)
        get!(branch_index, (b, a), idx)
    end
    Cumulative(cum_curr, max_curr, cum_branch_curr, cum_node_curr,
               pair_coords, branch_index, ReentrantLock())
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
        false, false, false, false, false
    )

    cfg = CSConfig(cs_cfg)

    flags = RasterFlags(
        true, false, true, false, false, false, Symbol("rmvsrc"),
        cfg.connect_four_neighbors_only, false, o
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
