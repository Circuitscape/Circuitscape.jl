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

"""
    compute_omniscape_current(conductance, source, ground, cs_cfg) -> Matrix{Float64}

Advanced-mode solve on in-memory arrays; Omniscape's entry point for every
moving-window solve. No files are read or written.

Arguments:
- `conductance::Matrix{T}`, `T <: Union{Float32, Float64}`: per-cell
  conductance (not resistance). Cells with a value of `0` are dropped from the
  graph; every other cell becomes a node connected to its eight neighbours
  (four if `connect_four_neighbors_only = True` in `cs_cfg`) with the average
  conductance of the two cells as edge weight.
- `source::Matrix{T}`: current injected at each cell, in amps; `0` where there
  is no source.
- `ground::Matrix{T}`: conductance to ground at each cell (`Inf` for a direct
  ground, `0` where there is no ground). Where a cell is both a source and a
  ground the source is removed (`remove_src_or_gnd = rmvsrc`).
- `cs_cfg::Dict{String,String}`: INI keys and values as from
  [`init_config`](@ref) merged with [`update!`](@ref). It is parsed with the
  same strict rules as an INI file, so `scenario` must be set and values must
  be valid. Of it, only `solver`, `cholmod_batch_size` and
  `residual_tolerance` (solver settings) and `connect_four_neighbors_only`
  affect the result; `data_type`, `scenario`, `ground_file_is_resistances`,
  `remove_src_or_gnd`, `connect_using_avg_resistances` and every map-writing
  option are overridden. The arithmetic precision is `T`, the element type of
  `conductance`, not `cs_cfg["precision"]`; node indices are always `Int64`.

Returns a `Matrix{Float64}` of `size(conductance)`: the current through each
cell (`max(inflow, outflow)` at its node), summed over every connected
component that contains both a source and a ground. Cells with no node, or
in components without a source or a ground, are `0`. The solve of each
component is accepted only if its residual passes the gate described at
[`residual_tolerance`](@ref).
"""
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

    # Omniscape drives a raster advanced-mode solve on in-memory arrays with
    # the sources removed at conflicts and no maps written; only the solver
    # settings and connect_four_neighbors_only are taken from the caller.
    cfg = CSConfig(CSConfig(cs_cfg);
                   data_type = dt_raster,
                   scenario = sc_advanced,
                   ground_file_is_resistances = false,
                   remove_src_or_gnd = rp_rmvsrc,
                   connect_using_avg_resistances = false,
                   write_volt_maps = false,
                   write_cur_maps = false,
                   write_cum_cur_map_only = false,
                   write_max_cur_maps = false,
                   set_null_currents_to_nodata = false,
                   set_null_voltages_to_nodata = false,
                   set_focal_node_currents_to_zero = false,
                   compress_grids = false,
                   log_transform_maps = false)

    # The kernel accumulates the current map for us even though no map is
    # written; that map is Omniscape's result. Omniscape calls this from many
    # tasks at once, and a TimerOutput is not safe to share between them, so
    # each call times into a private timer instead of the default CSTIMER.
    with(CSTIMER => TimerOutput()) do
        prob = advanced_problem(rasterdata, cfg)
        _, outcurr, _ = advanced_kernel(prob, cfg; accumulate_currents = true)
        outcurr
    end
end
