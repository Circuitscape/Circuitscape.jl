"""
Primary driver for network pairwise. 
"""
function network_pairwise(T, V, cfg)::Matrix{T}
    
    # Get input
    networkdata = get_network_data(T, V, cfg)

    # Get compute flags
    flags = get_network_flags(cfg)

    # Compute graph data based on compute flags
    graphdata = compute_graph_data(networkdata)

    # Send to main kernel
    single_ground_all_pairs(graphdata, flags, cfg)
end

function compute_graph_data(data::NetworkData{T,V})::GraphData{T,V} where {T,V}

    i,j,v = data.coords

    idx = findfirst(x -> x < 1, i)
    idx != nothing && throw("Indices no good")
    idx = findfirst(x -> x < 1, j)
    idx != nothing && throw("Indices no good")
    
    m = max(i[end], j[end])
    A = sparse(i,j,v,m,m)
    A = A + A'

    cc = connected_components(SimpleWeightedGraph(A))

    t = @elapsed G = laplacian(A)
    csinfo("Time taken to construct graph laplacian = $t")

    # T = eltype(i)
    exclude_pairs = Tuple{V,V}[]

    nodemap = Matrix{V}(undef,0,0)
    polymap = Matrix{V}(undef,0,0)
    hbmeta = RasterMeta()
    cellmap = Matrix{T}(undef,0,0)

    cum = initialize_cum_vectors(v)

    GraphData(G, cc, data.fp, data.fp, 
                exclude_pairs, nodemap, polymap, hbmeta, cellmap, cum)
end

function get_network_flags(cfg)
    
    # Computation flags
    is_raster = false
    is_advanced = cfg["scenario"] in ADVANCED
    is_alltoone = false
    is_onetoall = false
    grnd_file_is_res = cfg["ground_file_is_resistances"] in TRUELIST
    policy = Symbol(cfg["remove_src_or_gnd"])
    solver = cfg["solver"]

    # Output flags
    write_volt_maps = cfg["write_volt_maps"] in TRUELIST
    write_cur_maps = cfg["write_cur_maps"] in TRUELIST
    write_cum_cur_maps_only = cfg["write_cum_cur_map_only"] in TRUELIST
    write_max_cur_maps = cfg["write_max_cur_maps"] in TRUELIST
    set_null_currents_to_nodata = cfg["set_null_currents_to_nodata"] in TRUELIST
    set_null_voltages_to_nodata = cfg["set_null_voltages_to_nodata"] in TRUELIST
    compress_grids = cfg["compress_grids"] in TRUELIST
    log_transform_maps = cfg["log_transform_maps"] in TRUELIST

    o = OutputFlags(write_volt_maps, write_cur_maps,
                    write_cum_cur_maps_only, write_max_cur_maps,
                    set_null_currents_to_nodata, set_null_voltages_to_nodata,
                    compress_grids, log_transform_maps)
    
    NetworkFlags(is_raster, is_advanced, is_alltoone, is_onetoall, 
                grnd_file_is_res, policy, solver, o)
end

struct NetworkFlags
    is_raster::Bool
    is_advanced::Bool
    is_alltoone::Bool
    is_onetoall::Bool
    grnd_file_is_res::Bool
    policy::Symbol
    solver::String
    outputflags::OutputFlags
end

