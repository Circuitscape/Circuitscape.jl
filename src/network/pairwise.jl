"""
Primary driver for network pairwise. 
"""
function network_pairwise(T, cfg)
    
    # Get input
    networkdata = get_network_data(T, cfg)

    # Get compute flags
    flags = get_network_flags(cfg)

    # Compute graph data based on compute flags
    graphdata = compute_graph_data(networkdata)

    # Send to main kernel
    single_ground_all_pairs(graphdata, flags, cfg)
end

function compute_graph_data(data::NetworkData)

    i,j,v = data.coords

    idx = findfirst(x -> x < 1, i)
    idx != 0 && throw("Indices no good")
    idx = findfirst(x -> x < 1, j)
    idx != 0 && throw("Indices no good")
    
    m = max(i[end], j[end])
    A = sparse(i,j,v,m,m)
    A = A + A'

    cc = connected_components(SimpleWeightedGraph(A))

    t = @elapsed G = laplacian(A)
    csinfo("Time taken to construct graph laplacian = $t")

    T = eltype(i)
    exclude_pairs = Tuple{T,T}[]

    nodemap = Matrix{T}(0,0)
    polymap = Matrix{T}(0,0)
    hbmeta = RasterMeta()

    GraphData(G, cc, data.fp, data.fp, 
                exclude_pairs, nodemap, polymap, hbmeta)
end

function get_network_flags(cfg)
    
    # Computation flags
    is_raster = false
    solver = cfg["solver"]

    # Output flags
    write_volt_maps = cfg["write_volt_maps"] in truelist
    write_cur_maps = cfg["write_cur_maps"] in truelist
    write_cum_cur_maps_only = cfg["write_cum_cur_map_only"] in truelist
    write_max_cur_maps = cfg["write_max_cur_maps"] in truelist
    set_null_currents_to_nodata = cfg["set_null_currents_to_nodata"] in truelist
    set_null_voltages_to_nodata = cfg["set_null_voltages_to_nodata"] in truelist
    compress_grids = cfg["compress_grids"] in truelist
    log_transform_maps = cfg["log_transform_maps"] in truelist

    o = OutputFlags(write_volt_maps, write_cur_maps,
                    write_cum_cur_maps_only, write_max_cur_maps,
                    set_null_currents_to_nodata, set_null_voltages_to_nodata,
                    compress_grids, log_transform_maps)
    
    NetworkFlags(is_raster, solver, o)
end

struct NetworkFlags
    is_raster::Bool
    solver::String
    outputflags::OutputFlags
end

