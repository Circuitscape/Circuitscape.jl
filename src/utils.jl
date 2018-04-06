export  model_problem,
        test_problem

 """
 Construct nodemap specific to a connected component
 """
function construct_local_node_map(nodemap, component, polymap)
    local_nodemap = zeros(eltype(nodemap), size(nodemap))
    idx = findin(nodemap, component)
    local_nodemap[idx] = nodemap[idx]
    if nodemap == local_nodemap
        return local_nodemap
    end
    _construct_local_nodemap(local_nodemap, polymap, idx)
end

function _construct_local_nodemap(local_nodemap, polymap, idx)
    if isempty(polymap)
        i = find(local_nodemap)
        local_nodemap[i] = 1:length(i)
        return local_nodemap
    else
        local_polymap = zeros(eltype(local_nodemap), size(local_nodemap))
        local_polymap[idx] = polymap[idx]
        return construct_node_map(local_nodemap, local_polymap)
    end
end

"""
Define model circuitscape problem - helps in tests
"""
function model_problem(T, s)
    
    # Cell map is the uniform 1 across the grid
    cellmap = ones(T, s, s)
    
    # Nodemap is just 1:endof(cellmap)
    nodemap = reshape(1:endof(cellmap), s, s) |> collect

    # Construct graph with 4 neighbors
    G = construct_graph(cellmap, nodemap, true, true)

    # Return laplacian
    laplacian(G)
end
model_problem(s::Integer) = model_problem(Float64, s)

# Testing utilities 

function test_problem(str)
    base_path = joinpath(Pkg.dir("Circuitscape"), "test", "input")
    str2 = replace(str, ".ini", "")
    if contains(str, "sgVerify")
        config_path = joinpath(base_path, "raster", "pairwise", replace(str2, "sgVerify", ""), str)
    elseif contains(str, "mgVerify")
        config_path = joinpath(base_path, "raster", "advanced", replace(str2, "mgVerify", ""), str)
    elseif contains(str, "oneToAll")
        config_path = joinpath(base_path, "raster", "one_to_all", replace(str2, "oneToAllVerify", ""), str)
    elseif contains(str, "allToOne")
        config_path = joinpath(base_path, "raster", "all_to_one", replace(str2, "allToOneVerify", ""), str)
    elseif contains(str, "sgNetworkVerify")
        config_path = joinpath(base_path, "network", str)
    elseif contains(str, "mgNetworkVerify")
        config_path = joinpath(base_path, "network", str)
    end
end    

# Utility funciton for calling Circuitscape with Cholmod
function compute_cholmod(str)
    cfg = parse_config(str)
    T = cfg["precision"] in SINGLE ? Float32 : Float64
    if T == Float32
        cswarn("Cholmod supports only double precision. Implicit conversion may occur")
    end
    cfg["solver"] = "cholmod"
    _compute(T, cfg)
end

function compute_single(str)
    cfg = parse_config(str)
    cfg["precision"] = "single"
    _compute(Float32, cfg)
end

function get_output_flags(cfg)

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
end
