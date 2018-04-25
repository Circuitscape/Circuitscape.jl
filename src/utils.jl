export  model_problem,
        test_problem, 
        accumulate_current_maps,
        calculate_cum_current_maps,
        calculate_max_current_maps

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

# Helps start new processes from the INI file
function myaddprocs(n)
    addprocs(n)
    @everywhere eval(:(using Circuitscape))
end

# Reads the directory with the current maps 
# and accumulates all current maps
function accumulate_current_maps(path, f)
    dir = dirname(path)
    base = basename(path)
    
    # If base file has a dot 
    name = split(base, '.')[1]

    cmap_list = readdir(dir) |> 
                    x -> filter(y -> startswith(y, "$(name)_"), x) |>
                    x -> filter(y -> contains(y, "_curmap_"), x)
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
        csinfo("Accumulating $file")
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
    csinfo("Writing to $accum_path")
    open(accum_path, "w") do f
        write(f, headers)
        writedlm(f, round.(accum, 8), ' ')
    end

end

function f_in_place!(accum, cmap, f)
    accum .= f.(accum, cmap)
end

calculate_cum_current_map(path) = accumulate_current_maps(path, +)
calculate_max_current_map(path) = accumulate_current_maps(path, max)
