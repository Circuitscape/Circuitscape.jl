"""
    `compute(path::String)`

Call the `compute` function on the configuration file. 

Inputs:
======

* path::String - Path to configuration file

"""
function compute(path::String)
    cfg = parse_config(path)
    data_type = cfg["data_type"]
    if data_type == "network"
        result = compute_network(cfg)
        return result
    else
        result = compute_raster(cfg)
        return result
    end
end
