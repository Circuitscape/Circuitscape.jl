"""
    `compute(path::String)`

Call the `compute` function on the configuration file. 

Inputs:
======

* path::String - Path to configuration file

"""
function compute(path::String)
    cfg = ConfigFile(path)
    data_type = get(cfg, "Circuitscape mode", "data_type")
    if data_type == "network"
        result = compute_network(cfg)
        return result
    else
        result = compute_raster(cfg)
        return result
    end
end
