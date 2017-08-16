"""
    `compute(path::String)`

Call the `compute` function on the configuration file.

Inputs:
======

* path::String - Path to configuration file

"""
function compute(path::String)
    cfg = parse_config(path)
    data_type = cfg["data_type"] == "network" ? Network() : Raster()
    compute(data_type, cfg)
end
