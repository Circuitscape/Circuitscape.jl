export compute

"""
`compute(path::String)`

Call the `compute` function on the configuration file.

Inputs:
======

* path::String - Path to configuration file

"""
function compute(path::String)
    cfg = parse_config(path)
    update_logging!(cfg)
    T = cfg["precision"] in SINGLE ? Float32 : Float64
    csinfo("Precision used: $(cfg["precision"])")
    t = @elapsed r = _compute(T, cfg)
    csinfo("Time taken to complete job = $t")
    r
end

function _compute(T, cfg)
    is_raster = cfg["data_type"] == "raster"
    scenario = cfg["scenario"]
    if is_raster
        if scenario == "pairwise"
            raster_pairwise(T, cfg)
        elseif scenario == "advanced"
            raster_advanced(T, cfg)
        elseif scenario == "one-to-all"
            raster_one_to_all(T, cfg)
        else
            raster_one_to_all(T, cfg)
        end
    else
        if scenario == "pairwise"
            network_pairwise(T, cfg)
        else
            network_advanced(T, cfg)
        end
    end
end

function compute(dict)
    cfg = init_config()
    update!(cfg, dict)
    update_logging!(cfg)
    T = cfg["precision"] in SINGLE ? Float32 : Float64
    csinfo("Precision used: $(cfg["precision"])")
    t = @elapsed r = _compute(T, cfg)
    csinfo("Time taken to complete job = $t")
    r
end
