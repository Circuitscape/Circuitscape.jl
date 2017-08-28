"""
    `compute(path::String)`

Call the `compute` function on the configuration file.

Inputs:
======

* path::String - Path to configuration file

"""
function compute(path::String)
    cfg = parse_config(path)
    T = parse_mode(cfg["data_type"], cfg["scenario"])
    compute(T, cfg)
end
function parse_mode(dt, scen)
    d = dt == "network" ? :Network : :Raster
    s = scen == "pairwise" ? :Pairwise : :Advanced
    eval(Expr(:call, Expr(:curly, d, s)))
end
