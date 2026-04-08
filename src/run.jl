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
    _run(cfg)
end

function compute(dict)
    cfg_dict = init_config()
    update!(cfg_dict, dict)
    cfg = CSConfig(cfg_dict)
    _run(cfg)
end

function _run(cfg)
    update_logging!(cfg)
    write_config(cfg)
    T = cfg.precision == pr_single ? Float32 : Float64
    if T == Float32 && cfg.solver == st_pardiso
        @warn("Pardiso solver works only in double precision. Switching precision to double.")
        T = Float64
    end
    V = cfg.use_64bit_indexing ? Int64 : Int32
    @info("Precision used: $(_precision_str(cfg.precision))")
    if cfg.parallelize
        @info("Starting up Circuitscape to use $(Threads.nthreads()) threads in parallel")
    end
    reset_timer!(CSTIMER)
    r = @timeit CSTIMER "complete job" _compute(T, V, cfg)
    if cfg.log_level == Logging.Debug
        @info("\n", CSTIMER)
    end
    r
end

function _compute(T, V, cfg)
    is_raster = cfg.data_type == dt_raster
    scenario = cfg.scenario
    if is_raster
        if scenario == sc_pairwise
            raster_pairwise(T, V, cfg)
        elseif scenario == sc_advanced
            raster_advanced(T, V, cfg)
        elseif scenario == sc_onetoall
            raster_one_to_all(T, V, cfg)
        else
            raster_one_to_all(T, V, cfg)
        end
    else
        if scenario == sc_pairwise
            network_pairwise(T, V, cfg)
        else
            network_advanced(T, V, cfg)
        end
    end
end
