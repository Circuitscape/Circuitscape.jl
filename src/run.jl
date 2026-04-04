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
    write_config(cfg)
    T = cfg.precision == pr_single ? Float32 : Float64
    if T == Float32 && (cfg.solver == st_cholmod || cfg.solver == st_pardiso)
        cswarn("Cholmod & Pardiso solver modes work only in double precision. Switching precision to double.")
        T = Float64
    end
    V = cfg.use_64bit_indexing ? Int64 : Int32
    csinfo("Precision used: $(_precision_str(cfg.precision))", cfg.suppress_messages)
    is_parallel = cfg.parallelize
    if is_parallel
        n = cfg.max_parallel
        csinfo("Starting up Circuitscape to use $n processes in parallel", cfg.suppress_messages)
        myaddprocs(n)
    end
    t = @elapsed r = _compute(T, V, cfg)

    csinfo("Time taken to complete job = $t", cfg.suppress_messages)
    is_parallel && rmprocs(workers())
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

function compute(dict)
    cfg_dict = init_config()
    update!(cfg_dict, dict)
    cfg = CSConfig(cfg_dict)
    update_logging!(cfg)
    write_config(cfg)
    T = cfg.precision == pr_single ? Float32 : Float64
    V = cfg.use_64bit_indexing ? Int64 : Int32
    csinfo("Precision used: $(_precision_str(cfg.precision))", cfg.suppress_messages)
    is_parallel = cfg.parallelize
    if is_parallel
        n = cfg.max_parallel
        csinfo("Starting up Circuitscape to use $n processes in parallel", cfg.suppress_messages)
        myaddprocs(n)
    end
    t = @elapsed r = _compute(T, V, cfg)

    csinfo("Time taken to complete job = $t", cfg.suppress_messages)
    is_parallel && rmprocs(workers())

    r
end
