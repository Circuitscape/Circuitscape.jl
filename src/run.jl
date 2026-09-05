
"""
    Circuitscape.compute(path::String) -> Matrix
    Circuitscape.compute(cfg::Dict{String,String}) -> Matrix

Run a Circuitscape job. `path` names an INI file; alternatively pass a
dictionary of INI keys and string values, as returned by
[`Circuitscape.init_config`](@ref), with the entries you need changed. Keys
that are not given take their default value.

The configuration is validated before any data is read (see
[`Circuitscape.validate`](@ref)), written in INI form to `output_file`, and
run. Output files are written next to `output_file`; the return value is the
result matrix: pairwise resistances (first row and column hold the focal node
IDs), one-to-all resistances per focal node, node voltages in advanced mode.
Set `parallelize = True` and start Julia with several threads to run pair
solves in parallel.
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
    validate(cfg)
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
    with(CSTIMER => TimerOutput()) do
        r = @timeit CSTIMER[] "complete job" _compute(T, V, cfg)
        if cfg.log_level == Logging.Debug
            timings = CSTIMER[]
            @info("\n", timings)
        end
        r
    end
end

function _compute(T, V, cfg)
    data = @timeit CSTIMER[] "load data" load_data(T, V, cfg)
    if is_pairwise(cfg)
        run_pairwise(data, cfg)
    elseif is_advanced(cfg) || is_network(cfg)
        # A network has no one-to-all / all-to-one; any other scenario runs as
        # advanced, as it always has.
        run_advanced(build_problem(data, cfg), cfg)
    else
        run_onetoall(data, cfg)
    end
end
