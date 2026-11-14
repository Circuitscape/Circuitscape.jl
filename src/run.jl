
"""
    Circuitscape.compute(path::String) -> Matrix
    Circuitscape.compute(cfg::Dict{String,String}) -> Matrix

Run a Circuitscape job. `path` names an INI file; alternatively pass a
dictionary of INI keys and string values, as returned by
[`Circuitscape.init_config`](@ref), with the entries you need changed. Keys
that are not given take their default value.

The configuration is validated before any data is read (see
[`Circuitscape.validate`](@ref)), recorded in INI form in `<prefix>.ini`
(`output_file` without its extension; see [`Circuitscape.output_prefix`](@ref)),
and run. When `<prefix>.ini` is the INI file being run it is left as it is.
Output files are written next to `output_file` with that prefix; the return value is the
result matrix: pairwise resistances (first row and column hold the focal node
IDs), one-to-all resistances per focal node, node voltages in advanced mode.
Set `parallelize = True` and start Julia with several threads to run pair
solves in parallel.
"""
function compute(path::String)
    cfg = parse_config(path)
    _run(cfg, path)
end

function compute(dict)
    cfg_dict = init_config()
    update!(cfg_dict, dict)
    cfg = CSConfig(cfg_dict)
    _run(cfg)
end

function _run(cfg, input_ini = nothing)
    update_logging!(cfg)
    validate(cfg)
    record_config(cfg, input_ini)
    T = cfg.precision == pr_single ? Float32 : Float64
    if T == Float32 && cfg.solver == st_pardiso
        @warn("Pardiso solver works only in double precision. Switching precision to double.")
        T = Float64
    end
    @info("Precision used: $(_to_ini(cfg.precision))")
    if cfg.parallelize
        @info("Starting up Circuitscape to use $(Threads.nthreads()) threads in parallel")
    end
    with(CSTIMER => TimerOutput()) do
        r = @timeit CSTIMER[] "complete job" _compute(T, cfg)
        if cfg.log_level == Logging.Debug
            timings = CSTIMER[]
            @info("\n", timings)
        end
        r
    end
end

# Record the options of the run in `<prefix>.ini`, as Circuitscape 4 did.
# A common layout is `job.ini` with `output_file = job.out`: there the copy
# would replace the user's own file (comments, key order and any option the
# writer does not know about included), so the INI being run is left alone.
function record_config(cfg, input_ini)
    ini = output_prefix(cfg) * ".ini"
    if input_ini !== nothing && isfile(ini) && samefile(ini, input_ini)
        @info("Not writing the configuration to $ini: it is the INI file being run")
        return nothing
    end
    write_config(cfg, ini)
end

function _compute(T, cfg)
    data = @timeit CSTIMER[] "load data" load_data(T, cfg)
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
