@enum DataType dt_raster dt_network
@enum Scenario sc_pairwise sc_advanced sc_onetoall sc_alltoone
@enum SolverType st_cg_amg st_cholmod st_pardiso st_accelerate
@enum Precision pr_single pr_double
@enum LogLevel ll_none ll_info ll_debug ll_warning ll_critical
@enum RemovePolicy rp_keepall rp_rmvsrc rp_rmvgnd rp_rmvall

Base.@kwdef struct CSConfig
    version::String = "unknown"
    data_type::DataType = dt_raster
    scenario::Scenario = sc_pairwise
    habitat_file::String = ""
    habitat_map_is_resistances::Bool = true
    connect_four_neighbors_only::Bool = false
    connect_using_avg_resistances::Bool = false
    use_polygons::Bool = false
    polygon_file::String = ""
    source_file::String = ""
    ground_file::String = ""
    ground_file_is_resistances::Bool = true
    use_unit_currents::Bool = false
    use_direct_grounds::Bool = false
    remove_src_or_gnd::RemovePolicy = rp_keepall
    use_mask::Bool = false
    mask_file::String = ""
    solver::SolverType = st_cg_amg
    parallelize::Bool = false
    max_parallel::Int = 0
    precision::Precision = pr_double
    use_64bit_indexing::Bool = true
    cholmod_batch_size::Int = 1000
    low_memory_mode::Bool = false
    preemptive_memory_release::Bool = false
    use_variable_source_strengths::Bool = false
    variable_source_file::String = ""
    use_included_pairs::Bool = false
    included_pairs_file::String = ""
    point_file::String = ""
    use_reclass_table::Bool = false
    reclass_file::String = ""
    output_file::String = ""
    write_cur_maps::Bool = false
    write_volt_maps::Bool = false
    write_cum_cur_map_only::Bool = false
    write_max_cur_maps::Bool = false
    set_null_currents_to_nodata::Bool = false
    set_null_voltages_to_nodata::Bool = false
    set_focal_node_currents_to_zero::Bool = false
    compress_grids::Bool = false
    log_transform_maps::Bool = false
    write_as_tif::Bool = false
    log_file::String = ""
    log_level::LogLevel = ll_info
    suppress_messages::Bool = false
end

function _parse_bool(dict, key, default="false")
    get(dict, key, default) in ("True", "true", "1")
end

_parse_data_type(s) = s in RASTER ? dt_raster : dt_network

function _parse_scenario(s)
    s in PAIRWISE ? sc_pairwise :
    s in ADVANCED ? sc_advanced :
    s in ONETOALL ? sc_onetoall :
    s in ALLTOONE ? sc_alltoone : sc_pairwise
end

function _parse_solver(s)
    s in AMG ? st_cg_amg :
    s in CHOLMOD ? st_cholmod :
    s in PARDISO ? st_pardiso :
    s in ACCELERATE ? st_accelerate : st_cg_amg
end

_parse_precision(s) = s in SINGLE ? pr_single : pr_double

function _parse_log_level(s)
    s in NONE ? ll_none :
    s in INFO ? ll_info :
    s in DEBUG ? ll_debug :
    s in WARNING ? ll_warning :
    s in CRITICAL ? ll_critical : ll_info
end

function _parse_remove_policy(s)
    s == "rmvsrc" ? rp_rmvsrc :
    s == "rmvgnd" ? rp_rmvgnd :
    s == "rmvall" ? rp_rmvall : rp_keepall
end

function CSConfig(dict::Dict{String,String})
    CSConfig(
        version = get(dict, "version", "unknown"),
        data_type = _parse_data_type(get(dict, "data_type", "raster")),
        scenario = _parse_scenario(get(dict, "scenario", "not entered")),
        habitat_file = get(dict, "habitat_file", ""),
        habitat_map_is_resistances = _parse_bool(dict, "habitat_map_is_resistances", "True"),
        connect_four_neighbors_only = _parse_bool(dict, "connect_four_neighbors_only"),
        connect_using_avg_resistances = _parse_bool(dict, "connect_using_avg_resistances"),
        use_polygons = _parse_bool(dict, "use_polygons"),
        polygon_file = get(dict, "polygon_file", ""),
        source_file = get(dict, "source_file", ""),
        ground_file = get(dict, "ground_file", ""),
        ground_file_is_resistances = _parse_bool(dict, "ground_file_is_resistances", "True"),
        use_unit_currents = _parse_bool(dict, "use_unit_currents"),
        use_direct_grounds = _parse_bool(dict, "use_direct_grounds"),
        remove_src_or_gnd = _parse_remove_policy(get(dict, "remove_src_or_gnd", "keepall")),
        use_mask = _parse_bool(dict, "use_mask"),
        mask_file = get(dict, "mask_file", ""),
        solver = _parse_solver(get(dict, "solver", "cg+amg")),
        parallelize = _parse_bool(dict, "parallelize"),
        max_parallel = parse(Int, get(dict, "max_parallel", "0")),
        precision = _parse_precision(get(dict, "precision", "Double")),
        use_64bit_indexing = _parse_bool(dict, "use_64bit_indexing", "true"),
        cholmod_batch_size = parse(Int, get(dict, "cholmod_batch_size", "1000")),
        low_memory_mode = _parse_bool(dict, "low_memory_mode"),
        preemptive_memory_release = _parse_bool(dict, "preemptive_memory_release"),
        use_variable_source_strengths = _parse_bool(dict, "use_variable_source_strengths"),
        variable_source_file = get(dict, "variable_source_file", ""),
        use_included_pairs = _parse_bool(dict, "use_included_pairs"),
        included_pairs_file = get(dict, "included_pairs_file", ""),
        point_file = get(dict, "point_file", ""),
        use_reclass_table = _parse_bool(dict, "use_reclass_table"),
        reclass_file = get(dict, "reclass_file", ""),
        output_file = get(dict, "output_file", ""),
        write_cur_maps = _parse_bool(dict, "write_cur_maps"),
        write_volt_maps = _parse_bool(dict, "write_volt_maps"),
        write_cum_cur_map_only = _parse_bool(dict, "write_cum_cur_map_only"),
        write_max_cur_maps = _parse_bool(dict, "write_max_cur_maps"),
        set_null_currents_to_nodata = _parse_bool(dict, "set_null_currents_to_nodata"),
        set_null_voltages_to_nodata = _parse_bool(dict, "set_null_voltages_to_nodata"),
        set_focal_node_currents_to_zero = _parse_bool(dict, "set_focal_node_currents_to_zero"),
        compress_grids = _parse_bool(dict, "compress_grids"),
        log_transform_maps = _parse_bool(dict, "log_transform_maps"),
        write_as_tif = _parse_bool(dict, "write_as_tif"),
        log_file = let v = get(dict, "log_file", "None"); v == "None" ? "" : v end,
        log_level = _parse_log_level(get(dict, "log_level", "INFO")),
        suppress_messages = _parse_bool(dict, "suppress_messages"),
    )
end

# String converters for CSConfig fields
_bool_str(v::Bool) = v ? "True" : "False"

function _data_type_str(v::DataType)
    v == dt_raster ? "raster" : "network"
end

function _scenario_str(v::Scenario)
    v == sc_pairwise ? "pairwise" :
    v == sc_advanced ? "advanced" :
    v == sc_onetoall ? "one-to-all" :
    v == sc_alltoone ? "all-to-one" : "pairwise"
end

function _solver_str(v::SolverType)
    v == st_cg_amg ? "cg+amg" :
    v == st_cholmod ? "cholmod" :
    v == st_pardiso ? "mklpardiso" :
    v == st_accelerate ? "accelerate" : "cg+amg"
end

_precision_str(v::Precision) = v == pr_single ? "single" : "double"

function _log_level_str(v::LogLevel)
    v == ll_none ? "NONE" :
    v == ll_info ? "INFO" :
    v == ll_debug ? "DEBUG" :
    v == ll_warning ? "WARNING" :
    v == ll_critical ? "CRITICAL" : "INFO"
end

function _remove_policy_str(v::RemovePolicy)
    v == rp_keepall ? "keepall" :
    v == rp_rmvsrc ? "rmvsrc" :
    v == rp_rmvgnd ? "rmvgnd" :
    v == rp_rmvall ? "rmvall" : "keepall"
end

function _remove_policy_symbol(v::RemovePolicy)
    v == rp_keepall ? :keepall :
    v == rp_rmvsrc ? :rmvsrc :
    v == rp_rmvgnd ? :rmvgnd :
    v == rp_rmvall ? :rmvall : :keepall
end

function Base.Dict{String,String}(cfg::CSConfig)
    a = Dict{String,String}()
    a["version"] = cfg.version
    a["data_type"] = _data_type_str(cfg.data_type)
    a["scenario"] = _scenario_str(cfg.scenario)
    a["habitat_file"] = cfg.habitat_file
    a["habitat_map_is_resistances"] = _bool_str(cfg.habitat_map_is_resistances)
    a["connect_four_neighbors_only"] = _bool_str(cfg.connect_four_neighbors_only)
    a["connect_using_avg_resistances"] = _bool_str(cfg.connect_using_avg_resistances)
    a["use_polygons"] = _bool_str(cfg.use_polygons)
    a["polygon_file"] = cfg.polygon_file
    a["source_file"] = cfg.source_file
    a["ground_file"] = cfg.ground_file
    a["ground_file_is_resistances"] = _bool_str(cfg.ground_file_is_resistances)
    a["use_unit_currents"] = _bool_str(cfg.use_unit_currents)
    a["use_direct_grounds"] = _bool_str(cfg.use_direct_grounds)
    a["remove_src_or_gnd"] = _remove_policy_str(cfg.remove_src_or_gnd)
    a["use_mask"] = _bool_str(cfg.use_mask)
    a["mask_file"] = cfg.mask_file
    a["solver"] = _solver_str(cfg.solver)
    a["parallelize"] = _bool_str(cfg.parallelize)
    a["max_parallel"] = string(cfg.max_parallel)
    a["precision"] = _precision_str(cfg.precision)
    a["use_64bit_indexing"] = _bool_str(cfg.use_64bit_indexing)
    a["cholmod_batch_size"] = string(cfg.cholmod_batch_size)
    a["low_memory_mode"] = _bool_str(cfg.low_memory_mode)
    a["preemptive_memory_release"] = _bool_str(cfg.preemptive_memory_release)
    a["use_variable_source_strengths"] = _bool_str(cfg.use_variable_source_strengths)
    a["variable_source_file"] = cfg.variable_source_file
    a["use_included_pairs"] = _bool_str(cfg.use_included_pairs)
    a["included_pairs_file"] = cfg.included_pairs_file
    a["point_file"] = cfg.point_file
    a["use_reclass_table"] = _bool_str(cfg.use_reclass_table)
    a["reclass_file"] = cfg.reclass_file
    a["output_file"] = cfg.output_file
    a["write_cur_maps"] = _bool_str(cfg.write_cur_maps)
    a["write_volt_maps"] = _bool_str(cfg.write_volt_maps)
    a["write_cum_cur_map_only"] = _bool_str(cfg.write_cum_cur_map_only)
    a["write_max_cur_maps"] = _bool_str(cfg.write_max_cur_maps)
    a["set_null_currents_to_nodata"] = _bool_str(cfg.set_null_currents_to_nodata)
    a["set_null_voltages_to_nodata"] = _bool_str(cfg.set_null_voltages_to_nodata)
    a["set_focal_node_currents_to_zero"] = _bool_str(cfg.set_focal_node_currents_to_zero)
    a["compress_grids"] = _bool_str(cfg.compress_grids)
    a["log_transform_maps"] = _bool_str(cfg.log_transform_maps)
    a["write_as_tif"] = _bool_str(cfg.write_as_tif)
    a["log_file"] = cfg.log_file == "" ? "None" : cfg.log_file
    a["log_level"] = _log_level_str(cfg.log_level)
    a["suppress_messages"] = _bool_str(cfg.suppress_messages)
    a
end

function parse_config(path::String)
    cf = init_config()
    f = open(path, "r")
    for i in eachline(f, keep = true)
        if first(i) == '['
            continue
        end
        idx = something(findfirst(isequal('='), i), 0)
        var = rstrip(i[1:idx-1])
        val = strip(i[idx+1:end])
        cf[var] = val
    end
    close(f)
    CSConfig(cf)
end

# Defaults - kept for backward compatibility (returns Dict{String,String})
function init_config()
    a = Dict{String, String}()

    a["version"] = "unknown"
    a["connect_four_neighbors_only"] =  "False"
    a["connect_using_avg_resistances"] = "False"
    a["use_polygons"] =  "False"
    a["polygon_file"] = "(Browse for a short-circuit region file)"
    a["source_file"] = "(Browse for a current source file)"
    a["ground_file"]  = "(Browse for a ground point file)"
    a["ground_file_is_resistances"] = "True"
    a["use_unit_currents"] = "False"
    a["use_direct_grounds"] = "False"
    a["remove_src_or_gnd"]  = "keepall"
    a["mask_file"] =  "None"
    a["use_mask"] = "False"
    a["preemptive_memory_release"] = "False"
    a["low_memory_mode"] = "False"
    a["parallelize"] =  "False"
    a["max_parallel"] =  "0"
    a["print_timings"] = "False"
    a["print_rusages"] = "False"
    a["solver"] =  "cg+amg"
    a["use_variable_source_strengths"] = "False"
    a["variable_source_file"] =  "None"
    a["set_null_currents_to_nodata"] = "False"
    a["output_file"] = "(Choose a base name for output files)"
    a["write_cum_cur_map_only"] =  "False"
    a["log_transform_maps"] = "False"
    a["write_max_cur_maps"] = "False"
    a["compress_grids"] = "False"
    a["set_null_voltages_to_nodata"] = "False"
    a["set_focal_node_currents_to_zero"] = "False"
    a["write_volt_maps"] = "False"
    a["write_cur_maps"] = "False"
    a["habitat_map_is_resistances"] = "True"
    a["habitat_file"] = "(Browse for a resistance file)"
    a["scenario"] = "not entered"
    a["data_type"] = "raster"
    a["use_included_pairs"] = "False"
    a["included_pairs_file"] = "(Browse for a file with pairs to include or exclude)"
    a["point_file"] = "(Browse for file with locations of focal points or regions)"
    a["use_reclass_table"] = "False"
    a["reclass_file"] = "(Browse for file with reclassification data)"
    a["profiler_log_file"]  = "None"
    a["log_file"] = "None"
    a["log_level"] = "INFO"
    a["screenprint_log"] = "False"
    a["precision"] = "Double"
    a["log_file"] = "None"
    a["log_level"] = "INFO"
    a["cholmod_batch_size"] = "1000"
    a["use_64bit_indexing"] = "true"
    a["write_as_tif"] = "false"
    a["suppress_messages"] = "false"
    a
end

function update!(cfg, new)
    for (key,val) in new
        cfg[key] = val
    end
end

function write_config(cfg::CSConfig)
    open(cfg.output_file, "w") do f
        write(f, """
        [Circuitscape Mode]
        data_type = $(_data_type_str(cfg.data_type))
        scenario = $(_scenario_str(cfg.scenario))

        [Version]
        version = 5.0.0

        [Habitat raster or graph]
        habitat_file = $(cfg.habitat_file)
        habitat_map_is_resistances = $(cfg.habitat_map_is_resistances)

        [Connection Scheme for raster habitat data]
        connect_four_neighbors_only = $(cfg.connect_four_neighbors_only)
        connect_using_avg_resistances = $(cfg.connect_using_avg_resistances)

        [Short circuit regions (aka polygons)]
        use_polygons = $(cfg.use_polygons)
        polygon_file = $(cfg.polygon_file)

        [Options for advanced mode]
        ground_file_is_resistances = $(cfg.ground_file_is_resistances)
        source_file = $(cfg.source_file)
        remove_src_or_gnd = $(_remove_policy_str(cfg.remove_src_or_gnd))
        ground_file = $(cfg.ground_file)
        use_unit_currents = $(cfg.use_unit_currents)
        use_direct_grounds = $(cfg.use_direct_grounds)

        [Mask file]
        use_mask = $(cfg.use_mask)
        mask_file = $(cfg.mask_file)

        [Options for one-to-all and all-to-one modes]
        use_variable_source_strengths = $(cfg.use_variable_source_strengths)
        variable_source_file = $(cfg.variable_source_file)

        [Options for pairwise and one-to-all and all-to-one modes]
        included_pairs_file = $(cfg.included_pairs_file)
        use_included_pairs = $(cfg.use_included_pairs)
        point_file = $(cfg.point_file)

        [Calculation options]
        solver = $(_solver_str(cfg.solver))

        [Output options]
        write_cum_cur_map_only = $(cfg.write_cum_cur_map_only)
        log_transform_maps = $(cfg.log_transform_maps)
        output_file = $(cfg.output_file)
        write_max_cur_maps = $(cfg.write_max_cur_maps)
        write_volt_maps = $(cfg.write_volt_maps)
        set_null_currents_to_nodata = $(cfg.set_null_currents_to_nodata)
        set_null_voltages_to_nodata = $(cfg.set_null_voltages_to_nodata)
        compress_grids = $(cfg.compress_grids)
        write_cur_maps = $(cfg.write_cur_maps)
        """)
    end
end

# Keep backward compat for Dict-based write_config
function write_config(cfg::Dict{String,String})
    write_config(CSConfig(cfg))
end
