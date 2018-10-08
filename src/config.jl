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
    cf
end

# Defaults
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
    a["use_64bit_indexing"] = "false"

    a
end

function update!(cfg, new)
    for (key,val) in new
        cfg[key] = val
    end
end

function write_config(cfg)
    open(cfg["output_file"], "w") do f
        write(f, """
        [Circuitscape Mode]
        data_type = $(cfg["data_type"])
        scenario = $(cfg["scenario"])

        [Version]
        version = 5.0.0

        [Habitat raster or graph]
        habitat_file = $(cfg["habitat_file"])
        habitat_map_is_resistances = $(cfg["habitat_map_is_resistances"])

        [Connection Scheme for raster habitat data]
        connect_four_neighbors_only = $(cfg["connect_four_neighbors_only"] in TRUELIST)
        connect_using_avg_resistances = $(cfg["connect_using_avg_resistances"] in TRUELIST)

        [Short circuit regions (aka polygons)]
        use_polygons = $(cfg["use_polygons"] in TRUELIST)
        polygon_file = $(cfg["use_polygons"])

        [Options for advanced mode]
        ground_file_is_resistances = $(cfg["ground_file_is_resistances"] in TRUELIST)
        source_file = $(cfg["source_file"])
        remove_src_or_gnd = $(cfg["remove_src_or_gnd"])
        ground_file = $(cfg["ground_file"])
        use_unit_currents = $(cfg["use_unit_currents"] in TRUELIST)
        use_direct_grounds = $(cfg["use_direct_grounds"] in TRUELIST)

        [Mask file]
        use_mask = $(cfg["use_mask"] in TRUELIST)
        mask_file = $(cfg["mask_file"])

        [Options for one-to-all and all-to-one modes]
        use_variable_source_strengths = $(cfg["use_variable_source_strengths"] in TRUELIST)
        variable_source_file = $(cfg["variable_source_file"])

        [Options for pairwise and one-to-all and all-to-one modes]
        included_pairs_file = $(cfg["included_pairs_file"])
        use_included_pairs = $(cfg["use_included_pairs"] in TRUELIST)
        point_file = $(cfg["point_file"])

        [Calculation options]
        solver = cg+amg

        [Output options]
        write_cum_cur_map_only = $(cfg["write_cum_cur_map_only"] in TRUELIST)
        log_transform_maps = $(cfg["log_transform_maps"] in TRUELIST)
        output_file = $(cfg["output_file"])
        write_max_cur_maps = $(cfg["write_max_cur_maps"] in TRUELIST)
        write_volt_maps = $(cfg["write_volt_maps"] in TRUELIST)
        set_null_currents_to_nodata = $(cfg["set_null_currents_to_nodata"] in TRUELIST)
        set_null_voltages_to_nodata = $(cfg["set_null_voltages_to_nodata"] in TRUELIST)
        compress_grids = $(cfg["compress_grids"] in TRUELIST)
        write_cur_maps = $(cfg["write_cur_maps"] in TRUELIST)
        """)
    end
end
    
