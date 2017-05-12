function parse_config(path::String)
    cf = init_config()
    f = open(path, "r")
    for i in EachLine(f)
        if first(i) == '['
            continue
        end
        idx = search(i, '=')
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

    a
end
