using IniFile

immutable Config
    version::VersionNumber
    connect_four_neighbors_only::Bool
    connect_using_avg_resistances::Bool
    use_polygons::Bool
    polygon_file::String
    source_file::String
    ground_file::String
    ground_file_is_resistances::Bool
    use_unit_currents::Bool
    use_direct_grounds::Bool
    remove_src_or_gnd::Symbol
    mask_file::String
    use_mask::Bool
    preemptive_memory_release::Bool
    low_memory_mode::Bool
    parallelize::Bool
    max_parallel::Int
    print_timings::Bool
    print_rusages::Bool
    solver::Symbol
    use_variable_source_strengths::Bool
    variable_source_file::String
    set_null_currents_to_nodata::Bool
    output_file::String
    write_cum_cur_map_only::Bool
    log_transform_maps::Bool
    write_max_cur_maps::Bool
    compress_grids::Bool
    set_null_voltages_to_nodata::Bool
    set_focal_node_currents_to_zero::Bool
    write_volt_maps::Bool
    write_cur_maps::Bool
    habitat_map_is_resistances::Bool
    habitat_file::String
    scenario::Symbol
    data_type::Symbol
    use_included_pairs::Bool
    included_pairs_file::String
    point_file::String
    use_reclass_table::Bool
    reclass_file::String
    profiler_log_file::String
    log_file::String
    log_level::Symbol
    screenprint_log::Bool
end 

function Config(;version = v"0", connect_four_neighbors_only = false, 
                connect_using_avg_resistances =  false, use_polygons = false,
            polygon_file = "", source_file = "", ground_file = "", 
            ground_file_is_resistances = true, use_unit_currents = false, 
            use_direct_grounds = false, remove_src_or_gnd = :keepall,
            mask_file = "", use_mask = false, preemptive_memory_release = false,
            low_memory_mode = false, parallelize = false, max_parallel = 0,
            print_timings = false, print_rusages =  false, solver =  :cgamg,
            use_variable_source_strengths = false, variable_source_file = "",
            set_null_currents_to_nodata = false, output_file = "",
            write_cum_cur_map_only = false, log_transform_maps =  false, 
            write_max_cur_maps = false, compress_grids = false, 
            set_null_voltages_to_nodata = false, set_focal_node_currents_to_zero = false, 
            write_volt_maps = false, write_cur_maps = false, 
            habitat_map_is_resistances = true, habitat_file = "",
            scenario = :null, data_type  = :raster, use_included_pairs = false, 
            included_pairs_file = "", point_file = "",
            use_reclass_table = false, reclass_file = "", profiler_log_file = "",
            log_file = "", log_level = :info, screenprint_log = false)
            
    Config(version, connect_four_neighbors_only, connect_using_avg_resistances, 
        use_polygons, polygon_file, source_file, ground_file, ground_file_is_resistances, 
        use_unit_currents, use_direct_grounds, remove_src_or_gnd, mask_file, 
        use_mask, preemptive_memory_release, low_memory_mode, parallelize, 
        max_parallel, print_timings, print_rusages, solver, use_variable_source_strengths, 
        variable_source_file, set_null_currents_to_nodata, output_file,
        write_cum_cur_map_only, log_transform_maps, write_max_cur_maps, 
        compress_grids, set_null_voltages_to_nodata, set_focal_node_currents_to_zero,
        write_volt_maps, write_cur_maps, habitat_map_is_resistances, habitat_file,
        scenario, data_type, use_included_pairs, included_pairs_file, point_file,
        use_reclass_table, reclass_file, profiler_log_file, log_file, log_level, 
        screenprint_log)

end

function ConfigFile(path::String)
    a = Inifile()
    read(a, path)
    a
end
