# __precompile__(false)
module Circuitscape
using AlgebraicMultigrid
using ArchGDAL
using Graphs
using SimpleWeightedGraphs
using Krylov
using GZip

using LinearAlgebra
using SparseArrays
using SharedArrays
using Distributed
using DelimitedFiles
using Logging
using Dates
using SuiteSparse

include("consts.jl")
include("config.jl")
include("logging.jl")
include("utils.jl")
include("io.jl")
include("out.jl")
include("core.jl")
include("network/pairwise.jl")
include("raster/pairwise.jl")
include("raster/advanced.jl")
include("network/advanced.jl")
include("raster/onetoall.jl")
include("run.jl")
include("INIBuilder/INIBuilder.jl")

using .INIBuilder
export start

using PrecompileTools
@setup_workload begin
    workload_dir = mktempdir()
    # 5x5 raster grid
    cellmap = joinpath(workload_dir, "cellmap.asc")
    write(cellmap, """ncols         5
nrows         5
xllcorner     100
yllcorner     100
cellsize      10
NODATA_value  -9999
1\t1\t-9999\t1\t1
1\t-9999\t-9999\t-9999\t-9999
1\t-9999\t-9999\t1\t1
1\t-9999\t-9999\t1\t1
1\t-9999\t-9999\t1\t1
""")
    points = joinpath(workload_dir, "points.txt")
    write(points, """1\t105\t145
2\t115\t145
3\t105\t125
4\t105\t115
1\t105\t105
2\t145\t105
7\t135\t125
""")
    ini = joinpath(workload_dir, "job.ini")
    write(ini, """[Circuitscape mode]
data_type = raster
scenario = pairwise

[Habitat raster or graph]
habitat_file = $cellmap
habitat_map_is_resistances = True

[Options for pairwise and one-to-all and all-to-one modes]
point_file = $points
use_included_pairs = False
included_pairs_file = None

[Connection scheme for raster habitat data]
connect_four_neighbors_only = True
connect_using_avg_resistances = True

[Output options]
output_file = $(joinpath(workload_dir, "out.out"))
write_cur_maps = False
write_volt_maps = False
write_cum_cur_map_only = False
write_max_cur_maps = False
compress_grids = False
log_transform_maps = False
set_null_currents_to_nodata = False
set_null_voltages_to_nodata = False

[Calculation options]
solver = cg+amg
low_memory_mode = False

[Short circuit regions (aka polygons)]
use_polygons = False
polygon_file = None

[Options for advanced mode]
ground_file_is_resistances = False
source_file = None
remove_src_or_gnd = keepall
ground_file = None
use_unit_currents = False
use_direct_grounds = False

[Options for one-to-all and all-to-one modes]
use_variable_source_strengths = False
variable_source_file = None

[Mask file]
use_mask = False
mask_file = None

[Version]
version = unknown
""")
    @compile_workload begin
        Logging.disable_logging(Logging.Info)
        compute(ini)
        Logging.disable_logging(Logging.Debug)
    end
end

end
