# Issue 341: included pairs file should restrict which pairs are solved
# https://github.com/Circuitscape/Circuitscape.jl/issues/341

using Circuitscape, Test, DelimitedFiles

# Test 1: Simple synthetic - 3 focal points, include only pair (1,2)
let
    dir = mktempdir()
    write(joinpath(dir, "cell.asc"), """ncols         5
nrows         5
xllcorner     0
yllcorner     0
cellsize      1
NODATA_value  -9999
1 1 1 1 1
1 1 1 1 1
1 1 1 1 1
1 1 1 1 1
1 1 1 1 1
""")
    write(joinpath(dir, "pts.asc"), """ncols         5
nrows         5
xllcorner     0
yllcorner     0
cellsize      1
NODATA_value  -9999
1 0 0 0 2
0 0 0 0 0
0 0 0 0 0
0 0 0 0 0
3 0 0 0 0
""")
    write(joinpath(dir, "include.txt"), "mode\tinclude\n1\t2\n")
    ini = joinpath(dir, "job.ini")
    write(ini, """[Circuitscape mode]
data_type = raster
scenario = pairwise
[Habitat raster or graph]
habitat_file = $(joinpath(dir, "cell.asc"))
habitat_map_is_resistances = True
[Options for pairwise and one-to-all and all-to-one modes]
point_file = $(joinpath(dir, "pts.asc"))
use_included_pairs = True
included_pairs_file = $(joinpath(dir, "include.txt"))
[Connection scheme for raster habitat data]
connect_four_neighbors_only = True
connect_using_avg_resistances = True
[Output options]
output_file = $(joinpath(dir, "out.out"))
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
    r = compute(ini)
    # Only pair (1,2) included; point 3 is pruned.
    # Result: 3x3 (header + 2 points)
    @test size(r) == (3, 3)
    @test r[1, 2] == 1.0
    @test r[1, 3] == 2.0
    @test r[2, 3] > 0
end

# Test 2: Existing test case sgVerify17 — include pairs with focal points (txt list format)
let
    r = compute("input/raster/pairwise/17/sgVerify17.ini")
    x = readdlm("output_verify/sgVerify17_resistances.out")
    @test sum(abs2, x - r) < 1e-6
    # Verify excluded pairs are -1
    @test count(v -> v == -1, r[2:end, 2:end]) > 0
end

# Test 3: Existing test case sgVerify13 — include pairs with matrix format
let
    r = compute("input/raster/pairwise/13/sgVerify13.ini")
    x = readdlm("output_verify/sgVerify13_resistances.out")
    @test sum(abs2, x - r) < 1e-6
end

# Test 4: Focal regions (polygon path) — multiple cells per focal node
let
    dir = mktempdir()
    write(joinpath(dir, "cell.asc"), """ncols         6
nrows         6
xllcorner     0
yllcorner     0
cellsize      1
NODATA_value  -9999
1 1 1 1 1 1
1 1 1 1 1 1
1 1 1 1 1 1
1 1 1 1 1 1
1 1 1 1 1 1
1 1 1 1 1 1
""")
    # Focal regions: node 1 has 2 cells, node 2 has 2 cells, node 3 has 1 cell
    write(joinpath(dir, "pts.asc"), """ncols         6
nrows         6
xllcorner     0
yllcorner     0
cellsize      1
NODATA_value  -9999
1 1 0 0 2 2
0 0 0 0 0 0
0 0 0 0 0 0
0 0 0 0 0 0
0 0 0 0 0 0
3 0 0 0 0 0
""")
    write(joinpath(dir, "include.txt"), "mode\tinclude\n1\t2\n")
    ini = joinpath(dir, "job.ini")
    write(ini, """[Circuitscape mode]
data_type = raster
scenario = pairwise
[Habitat raster or graph]
habitat_file = $(joinpath(dir, "cell.asc"))
habitat_map_is_resistances = True
[Options for pairwise and one-to-all and all-to-one modes]
point_file = $(joinpath(dir, "pts.asc"))
use_included_pairs = True
included_pairs_file = $(joinpath(dir, "include.txt"))
[Connection scheme for raster habitat data]
connect_four_neighbors_only = True
connect_using_avg_resistances = True
[Output options]
output_file = $(joinpath(dir, "out.out"))
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
    r = compute(ini)
    # Point 3 pruned, result 3x3 (header + 2 nodes)
    @test size(r) == (3, 3)
    # Pair (1,2) solved
    @test r[2, 3] > 0
end

# Test 5: Exclude mode — exclude pair (1,3), solve (1,2) and (2,3)
let
    dir = mktempdir()
    write(joinpath(dir, "cell.asc"), """ncols         5
nrows         5
xllcorner     0
yllcorner     0
cellsize      1
NODATA_value  -9999
1 1 1 1 1
1 1 1 1 1
1 1 1 1 1
1 1 1 1 1
1 1 1 1 1
""")
    write(joinpath(dir, "pts.asc"), """ncols         5
nrows         5
xllcorner     0
yllcorner     0
cellsize      1
NODATA_value  -9999
1 0 0 0 2
0 0 0 0 0
0 0 0 0 0
0 0 0 0 0
3 0 0 0 0
""")
    write(joinpath(dir, "include.txt"), "mode\texclude\n1\t3\n")
    ini = joinpath(dir, "job.ini")
    write(ini, """[Circuitscape mode]
data_type = raster
scenario = pairwise
[Habitat raster or graph]
habitat_file = $(joinpath(dir, "cell.asc"))
habitat_map_is_resistances = True
[Options for pairwise and one-to-all and all-to-one modes]
point_file = $(joinpath(dir, "pts.asc"))
use_included_pairs = True
included_pairs_file = $(joinpath(dir, "include.txt"))
[Connection scheme for raster habitat data]
connect_four_neighbors_only = True
connect_using_avg_resistances = True
[Output options]
output_file = $(joinpath(dir, "out.out"))
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
    r = compute(ini)
    # All 3 nodes present, result 4x4 (header + 3)
    @test size(r) == (4, 4)
    # Pairs (1,2) and (2,3) should have resistance
    @test r[2, 3] > 0   # pair (1,2)
    @test r[3, 4] > 0   # pair (2,3)
    # Pair (1,3) should be excluded
    @test r[2, 4] == -1  # pair (1,3)
end

# Test 6: Exclude mode with multiple excluded pairs
let
    dir = mktempdir()
    write(joinpath(dir, "cell.asc"), """ncols         5
nrows         5
xllcorner     0
yllcorner     0
cellsize      1
NODATA_value  -9999
1 1 1 1 1
1 1 1 1 1
1 1 1 1 1
1 1 1 1 1
1 1 1 1 1
""")
    write(joinpath(dir, "pts.asc"), """ncols         5
nrows         5
xllcorner     0
yllcorner     0
cellsize      1
NODATA_value  -9999
1 0 0 0 2
0 0 0 0 0
0 0 0 0 0
0 0 0 0 0
3 0 0 0 4
""")
    # Exclude pairs (1,3) and (2,4) — only (1,2), (1,4), (2,3), (3,4) should solve
    write(joinpath(dir, "include.txt"), "mode\texclude\n1\t3\n2\t4\n")
    ini = joinpath(dir, "job.ini")
    write(ini, """[Circuitscape mode]
data_type = raster
scenario = pairwise
[Habitat raster or graph]
habitat_file = $(joinpath(dir, "cell.asc"))
habitat_map_is_resistances = True
[Options for pairwise and one-to-all and all-to-one modes]
point_file = $(joinpath(dir, "pts.asc"))
use_included_pairs = True
included_pairs_file = $(joinpath(dir, "include.txt"))
[Connection scheme for raster habitat data]
connect_four_neighbors_only = True
connect_using_avg_resistances = True
[Output options]
output_file = $(joinpath(dir, "out.out"))
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
    r = compute(ini)
    # All 4 nodes present, result 5x5 (header + 4)
    @test size(r) == (5, 5)
    # Solved pairs
    @test r[2, 3] > 0   # pair (1,2)
    @test r[2, 5] > 0   # pair (1,4)
    @test r[3, 4] > 0   # pair (2,3)
    @test r[4, 5] > 0   # pair (3,4)
    # Excluded pairs
    @test r[2, 4] == -1  # pair (1,3)
    @test r[3, 5] == -1  # pair (2,4)
end

# Test 7: Exclude mode with focal regions (polygon path)
let
    dir = mktempdir()
    write(joinpath(dir, "cell.asc"), """ncols         6
nrows         6
xllcorner     0
yllcorner     0
cellsize      1
NODATA_value  -9999
1 1 1 1 1 1
1 1 1 1 1 1
1 1 1 1 1 1
1 1 1 1 1 1
1 1 1 1 1 1
1 1 1 1 1 1
""")
    # Focal regions: node 1 has 2 cells, node 2 has 2 cells, node 3 has 1 cell
    write(joinpath(dir, "pts.asc"), """ncols         6
nrows         6
xllcorner     0
yllcorner     0
cellsize      1
NODATA_value  -9999
1 1 0 0 2 2
0 0 0 0 0 0
0 0 0 0 0 0
0 0 0 0 0 0
0 0 0 0 0 0
3 0 0 0 0 0
""")
    # Exclude pair (1,3) — pairs (1,2) and (2,3) should still solve
    write(joinpath(dir, "include.txt"), "mode\texclude\n1\t3\n")
    ini = joinpath(dir, "job.ini")
    write(ini, """[Circuitscape mode]
data_type = raster
scenario = pairwise
[Habitat raster or graph]
habitat_file = $(joinpath(dir, "cell.asc"))
habitat_map_is_resistances = True
[Options for pairwise and one-to-all and all-to-one modes]
point_file = $(joinpath(dir, "pts.asc"))
use_included_pairs = True
included_pairs_file = $(joinpath(dir, "include.txt"))
[Connection scheme for raster habitat data]
connect_four_neighbors_only = True
connect_using_avg_resistances = True
[Output options]
output_file = $(joinpath(dir, "out.out"))
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
    r = compute(ini)
    # All 3 nodes present, result 4x4 (header + 3)
    @test size(r) == (4, 4)
    # Solved pairs
    @test r[2, 3] > 0   # pair (1,2)
    @test r[3, 4] > 0   # pair (2,3)
    # Excluded pair
    @test r[2, 4] == -1  # pair (1,3)
end

# Test 8: The polygon path must count only the pairs it will solve, so the
# "Total number of pair solves" log reflects the include/exclude file rather
# than every combination of the nodes it mentions.
let
    pts = [1, 2, 3, 4]
    # Include only (1,2) and (1,3): everything else is excluded, both orders,
    # which is what generate_exclude_pairs produces in include mode.
    excluded = [(1,4), (4,1), (2,3), (3,2), (2,4), (4,2), (3,4), (4,3)]
    @test Circuitscape.pairs_to_solve(pts, excluded) == [(1,2), (1,3)]
    # No exclusions: all C(4,2) pairs, in row-major order
    @test Circuitscape.pairs_to_solve(pts, Tuple{Int,Int}[]) ==
        [(1,2), (1,3), (1,4), (2,3), (2,4), (3,4)]
    # A pair listed in only one order is still excluded
    @test Circuitscape.pairs_to_solve(pts, [(3,1)]) ==
        [(1,2), (1,4), (2,3), (2,4), (3,4)]
    # Node IDs need not be contiguous or sorted; indices are into `pts`
    @test Circuitscape.pairs_to_solve([10, 6, 1], [(6,10)]) == [(1,3), (2,3)]
end

# Shared tail for the raster ini files below.
const _INI_TAIL_341 = """
[Connection scheme for raster habitat data]
connect_four_neighbors_only = True
connect_using_avg_resistances = True
[Output options]
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
"""

_asc_hdr_341(n) = "ncols $n\nnrows $n\nxllcorner 0\nyllcorner 0\ncellsize 1\nNODATA_value -9999\n"

# Test 9: one-to-all must index the include matrix by node ID, not by position
# among the focal points present in the raster. The two differ as soon as the
# include file names an ID the point raster does not contain, which silently
# grounded the wrong nodes.
let
    dir = mktempdir()
    write(joinpath(dir, "cell.asc"),
          _asc_hdr_341(5) * join(fill(join(fill("1", 5), " "), 5), "\n") * "\n")
    # Raster contains nodes 2, 3, 4 only.
    write(joinpath(dir, "pts.asc"), _asc_hdr_341(5) *
          "2 0 0 0 3\n0 0 0 0 0\n0 0 0 0 0\n0 0 0 0 0\n4 0 0 0 0\n")
    # Node 1 is named in the file but absent from the raster, shifting every
    # include-matrix row by one relative to the surviving focal points.
    write(joinpath(dir, "include.txt"), "mode\tinclude\n1\t2\n3\t4\n")
    ini = joinpath(dir, "job.ini")
    write(ini, """[Circuitscape mode]
data_type = raster
scenario = one-to-all
[Habitat raster or graph]
habitat_file = $(joinpath(dir, "cell.asc"))
habitat_map_is_resistances = True
[Options for pairwise and one-to-all and all-to-one modes]
point_file = $(joinpath(dir, "pts.asc"))
use_included_pairs = True
included_pairs_file = $(joinpath(dir, "include.txt"))
output_file = $(joinpath(dir, "out.out"))
""" * _INI_TAIL_341)
    r = compute(ini)
    @test size(r) == (3, 2)
    @test r[:, 1] == [2.0, 3.0, 4.0]
    # Node 2's only listed partner (1) is absent, so it has nothing to ground to.
    @test r[1, 2] == -1
    # Nodes 3 and 4 are each other's only partner, so both see the same
    # single-source/single-ground resistance.
    @test r[2, 2] > 0
    @test r[2, 2] ≈ r[3, 2]
end

# Test 10: network mode must honour the include file too — it used to be read
# and then discarded, so every pair was solved.
let
    dir = mktempdir()
    # Path graph 1-2-3-4 with unit conductances.
    write(joinpath(dir, "graph.txt"), "1\t2\t1.0\n2\t3\t1.0\n3\t4\t1.0\n")
    write(joinpath(dir, "fp.txt"), "1\n2\n3\n4\n")
    write(joinpath(dir, "include.txt"), "mode\tinclude\n1\t2\n")
    ini = joinpath(dir, "job.ini")
    write(ini, """[Circuitscape mode]
data_type = network
scenario = pairwise
[Habitat raster or graph]
habitat_file = $(joinpath(dir, "graph.txt"))
habitat_map_is_resistances = False
[Options for pairwise and one-to-all and all-to-one modes]
point_file = $(joinpath(dir, "fp.txt"))
use_included_pairs = True
included_pairs_file = $(joinpath(dir, "include.txt"))
output_file = $(joinpath(dir, "out.out"))
""" * _INI_TAIL_341)
    r = compute(ini)
    # Nodes 3 and 4 are not named in the include file, so they are dropped.
    @test size(r) == (3, 3)
    @test r[1, 2] == 1.0
    @test r[1, 3] == 2.0
    # Adjacent nodes on a unit-conductance path graph.
    @test r[2, 3] ≈ 1.0
end

# Test 11: exclude_pairs_from derives the excluded node ID pairs from the matrix
let
    inc = Circuitscape.IncludeExcludePairs(:include, [1, 2, 3],
                                           [0 1 0; 1 0 0; 0 0 0])
    ex = Circuitscape.exclude_pairs_from(inc)
    # (1,2) is the only included pair, so every other combination is excluded.
    @test (1, 3) in ex && (3, 1) in ex
    @test (2, 3) in ex && (3, 2) in ex
    @test !((1, 2) in ex) && !((2, 1) in ex)

    exc = Circuitscape.IncludeExcludePairs(:exclude, [1, 2, 3],
                                           [0 1 0; 1 0 0; 0 0 0])
    ex2 = Circuitscape.exclude_pairs_from(exc)
    # Only the listed pair is excluded.
    @test (1, 2) in ex2 && (2, 1) in ex2
    @test !((1, 3) in ex2) && !((2, 3) in ex2)
end

# Test 12: network graphs numbered from 0 are shifted to 1-based on load, so the
# IDs in the include file must shift with them. Previously the unshifted file IDs
# were matched against shifted focal points, selecting the neighbouring pair.
let
    dir = mktempdir()
    # 0-based path graph with distinct resistances, so picking the wrong pair
    # gives a different number: 0-1 = 1, 1-2 = 5, 2-3 = 1.
    write(joinpath(dir, "graph.txt"), "0\t1\t1.0\n1\t2\t5.0\n2\t3\t1.0\n")
    write(joinpath(dir, "fp.txt"), "0\n1\n2\n3\n")
    write(joinpath(dir, "include.txt"), "mode\tinclude\n1\t2\n")
    ini = joinpath(dir, "job.ini")
    write(ini, """[Circuitscape mode]
data_type = network
scenario = pairwise
[Habitat raster or graph]
habitat_file = $(joinpath(dir, "graph.txt"))
habitat_map_is_resistances = True
[Options for pairwise and one-to-all and all-to-one modes]
point_file = $(joinpath(dir, "fp.txt"))
use_included_pairs = True
included_pairs_file = $(joinpath(dir, "include.txt"))
output_file = $(joinpath(dir, "out.out"))
""" * _INI_TAIL_341)
    r = compute(ini)
    @test size(r) == (3, 3)
    # Node IDs come back shifted to 1-based, as load_graph reports for 0-based input.
    @test r[1, 2] == 2.0
    @test r[1, 3] == 3.0
    # The requested pair is the 5.0 resistor, not the adjacent 1.0 one.
    @test r[2, 3] ≈ 5.0
end

