import Circuitscape: construct_node_map, compute_omniscape_current
using Circuitscape
import LinearAlgebra

# Omniscape moving window solve test 
# just checking syntax, other tests should be sufficient to ensure correctness.
let
        conductance = [
                1 5   1.;
                2 1   1;
                9 1   6
        ]
        source = [
                1 0 0.;
                0 0 0;
                0 1 0
        ]
        ground = [
                0 0 1.;
                0 0 0;
                0 0 0
        ]

        cs_cfg = Dict{String, String}()

        cs_cfg["ground_file_is_resistances"] = "True"
        cs_cfg["use_direct_grounds"] = "False"
        cs_cfg["output_file"] = "temp"
        cs_cfg["write_cum_cur_map_only"] = "False"
        cs_cfg["scenario"] = "Advanced"
        cs_cfg["suppress_messages"] = "True"
        cs_cfg["connect_four_neighbors_only"] = "False"
        cs_cfg["solver"] = "cholmod"
        cs_cfg["cholmod_batch_size"] = "1000"
        cs_cfg["data_type"] = "raster"

        current = compute_omniscape_current(
                conductance,
                source,
                ground,
                cs_cfg
        )
        
end
# Construct nodemap tests
let
        gmap = [0 1 2
                2 0 0
                2 0 2]
        nodemap = construct_node_map(gmap, Matrix{Int}(undef,0,0))
        @test nodemap == [0 3 4
                          1 0 0
                          2 0 5]
end

let
        gmap = [0 1 2
               2 0 0
               2 0 2]
        polymap = [1 0 1
                   2 1 0
                   0 0 2]
        nodemap = construct_node_map(gmap, polymap)
        @test nodemap == [4  3  4
                          1  4  0
                          2  0  1]
end

let
        gmap = [1 0 1
                0 1 0
                1 0 1]

        polymap = [1 0 1
                0 2 0
                2 0 0]

        r = construct_node_map(gmap, polymap)

        @test r == [1 0 1
                    0 2 0
                    2 0 3]
end

let
    polymap = [ 1.0  2.0  0.0  0.0  0.0
                0.0  0.0  0.0  0.0  0.0
                0.0  0.0  0.0  0.0  0.0
                0.0  0.0  0.0  0.0  0.0
                1.0  0.0  0.0  0.0  2.0]

    gmap = [0    0    0    1.0   1.0
            0    0    0    3.01  2.0
            1.0  2.0  2.0  1.0   1.0
            1.0  2.0  2.0  1.0   1.0
            1.0  2.0  2.0  0     1.0]

    nodemap = construct_node_map(gmap, polymap)

    @test nodemap == [ 3.0  18.0  0.0  10.0  14.0
                       0.0   0.0  0.0  11.0  15.0
                       1.0   4.0  7.0  12.0  16.0
                       2.0   5.0  8.0  13.0  17.0
                       3.0   6.0  9.0   0.0  18.0]
end

let

    println("pwd = $(pwd())")
    cfg = Circuitscape.parse_config("input/raster/one_to_all/11/oneToAllVerify11.ini")
    r = Circuitscape.load_raster_data(Float64, Int32, cfg)

    cellmap = r.cellmap
    polymap = r.polymap
    points_rc = r.points_rc
    point_map = [ 1.0  2.0  0.0  0.0  0.0
                  0.0  0.0  0.0  0.0  0.0
                  3.0  0.0  0.0  7.0  0.0
                  4.0  0.0  0.0  0.0  0.0
                  1.0  0.0  0.0  0.0  2.0 ]

    r = Circuitscape.create_new_polymap(cellmap, polymap, points_rc, 0, 0, point_map)

    @test r == [ 1.0  2.0  0.0  0.0  0.0
                 0.0  0.0  0.0  0.0  0.0
                 12.0  0.0  0.0  2.0  0.0
                 1.0  0.0  0.0  0.0  0.0
                 1.0  0.0  0.0  0.0  2.0 ]
end

import Circuitscape: resolve_conflicts

@test resolve_conflicts([1.,0.,0.], [1.,0.,0.], :rmvgnd) == ([1, 0, 0], [0, 0, 0], [1, 0, 0])
@test resolve_conflicts([1.,0.,0.], [1.,0.,0.], :rmvsrc) == ([0, 0, 0], [1, 0, 0], [1, 0, 0])
@test resolve_conflicts([1.,0.,0.], [1.,0.,0.], :keepall) == ([1, 0, 0], [1, 0, 0], [1, 0, 0])
@test resolve_conflicts([1.,0.,0.], [1.,0.,0.], :rmvall) == ([0, 0, 0], [1, 0, 0], [1, 0, 0])

# Construct graph
import Circuitscape: construct_graph
let
        gmap = Float64[0 1 2
                2 0 0
                2 0 2]
        nodemap = [0 3 4
                   1 0 0
                   2 0 5]
        A = construct_graph(gmap, nodemap, false, true)
        r = Matrix(A) - [0 2 0 0 0
                       2 0 0 0 0
                       0 0 0 1.5 0
                       0 0 1.5 0 0
                       0 0 0 0 0]
        @test sum(abs2, r) < 1e-6
        A = construct_graph(gmap, nodemap, true, true)
        r = Matrix(A) - [0 2 0 0 0
                       2 0 0 0 0
                       0 0 0 1.3333 0
                       0 0 1.33333 0 0
                       0 0 0 0 0]
        @test sum(abs2, r) < 1e-6
        A = construct_graph(gmap, nodemap, false, false)
        r = Matrix(A) - [0 2 1.06066 0 0
                       2 0 0 0 0
                       1.06066 0 0 1.5 0
                       0 0 1.5 0 0
                       0 0 0 0 0]
        @test sum(abs2, r) < 1e-6
        A = construct_graph(gmap, nodemap, true, false)
        r = Matrix(A) - [0 2 0.942809 0 0
                       2 0 0 0 0
                       0.942809 0 0 1.3333 0
                       0 0 1.3333 0 0
                       0 0 0 0 0]
        @test sum(abs2, r) < 1e-6

end

# 2D Model Problems

SIZE_2 =

[2.0  -1.0  -1.0   0.0
-1.0   2.0   0.0  -1.0
-1.0   0.0   2.0  -1.0
 0.0  -1.0  -1.0   2.0]

@test model_problem(2) == SIZE_2

SIZE_3 =

[2.0  -1.0   0.0  -1.0   0.0   0.0   0.0   0.0   0.0
-1.0   3.0  -1.0   0.0  -1.0   0.0   0.0   0.0   0.0
 0.0  -1.0   2.0   0.0   0.0  -1.0   0.0   0.0   0.0
-1.0   0.0   0.0   3.0  -1.0   0.0  -1.0   0.0   0.0
 0.0  -1.0   0.0  -1.0   4.0  -1.0   0.0  -1.0   0.0
 0.0   0.0  -1.0   0.0  -1.0   3.0   0.0   0.0  -1.0
 0.0   0.0   0.0  -1.0   0.0   0.0   2.0  -1.0   0.0
 0.0   0.0   0.0   0.0  -1.0   0.0  -1.0   3.0  -1.0
 0.0   0.0   0.0   0.0   0.0  -1.0   0.0  -1.0   2.0]

 @test model_problem(3) == SIZE_3

# Issue 151
# Issue 151
try
    Circuitscape.read_point_map(Int32, "samples.txt",
                                Circuitscape.RasterMeta(50, 50, 0.0, 0.0, 0.5, -9999.0, [0.0], ""))
catch e
    @test e == "At least one focal node location falls outside of habitat map"
end

# Users with dots in their names - issue #181
# Just check that this does not break
#compute("input/raster/extra.one/1/oneToAllVerify1.ini")

# Issue 158: cumulative current map should not be produced when write_cur_maps=False
let
    cum_file = "output/sgVerify12_cum_curmap.asc"
    isfile(cum_file) && rm(cum_file)
    compute("input/raster/pairwise/12/sgVerify12.ini")
    @test !isfile(cum_file)
end

# Configuration parsing and validation: typos and missing files must fail
# before any data is read, naming the INI key (#246, #341, #436).
@testset "Config validation" begin
    base = Circuitscape.init_config()
    base["scenario"] = "pairwise"

    # Unrecognised enum values are errors, not silent defaults
    for (k, v) in [("solver", "cholmodd"), ("scenario", "pairwsie"),
                   ("precision", "half"), ("write_cur_maps", "Ture"),
                   ("remove_src_or_gnd", "rmvboth"), ("log_level", "VERBOSE"),
                   ("data_type", "vector")]
        d = copy(base); d[k] = v
        @test_throws ArgumentError Circuitscape.CSConfig(d)
    end
    d = copy(base); d["scenario"] = "not entered"
    @test_throws ArgumentError Circuitscape.CSConfig(d)

    # Placeholders from init_config count as unset
    err = try Circuitscape.validate(Circuitscape.CSConfig(base)); nothing catch e; e end
    @test err isa ArgumentError
    @test occursin("habitat_file is not set", err.msg)
    @test occursin("point_file is not set", err.msg)
    @test occursin("output_file is not set", err.msg)

    good = copy(base)
    good["habitat_file"] = "input/raster/pairwise/1/cellmap.asc"
    good["point_file"] = "input/raster/pairwise/1/points.asc"
    good["output_file"] = "output/validate.out"
    @test Circuitscape.validate(Circuitscape.CSConfig(good)) isa Circuitscape.CSConfig

    d = copy(good); d["habitat_file"] = "missing.asc"
    err = try Circuitscape.validate(Circuitscape.CSConfig(d)); nothing catch e; e end
    @test occursin("habitat_file = \"missing.asc\" does not exist", err.msg)

    d = copy(good); d["output_file"] = "no_such_dir/x.out"
    @test_throws ArgumentError Circuitscape.validate(Circuitscape.CSConfig(d))

    # Only files that the run will use are checked
    d = copy(good); d["polygon_file"] = "missing.asc"; d["use_polygons"] = "False"
    @test Circuitscape.validate(Circuitscape.CSConfig(d)) isa Circuitscape.CSConfig
    d["use_polygons"] = "True"
    @test_throws ArgumentError Circuitscape.validate(Circuitscape.CSConfig(d))

    # Advanced mode needs sources and grounds, not a point file
    d = copy(good); d["scenario"] = "advanced"; d["point_file"] = "missing.asc"
    d["source_file"] = "input/raster/advanced/1/sources5x5.asc"
    d["ground_file"] = "input/raster/advanced/1/grounds5x5.asc"
    @test Circuitscape.validate(Circuitscape.CSConfig(d)) isa Circuitscape.CSConfig

    # Options accepted for INI compatibility but never implemented
    for k in ("low_memory_mode", "preemptive_memory_release")
        d = copy(good); d[k] = "True"
        err = try Circuitscape.validate(Circuitscape.CSConfig(d)); nothing catch e; e end
        @test occursin("$k is not implemented", err.msg)
    end

    # Unknown keys are ignored (with a warning), not fatal
    d = copy(good); d["max_parallel"] = "4"
    @test Circuitscape.CSConfig(d) isa Circuitscape.CSConfig
end

# Issue 342: set_focal_node_currents_to_zero
@testset "set_focal_node_currents_to_zero (#342)" begin
    dir = mktempdir()
    hdr = "ncols 5\nnrows 5\nxllcorner 0\nyllcorner 0\ncellsize 1\nNODATA_value -9999\n"
    write(joinpath(dir, "cell.asc"), hdr * join(fill(join(fill("1", 5), " "), 5), "\n") * "\n")
    # Node 1 is a two-cell region; nodes 2 and 3 are single cells.
    write(joinpath(dir, "pts.asc"), hdr *
          "1 1 0 0 2\n0 0 0 0 0\n0 0 0 0 0\n0 0 0 0 0\n3 0 0 0 0\n")
    focal = [(1,1), (1,2), (1,5), (5,1)]
    d = Circuitscape.init_config()
    d["data_type"] = "raster"; d["scenario"] = "pairwise"
    d["habitat_file"] = joinpath(dir, "cell.asc")
    d["point_file"] = joinpath(dir, "pts.asc")
    d["connect_four_neighbors_only"] = "True"
    d["write_cur_maps"] = "True"
    read_map(name) = readdlm(joinpath(dir, name), skipstart = 6)

    d["output_file"] = joinpath(dir, "plain.out")
    r_plain = compute(d)
    plain_pair = read_map("plain_curmap_1_2.asc")
    plain_cum = read_map("plain_cum_curmap.asc")
    @test all(plain_pair[i, j] > 0 for (i, j) in focal[1:3])

    d["set_focal_node_currents_to_zero"] = "True"
    d["output_file"] = joinpath(dir, "zeroed.out")
    r_zero = compute(d)
    zero_pair = read_map("zeroed_curmap_1_2.asc")
    zero_cum = read_map("zeroed_cum_curmap.asc")

    # Resistances are unaffected; only the maps change
    @test r_zero ≈ r_plain
    # The active pair (region 1, node 2) is zeroed; node 3 was not active
    @test all(zero_pair[i, j] == 0 for (i, j) in focal[1:3])
    @test zero_pair[5, 1] == plain_pair[5, 1]
    # Everything else is identical
    mask = trues(5, 5); for (i, j) in focal[1:3]; mask[i, j] = false; end
    @test zero_pair[mask] ≈ plain_pair[mask]
    # Every focal cell is active in some pair, so the cumulative map is
    # strictly smaller there and identical elsewhere
    @test all(zero_cum[i, j] < plain_cum[i, j] for (i, j) in focal)
    mask = trues(5, 5); for (i, j) in focal; mask[i, j] = false; end
    @test zero_cum[mask] ≈ plain_cum[mask]

    # One-to-all: all focal cells are zeroed in the cumulative map
    d["scenario"] = "one-to-all"
    d["output_file"] = joinpath(dir, "ota.out")
    compute(d)
    ota_cum = read_map("ota_cum_curmap.asc")
    @test all(ota_cum[i, j] == 0 for (i, j) in focal)
    @test any(ota_cum .> 0)
end

# Issue 470: the residual gate is configurable, and an unreachable tolerance
# fails with the retry history rather than a bare threshold.
@testset "Residual tolerance (#470)" begin
    d = Circuitscape.init_config(); d["scenario"] = "pairwise"
    cfg = Circuitscape.CSConfig(d)
    @test cfg.residual_tolerance == 0.0
    @test Circuitscape.residual_tolerance(cfg, Float64) == 1e-4
    @test Circuitscape.residual_tolerance(cfg, Float32) == 1e-3
    d["residual_tolerance"] = "1e-6"
    cfg = Circuitscape.CSConfig(d)
    @test Circuitscape.residual_tolerance(cfg, Float64) == 1e-6
    @test Circuitscape.residual_tolerance(cfg, Float32) == 1e-6
    @test Dict{String,String}(cfg)["residual_tolerance"] == "1.0e-6"

    G = model_problem(4)
    G.nzval .+= eps(Float64) * LinearAlgebra.norm(G.nzval)
    M = Circuitscape.aspreconditioner(Circuitscape.smoothed_aggregation(G))
    b = zeros(16); b[1] = -1; b[16] = 1
    v = Circuitscape.solve_linear_system(G, b, M)
    @test LinearAlgebra.norm(G * v .- b) / LinearAlgebra.norm(b) < 1e-4
    err = try Circuitscape.solve_linear_system(G, b, M; tol = 0.0); nothing catch e; e end
    @test err isa ErrorException
    @test occursin("after 3 attempts", err.msg)
end

# Issue 231: reclassify habitat raster values through a lookup table
@testset "Reclass table (#231)" begin
    # One pass: chained rows must not compound (1 -> 2 -> 3)
    m = [1.0 2.0; 2.0 3.0]
    Circuitscape.reclassify!(m, [1.0 2.0; 2.0 3.0])
    @test m == [2.0 3.0; 3.0 3.0]

    # A table that is not two columns is rejected by name
    dir = mktempdir()
    write(joinpath(dir, "bad.txt"), "1 2 3\n")
    @test_throws ArgumentError Circuitscape.read_reclass_table(Float64, joinpath(dir, "bad.txt"))

    # End to end: doubling every resistance doubles every effective resistance
    hdr = "ncols 5\nnrows 5\nxllcorner 0\nyllcorner 0\ncellsize 1\nNODATA_value -9999\n"
    write(joinpath(dir, "cell.asc"), hdr * join(fill(join(fill("2", 5), " "), 5), "\n") * "\n")
    write(joinpath(dir, "pts.asc"), hdr *
          "1 0 0 0 2\n0 0 0 0 0\n0 0 0 0 0\n0 0 0 0 0\n3 0 0 0 0\n")
    write(joinpath(dir, "reclass.txt"), "2\t4\n")
    d = Circuitscape.init_config()
    d["data_type"] = "raster"; d["scenario"] = "pairwise"
    d["habitat_file"] = joinpath(dir, "cell.asc")
    d["point_file"] = joinpath(dir, "pts.asc")
    d["output_file"] = joinpath(dir, "out.out")
    d["connect_four_neighbors_only"] = "True"
    base = compute(d)
    d["use_reclass_table"] = "True"
    d["reclass_file"] = joinpath(dir, "reclass.txt")
    reclassed = compute(d)
    @test reclassed[2:end, 2:end] ≈ 2 .* base[2:end, 2:end]
    # Cells the table does not mention are untouched: no-op table gives base back
    write(joinpath(dir, "reclass.txt"), "7\t9\n")
    @test compute(d) ≈ base
end

# Native ESRI ASCII grid I/O (no GDAL on the .asc path)
@testset "ASCII grid I/O" begin
    dir = mktempdir()
    a = [1.0 2.5 -9999.0; 4.0 5.0 6.0]
    tr = [10.0, 0.5, 0.0, 21.0, 0.0, -0.5]      # xll=10, cellsize=0.5, 2 rows -> yll=20

    # Round trip with and without a projection sidecar
    Circuitscape.write_raster(joinpath(dir, "r"), a, "PROJCS[\"x\"]", tr, "asc")
    b, wkt, tr2 = Circuitscape.read_raster(joinpath(dir, "r.asc"), Float64)
    @test b == a && tr2 == tr && wkt == "PROJCS[\"x\"]"
    @test isfile(joinpath(dir, "r.prj"))
    Circuitscape.write_raster(joinpath(dir, "q"), a, "", tr, "asc")
    @test !isfile(joinpath(dir, "q.prj"))
    @test Circuitscape.read_raster(joinpath(dir, "q.asc"), Float32)[1] == Float32.(a)
    @test occursin("cellsize", read(joinpath(dir, "q.asc"), String))

    # Header dialects GDAL accepts: center origin, dx/dy, no NODATA line,
    # tabs, CRLF, mixed case, custom nodata
    write(joinpath(dir, "c.asc"),
          "NCOLS 3\r\nNROWS 2\r\nxllcenter\t10.25\r\nyllcenter\t20.25\r\ncellsize 0.5\r\n1 2 3\r\n4 5 6\r\n")
    b, _, tr2 = Circuitscape.read_raster(joinpath(dir, "c.asc"), Float64)
    @test b == [1.0 2 3; 4 5 6] && tr2 == tr
    write(joinpath(dir, "d.asc"),
          "ncols 3\nnrows 2\nxllcorner 0\nyllcorner 0\ndx 1\ndy 2\nNODATA_value -1\n1 -1 3\n4 5 NaN\n")
    b, _, tr2 = Circuitscape.read_raster(joinpath(dir, "d.asc"), Float64)
    @test b == [1.0 -9999 3; 4 5 -9999] && tr2 == [0.0, 1.0, 0.0, 4.0, 0.0, -2.0]
    Circuitscape.write_raster(joinpath(dir, "e"), b, "", tr2, "asc")
    @test occursin("dx", read(joinpath(dir, "e.asc"), String))
    @test Circuitscape.read_raster(joinpath(dir, "e.asc"), Float64)[3] == tr2

    # Single-row grid and gzipped input
    write(joinpath(dir, "one.asc"), "ncols 2\nnrows 1\nxllcorner 0\nyllcorner 0\ncellsize 1\n7 8\n")
    @test Circuitscape.read_raster(joinpath(dir, "one.asc"), Float64)[1] == [7.0 8.0]
    Circuitscape.GZip.open(joinpath(dir, "z.asc.gz"), "w") do io
        write(io, read(joinpath(dir, "r.asc")))
    end
    @test Circuitscape.read_raster(joinpath(dir, "z.asc.gz"), Float64)[1] == a

    # Bad headers are reported by file
    write(joinpath(dir, "bad.asc"), "ncols 3\nnrows 3\nxllcorner 0\nyllcorner 0\ncellsize 1\n1 2 3\n4 5 6\n")
    @test_throws ErrorException Circuitscape.read_raster(joinpath(dir, "bad.asc"), Float64)

    # GDAL is optional; .asc never needs it, GeoTIFF does
    @test !Circuitscape.needs_gdal(joinpath(dir, "r.asc"))
    @test Circuitscape.needs_gdal("input/raster/pairwise/1/polygons.tif")
    @test occursin("Pkg.add(\"ArchGDAL\")", Circuitscape.gdal_missing_message("x"))
    @test Circuitscape.gdal_loaded() == GDAL_AVAILABLE

    if GDAL_AVAILABLE
        # Every .asc shipped with the tests reads the same as through GDAL.
        # GDAL's AAIGrid driver stores these values as Float32
        # (2.002 -> 2.0020000934); the native parser is exact, hence the
        # Float32-level tolerance.
        files = String[]
        for (root, _, fs) in walkdir("input"), f in fs
            endswith(f, ".asc") && push!(files, joinpath(root, f))
        end
        @test length(files) > 100
        for f in files
            a1, w1, t1 = Circuitscape.read_asc(f, Float64)
            a2, w2, t2 = Circuitscape.read_raster_gdal(f, Float64)
            @test isapprox(a1, a2; rtol = 1e-6) && t1 ≈ t2 && w1 == w2
        end
    else
        # Without ArchGDAL, a GeoTIFF is refused before any data is read
        d = Circuitscape.init_config()
        d["scenario"] = "pairwise"
        d["habitat_file"] = "input/raster/pairwise/1/cellmap.asc"
        d["point_file"] = "input/raster/pairwise/1/points.asc"
        d["use_polygons"] = "True"
        d["polygon_file"] = "input/raster/pairwise/1/polygons.tif"
        d["output_file"] = "output/gdal.out"
        err = try Circuitscape.validate(Circuitscape.CSConfig(d)); nothing catch e; e end
        @test err isa ArgumentError && occursin("polygon_file", err.msg) && occursin("ArchGDAL", err.msg)
        d["use_polygons"] = "False"; d["write_as_tif"] = "True"
        err = try Circuitscape.validate(Circuitscape.CSConfig(d)); nothing catch e; e end
        @test err isa ArgumentError && occursin("write_as_tif", err.msg)
    end
end
