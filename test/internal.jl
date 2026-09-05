import Circuitscape: construct_node_map, compute_omniscape_current
using Circuitscape
import LinearAlgebra
import Graphs, SparseArrays, Random

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

# construct_graph must return exactly the matrix the original push!/`a + a'`
# implementation did (same sparsity pattern, same values). The old version is
# kept here as the reference.
let
    function construct_graph_ref(gmap, nodemap::Matrix{S}, avg_res, four_neighbors) where S
        f1 = avg_res ? Circuitscape.res_avg : Circuitscape.cond_avg
        f2 = avg_res ? Circuitscape.weirder_avg : Circuitscape.weird_avg
        I = Vector{S}()
        J = Vector{S}()
        V = Vector{eltype(gmap)}()
        for j = 1:size(gmap, 2)
            for i = 1:size(gmap, 1)
                if nodemap[i,j] == 0
                    continue
                else
                    if j != size(gmap, 2) && nodemap[i,j+1] != 0
                        push!(I, nodemap[i,j]); push!(J, nodemap[i,j+1])
                        push!(V, f1(gmap[i,j], gmap[i,j+1]))
                    end
                    if i != size(gmap, 1) && nodemap[i+1, j] != 0
                        push!(I, nodemap[i,j]); push!(J, nodemap[i+1,j])
                        push!(V, f1(gmap[i,j], gmap[i+1,j]))
                    end
                    if !four_neighbors
                        if i != size(gmap, 1) && j != size(gmap, 2) && nodemap[i+1, j+1] != 0
                            push!(I, nodemap[i,j]); push!(J, nodemap[i+1,j+1])
                            push!(V, f2(gmap[i,j], gmap[i+1,j+1]))
                        end
                        if i != 1 && j != size(gmap, 2) && nodemap[i-1, j+1] != 0
                            push!(I, nodemap[i,j]); push!(J, nodemap[i-1,j+1])
                            push!(V, f2(gmap[i,j], gmap[i-1,j+1]))
                        end
                    end
                end
            end
        end
        m = maximum(nodemap)
        a = SparseArrays.sparse(I,J,V, m, m)
        a = a + a'
        a
    end

    # Number nodes column-major over the non-zero cells of a mask
    function mask_nodemap(mask)
        nodemap = zeros(Int, size(mask))
        k = 0
        for j in 1:size(mask, 2), i in 1:size(mask, 1)
            if mask[i,j]
                k += 1
                nodemap[i,j] = k
            end
        end
        nodemap
    end

    function check_same(gmap, nodemap)
        for avg_res in (true, false), four_neighbors in (true, false)
            A = Circuitscape.construct_graph(gmap, nodemap, avg_res, four_neighbors)
            R = construct_graph_ref(gmap, nodemap, avg_res, four_neighbors)
            @test A isa SparseArrays.SparseMatrixCSC
            @test size(A) == size(R)
            @test SparseArrays.nnz(A) == SparseArrays.nnz(R)
            @test A.colptr == R.colptr
            @test A.rowval == R.rowval
            @test A.nzval == R.nzval
            @test A == R
            @test LinearAlgebra.issymmetric(A)
        end
    end

    # 5x5 all-ones grid
    check_same(ones(5, 5), mask_nodemap(trues(5, 5)))

    # Non-uniform resistances, so that f1/f2 are exercised with distinct values
    check_same(Float64[1 + (i * j) % 7 for i in 1:5, j in 1:5], mask_nodemap(trues(5, 5)))

    # NODATA holes in the interior and on the edges
    mask = trues(6, 7)
    mask[3, 4] = false   # interior
    mask[4, 4] = false   # adjacent interior
    mask[1, 2] = false   # top edge
    mask[6, 7] = false   # bottom-right corner
    mask[3, 1] = false   # left edge
    mask[1, 7] = false   # top-right corner
    gmap = Float64[1 + (i + 2j) % 5 for i in 1:6, j in 1:7]
    gmap[.!mask] .= 0
    check_same(gmap, mask_nodemap(mask))

    # 1xN and Nx1 grids
    check_same(ones(1, 9), mask_nodemap(trues(1, 9)))
    check_same(ones(9, 1), mask_nodemap(trues(9, 1)))
    check_same(Float64[1 2 3 4 5 6 7], mask_nodemap(trues(1, 7)))
    check_same(Float64[1 2 3 4 5 6 7]', mask_nodemap(trues(7, 1)))

    # Polygons: several adjacent cells share one node, producing duplicate
    # (i, j) pairs and self-loops which `sparse` must sum exactly as `a + a'` did
    nodemap = [1 1 2
               1 3 2
               4 3 5]
    check_same(Float64[1 2 3; 4 5 6; 7 8 9], nodemap)

    # Original test fixture
    check_same(Float64[0 1 2; 2 0 0; 2 0 2], [0 3 4; 1 0 0; 2 0 5])
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

    # CSConfig is the only representation past parsing: predicates read it
    # directly and a copy-with constructor replaces selected fields.
    cfg = Circuitscape.CSConfig(good)
    @test Circuitscape.is_raster(cfg) && Circuitscape.is_pairwise(cfg)
    @test !Circuitscape.is_advanced(cfg) && !Circuitscape.is_onetoall(cfg)
    @test Circuitscape.policy(cfg) == :keepall
    cfg2 = Circuitscape.CSConfig(cfg; scenario = Circuitscape.sc_advanced,
                                 remove_src_or_gnd = Circuitscape.rp_rmvsrc)
    @test Circuitscape.is_advanced(cfg2) && Circuitscape.policy(cfg2) == :rmvsrc
    @test cfg2.habitat_file == cfg.habitat_file && cfg2.output_file == cfg.output_file
    @test Dict{String,String}(cfg2)["scenario"] == "advanced"
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

# Connected components on the CSC structure must match what Graphs returned:
# same components, same order (by smallest vertex, ascending within).
@testset "connected_components" begin
    ref(A) = Graphs.connected_components(Graphs.SimpleGraph(A))
    same(A) = map(c -> Int.(c), Circuitscape.connected_components(A)) == ref(A)

    # Grid Laplacians (what the solvers pass), including a disconnected one
    @test same(model_problem(4))
    G = model_problem(5); G[13, :] .= 0; G[:, 13] .= 0; SparseArrays.dropzeros!(G)
    @test same(G)   # isolated vertex is its own component
    @test length(Circuitscape.connected_components(G)) == 2

    # Random symmetric graphs: many components, isolated vertices, Int32 indices
    Random.seed!(7)
    for trial in 1:200
        n = rand(1:60)
        B = SparseArrays.sprand(n, n, rand() * 0.15)
        A = B + B'
        @test same(A)
        @test same(convert(SparseArrays.SparseMatrixCSC{Float64,Int32}, A))
    end

    # Stored zeros are not edges, exactly as SimpleGraph treats them
    A = SparseArrays.sparse([1, 2, 2, 3], [2, 1, 3, 2], [0.0, 0.0, 1.0, 1.0], 3, 3)
    @test SparseArrays.nnz(A) == 4
    @test same(A) && length(Circuitscape.connected_components(A)) == 2

    @test_throws ArgumentError Circuitscape.connected_components(SparseArrays.sprand(3, 4, 0.5))
end

# The INI builder is an extension on the REPL stdlib
@testset "INI builder extension" begin
    # start() loads REPL itself, so a script need not; afterwards the
    # extension is active and the zero-argument method exists
    ext = Circuitscape.load_repl_extension()
    @test ext !== nothing
    @test Base.get_extension(Circuitscape, :CircuitscapeREPLExt) === ext
    @test hasmethod(Circuitscape.start, Tuple{})
    @test_throws MethodError Circuitscape.start(1)
end

# Network cumulative currents: the branch index must resolve either
# orientation of an edge to one slot, keep the first slot for a duplicated
# edge, and the accumulated result must not depend on how edges are listed.
@testset "network cumulative branch index" begin
    coords = ([1, 3, 2, 3, 4], [2, 2, 1, 4, 3], [1.0, 1.0, 1.0, 1.0, 1.0])  # (2,1) duplicates (1,2); (4,3) duplicates (3,4)
    cum = Circuitscape.initialize_cum_vectors(coords, 4)
    @test cum.branch_index[(1, 2)] == 1 && cum.branch_index[(2, 1)] == 1
    @test cum.branch_index[(3, 2)] == 2 && cum.branch_index[(2, 3)] == 2
    @test cum.branch_index[(3, 4)] == 4 && cum.branch_index[(4, 3)] == 4
    @test length(cum.cum_branch_curr) == 5

    # Path 1-2-3-4-5 with rows deliberately listed in mixed orientation. One
    # focal pair, so the cumulative branch currents equal that pair's own.
    dir = mktempdir()
    write(joinpath(dir, "g.txt"), "2\t1\t1.0\n2\t3\t2.0\n4\t3\t1.0\n4\t5\t3.0\n")
    write(joinpath(dir, "fp.txt"), "1\n5\n")
    d = Circuitscape.init_config()
    d["data_type"] = "network"; d["scenario"] = "pairwise"
    d["habitat_file"] = joinpath(dir, "g.txt"); d["habitat_map_is_resistances"] = "True"
    d["point_file"] = joinpath(dir, "fp.txt"); d["output_file"] = joinpath(dir, "o.out")
    d["write_cur_maps"] = "True"
    r = compute(d)
    @test r[2, 3] ≈ 7.0   # series resistances 1 + 2 + 1 + 3
    keyed(m) = Dict(minmax(Int(m[i, 1]), Int(m[i, 2])) => m[i, 3] for i in axes(m, 1))
    pair = keyed(readdlm(joinpath(dir, "o_branch_currents_1_5.txt")))
    cum = keyed(readdlm(joinpath(dir, "o_branch_currents_cum.txt")))
    @test length(cum) == 4 && keys(cum) == keys(pair)
    @test all(cum[k] ≈ pair[k] for k in keys(cum))
    @test all(v ≈ 1.0 for v in values(cum))   # unit current through every series branch
end

# The pairwise driver enumerates pairs with O(log n) lookups into the sorted
# component and a node -> positions map. The exact list of jobs, in order, is
# what the log numbering and the gold resistance files depend on.
@testset "pair_jobs" begin
    # Positions 1 and 3 of `points` (user ids 1 and 5) both map to node 10;
    # user ids 2 and 3 (nodes 20 and 30) are an excluded pair; node 50 in the
    # component has no focal point; node 99 is a focal point outside it.
    points   = [10, 20, 10, 30, 40, 99]
    orig_pts = [1, 2, 5, 3, 4, 9]
    comp     = [5, 10, 20, 30, 40, 50]
    exclude  = [(2, 3)]
    csub     = Circuitscape.component_points(comp, points)
    @test csub == [10, 20, 30, 40]

    positions = Circuitscape.point_positions(points)
    @test positions == Dict(10 => [1, 3], 20 => [2], 30 => [4], 40 => [5], 99 => [6])

    tup(j) = (j.src_node, j.dst_node, j.comp_i, j.comp_j, j.src_indices, j.dst_indices)
    groups = Circuitscape.pair_jobs(csub, comp, points, orig_pts, exclude, false)
    @test map(g -> map(tup, g), groups) == [
        [(10, 20, 2, 3, [1, 3], [2]),
         (10, 30, 2, 4, [1, 3], [4]),
         (10, 40, 2, 5, [1, 3], [5])],
        [(20, 40, 3, 5, [2], [5])],           # (20, 30) excluded via user ids (2, 3)
        [(30, 40, 4, 5, [4], [5])],
        [],                                   # last source: no destinations
    ]
    @test eltype(groups) == Vector{Circuitscape.PairJob{Int}}

    # Resistance shortcut: only the first source, still one group
    groups = Circuitscape.pair_jobs(csub, comp, points, orig_pts, exclude, true)
    @test length(groups) == 1
    @test map(tup, groups[1]) == [(10, 20, 2, 3, [1, 3], [2]),
                                  (10, 30, 2, 4, [1, 3], [4]),
                                  (10, 40, 2, 5, [1, 3], [5])]

    # A shared node keeps the pair while one index combination is not
    # excluded: (1, 2) is, (5, 2) is not
    groups = Circuitscape.pair_jobs([10, 20], comp, points, orig_pts, [(1, 2)], false)
    @test map(tup, groups[1]) == [(10, 20, 2, 3, [1, 3], [2])]

    # Excluding every combination drops the pair
    groups = Circuitscape.pair_jobs([20, 30], comp, points, orig_pts, [(2, 3)], false)
    @test groups == [Circuitscape.PairJob{Int}[], Circuitscape.PairJob{Int}[]]

    # Lookups agree with the linear scans they replace, on Int32 too
    comp32 = Int32[3, 7, 8, 12]
    @test Circuitscape.component_index(comp32, Int32(8)) == findfirst(isequal(8), comp32)
    @test Circuitscape.component_index(comp32, Int32(9)) === nothing
    @test Circuitscape.component_index(comp32, Int32(13)) === nothing
    @test Circuitscape.in_component(comp32, Int32(12)) && !Circuitscape.in_component(comp32, Int32(1))

    # A focal node missing from the component is an error, as before
    @test_throws ErrorException Circuitscape.pair_jobs([10, 99], comp, points, orig_pts, exclude, false)
    # The sorted-component contract is checked, not assumed
    @test_throws ArgumentError Circuitscape.pair_jobs([10], [10, 5], points, orig_pts, exclude, false)
end

# Node currents are computed in a single sweep over the CSC structure; the
# previous implementation (six sparse temporaries per pos/neg pair) is kept
# here as the reference.
@testset "get_node_currents matches sparse reference" begin
    function ref_branch_posneg(G::SparseArrays.SparseMatrixCSC{T}, v, pos) where {T}
        b = T[]
        for i = 1:size(G, 1), j in SparseArrays.nzrange(G, i)
            row = G.rowval[j]; val = G.nzval[j]
            i > row && push!(b, pos ? abs(val) * (v[row] - v[i]) : abs(val) * (v[i] - v[row]))
        end
        maxcur = maximum(b)
        for i in eachindex(b)
            abs(b[i] / maxcur) < 1e-8 && (b[i] = 0)
        end
        b
    end
    function ref_branch(G::SparseArrays.SparseMatrixCSC{T,V}, v, pos) where {T,V}
        b = ref_branch_posneg(G, v, pos)
        N = size(G, 1)
        I = V[]; J = V[]
        for i = 1:N, j in SparseArrays.nzrange(G, i)
            row = G.rowval[j]
            i > row && (push!(I, row); push!(J, i))
        end
        SparseArrays.sparse(I, J, b, N, N)
    end
    function ref_posneg(G, v, finitegrounds, pos)
        bc = ref_branch(G, v, pos)
        bc = bc - bc'
        for j in eachindex(bc.nzval)
            bc.nzval[j] < 0 && (bc.nzval[j] = 0)
        end
        SparseArrays.dropzeros!(bc)
        if finitegrounds[1] != -9999
            fgc = finitegrounds .* v
            if pos
                map!(x -> x < 0 ? -x : 0, fgc, fgc)
            else
                map!(x -> x > 0 ? x : 0, fgc, fgc)
            end
            bc = bc + SparseArrays.spdiagm(0 => fgc)
        end
        vec(sum(bc, dims = 1))
    end
    ref_node_currents(G, v, fg) = map((x, y) -> x > y ? x : y,
                                      ref_posneg(G, v, fg, true), ref_posneg(G, v, fg, false))

    function check(G::SparseArrays.SparseMatrixCSC{T}, v, fg) where {T}
        new = Circuitscape.get_node_currents(G, v, fg)
        old = ref_node_currents(G, v, fg)
        @test new isa Vector{T}
        @test length(new) == length(old)
        tol = T == Float64 ? 1e-12 : 1e-6
        @test all(isapprox.(new, old, rtol = tol, atol = tol * maximum(abs, old)))
        # Something non-trivial was computed
        @test maximum(old) > 0
    end

    Random.seed!(11)
    for T in (Float64, Float32)
        # Small grid Laplacian, 8-neighbour connectivity
        G = model_problem(T, 7)
        n = size(G, 1)
        v = randn(T, n)
        check(G, v, T[-9999.])
        fg = T.(rand(n) .< 0.3) .* rand(T, n)
        check(G, v, fg)

        # Random sparse symmetric Laplacian, including isolated nodes
        n = 120
        B = SparseArrays.sprand(T, n, n, 0.05)
        A = B + B'
        A[SparseArrays.diagind(A)] .= 0
        SparseArrays.dropzeros!(A)
        L = SparseArrays.spdiagm(0 => vec(sum(A, dims = 1))) - A
        v = randn(T, n)
        check(L, v, T[-9999.])
        check(L, v, T.(rand(n) .< 0.5) .* rand(T, n))

        # Sentinel path with a Float64 sentinel on a Float32 problem, as core.jl passes
        check(G, randn(T, size(G, 1)), [-9999.])
    end

    # Hand-checked: a 3-node path 1 - 2 - 3 with unit conductances and
    # voltages 1, 0.5, 0: half an ampere flows in at 1, through 2, out at 3.
    L = SparseArrays.sparse([1, 1, 2, 2, 2, 3, 3], [1, 2, 1, 2, 3, 2, 3],
                            [1.0, -1, -1, 2, -1, -1, 1], 3, 3)
    @test Circuitscape.get_node_currents(L, [1.0, 0.5, 0.0], [-9999.]) ≈ [0.5, 0.5, 0.5]
    # A finite ground at node 3 carries 2*0 = 0; at node 1 it is a source
    # of 3*1 = 3 flowing in, which beats the 0.5 flowing out to node 2.
    @test Circuitscape.get_node_currents(L, [1.0, 0.5, 0.0], [-3.0, 0.0, 2.0]) ≈ [3.0, 0.5, 0.5]
end
