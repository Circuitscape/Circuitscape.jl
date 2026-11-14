import Circuitscape: construct_node_map, compute_omniscape_current
using Circuitscape
using Circuitscape: compute
import LinearAlgebra, Logging
import Graphs, SparseArrays, Random
import Logging

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
        cs_cfg["cholmod_batch_size"] = "32"
        cs_cfg["data_type"] = "raster"

        current = compute_omniscape_current(
                conductance,
                source,
                ground,
                cs_cfg
        )

        # Omniscape calls compute_omniscape_current from many tasks at once.
        # Concurrent calls must neither error (the default CSTIMER is shared
        # and not thread-safe) nor perturb each other's results.
        cs_cfg["solver"] = "cg+amg"
        serial = compute_omniscape_current(conductance, source, ground, cs_cfg)
        results = Vector{Matrix{Float64}}(undef, 64)
        @sync for i in eachindex(results)
                Threads.@spawn results[i] = compute_omniscape_current(
                        conductance, source, ground, cs_cfg)
        end
        @test all(r -> r ≈ serial, results)
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

    # As in Circuitscape 4's get_poly_map_temp: polymap numbers run to 9, and
    # points 1, 2, 3, 4, 7 get ids 10..14 in turn, each absorbing the
    # polygons it overlaps. Point 4 at (4,1) lies on polygon 1, which region
    # 1 already absorbed, so region 1 and point 4 become one node (13);
    # point 7 at (3,4) lies on polygon 9, absorbed by point 2, so they become
    # one node (14).
    @test r == [ 13.0  14.0  0.0   0.0  0.0
                  0.0   0.0  0.0   0.0  0.0
                 12.0   0.0  0.0  14.0  0.0
                 13.0   0.0  0.0   0.0  0.0
                 13.0   0.0  0.0   0.0  14.0 ]
end

@testset "create_new_polymap: focal regions over polygons" begin
    # Rows, columns and ids of the focal cells, as in RasterData.points_rc.
    prc(cells, ids) = (first.(cells), last.(cells), ids)
    poly = [1 1 0 0
            0 0 0 0
            0 0 2 2
            0 0 0 0]
    gmap = ones(4, 4)

    # Region 5 touches no polygon: its cells get a new id, polygons untouched.
    r = Circuitscape.create_new_polymap(gmap, poly,
            prc([(2,1), (2,2), (4,4)], [5, 5, 6]), 5, 6)
    @test r == [1 1 0 0
                3 3 0 0
                0 0 2 2
                0 0 0 0]

    # Region 5 overlaps polygon 1 in one cell: the polygon and every cell
    # of the region become one node (this branch used to throw). The return
    # type must be inferred (a captured, reassigned `newpoly` made it `Any`).
    r = @inferred Circuitscape.create_new_polymap(gmap, poly,
            prc([(1,2), (2,2), (4,4)], [5, 5, 6]), 5, 6)
    @test r == [3 3 0 0
                0 3 0 0
                0 0 2 2
                0 0 0 0]

    # Region 5 overlaps polygons 1 and 2: both polygons and all of the
    # region's cells, including (2,3) which lies on no polygon, merge.
    r = Circuitscape.create_new_polymap(gmap, poly,
            prc([(1,2), (2,2), (2,3), (3,3), (4,4)], [5, 5, 5, 5, 6]), 5, 6)
    @test r == [3 3 0 0
                0 3 3 0
                0 0 3 3
                0 0 0 0]

    # Both regions touch polygon 1: they collapse into one node.
    r = Circuitscape.create_new_polymap(gmap, poly,
            prc([(1,2), (2,2), (1,1), (4,1)], [5, 5, 6, 6]), 5, 6)
    @test r == [4 4 0 0
                0 4 0 0
                0 0 2 2
                4 0 0 0]
end

@testset "create_new_polymap: one-to-all point_map form" begin
    # Mirrors Circuitscape 4's get_poly_map_temp: for each point id in turn
    # the polygons it overlaps and all of its cells get one new id.
    prc(cells, ids) = (first.(cells), last.(cells), ids)
    function pmap(points_rc, n)
        pm = zeros(Int, n, n)
        for k in eachindex(points_rc[3])
            pm[points_rc[1][k], points_rc[2][k]] = points_rc[3][k]
        end
        pm
    end
    gmap = ones(4, 4)

    # Polygon 5 on (1,1),(1,2); region 1 on (1,2),(2,3); point 2 at (4,4).
    # Region 1 and polygon 5 become one node (id 6), point 2 gets id 7.
    poly = zeros(Int, 4, 4); poly[1,1] = 5; poly[1,2] = 5
    points_rc = prc([(1,2), (2,3), (4,4)], [1, 1, 2])
    r = Circuitscape.create_new_polymap(gmap, poly, points_rc, 0, 0, pmap(points_rc, 4))
    @test r == [6 6 0 0
                0 0 6 0
                0 0 0 0
                0 0 0 7]
    nm = Circuitscape.construct_node_map(gmap, r)
    @test nm[1,1] == nm[1,2] == nm[2,3]
    @test nm[4,4] != nm[1,1]

    # A polygon elsewhere whose number equals a point id must not be merged
    # with that point's region.
    poly = zeros(Int, 4, 4); poly[1,1] = 5; poly[1,2] = 5; poly[4,1] = 2; poly[4,2] = 2
    points_rc = prc([(1,2), (2,2), (3,4)], [2, 2, 1])
    r = Circuitscape.create_new_polymap(gmap, poly, points_rc, 0, 0, pmap(points_rc, 4))
    @test r == [7 7 0 0
                0 7 0 0
                0 0 0 6
                2 2 0 0]

    # Two regions sharing a polygon collapse into one node.
    poly = zeros(Int, 4, 4); poly[1,1] = 1; poly[1,2] = 1
    points_rc = prc([(1,2), (2,2), (1,1), (4,1)], [5, 5, 6, 6])
    r = Circuitscape.create_new_polymap(gmap, poly, points_rc, 0, 0, pmap(points_rc, 4))
    @test r == [3 3 0 0
                0 3 0 0
                0 0 0 0
                3 0 0 0]

    # Single-cell points only: a point on a polygon joins it, others get
    # their own ids (the point_file_no_polygons branch).
    poly = zeros(Int, 4, 4); poly[1,1] = 1; poly[1,2] = 1
    points_rc = prc([(1,2), (3,3)], [1, 2])
    r = Circuitscape.create_new_polymap(gmap, poly, points_rc, 0, 0, pmap(points_rc, 4))
    nm = Circuitscape.construct_node_map(gmap, r)
    @test nm[1,1] == nm[1,2]
    @test nm[3,3] != nm[1,1]
end

@testset "construct_graph: polygons over nodata cells" begin
    # A polygon covering a nodata cell must not connect that cell's habitat
    # neighbours to the polygon node (Circuitscape 4 only links two cells
    # when both have nonzero conductance).
    g = [0.0 1.0 1.0
         1.0 1.0 1.0
         1.0 1.0 1.0]
    poly = zeros(Int, 3, 3); poly[1,1] = 1; poly[3,3] = 1
    nm = Circuitscape.construct_node_map(g, poly)
    @test nm[1,1] == nm[3,3] != 0
    for avg_res in (false, true), four in (false, true)
        A = Circuitscape.construct_graph(g, nm, avg_res, four)
        pn = nm[1,1]
        @test A[pn, nm[1,2]] == 0
        @test A[pn, nm[2,1]] == 0
        @test A[pn, nm[2,2]] ≈ (four ? 0.0 : 1 / √2)
        @test A[pn, nm[2,3]] == 1
        @test A[pn, nm[3,2]] == 1
        @test SparseArrays.nnz(A) == 2 * Circuitscape.count_graph_edges(g, nm, four)
        @test LinearAlgebra.issymmetric(A)
    end
    # Isolated polygon entirely over nodata: a node with no edges
    g2 = ones(3, 3); g2[1,1] = 0
    poly2 = zeros(Int, 3, 3); poly2[1,1] = 1
    nm2 = Circuitscape.construct_node_map(g2, poly2)
    @test nm2[1,1] == 0
end

@testset "read_cellmap: negative values are null" begin
    # runtests.jl disables warnings globally; turn them back on here so that
    # @test_logs can see the warning, and restore runtests.jl's setting after.
    Logging.disable_logging(Logging.BelowMinLevel)
    try
    d = mktempdir()
    hab = joinpath(d, "hab.asc")
    open(hab, "w") do io
        println(io, "ncols 3\nnrows 3\nxllcorner 0\nyllcorner 0\ncellsize 1\nNODATA_value -9999")
        println(io, "-2 1 1\n1 -9999 1\n1 1 4")
    end
    for is_res in (true, false)
        gmap, _ = @test_logs (:warn, r"1 cell\(s\) .* negative") Circuitscape.read_cellmap(hab, is_res, Float64)
        @test gmap[1,1] == 0
        @test gmap[2,2] == 0
        @test gmap[3,3] == (is_res ? 0.25 : 4.0)
        @test all(gmap .>= 0)
    end
    # A polygon covering the negative cell does not reach across it
    gmap, _ = @test_logs (:warn,) Circuitscape.read_cellmap(hab, true, Float64)
    poly = zeros(Int, 3, 3); poly[1,1] = 1; poly[3,3] = 1
    nm = Circuitscape.construct_node_map(gmap, poly)
    A = Circuitscape.construct_graph(gmap, nm, false, true)
    @test A[nm[1,1], nm[1,2]] == 0
    @test A[nm[1,1], nm[2,1]] == 0
    # No negative cells: no warning
    open(hab, "w") do io
        println(io, "ncols 2\nnrows 1\nxllcorner 0\nyllcorner 0\ncellsize 1\nNODATA_value -9999")
        println(io, "-9999 2")
    end
    gmap, _ = @test_logs Circuitscape.read_cellmap(hab, true, Float64)
    @test gmap == [0.0 0.5]
    finally
        Logging.disable_logging(Logging.Warn)
    end
end

@testset "construct_graph: non-positive in-memory cells inside a polygon" begin
    # API callers (Omniscape) pass conductance arrays directly, with -9999
    # for nodata. A polygon over such a cell must not link it either.
    for bad in (-9999.0, -1.0, NaN)
        g = [bad 1.0; 1.0 1.0]
        poly = [1 0; 0 1]
        nm = Circuitscape.construct_node_map(g, poly)
        @test nm[1,1] == nm[2,2] != 0
        for avg_res in (false, true), four in (false, true)
            A = Circuitscape.construct_graph(g, nm, avg_res, four)
            @test all(isfinite, SparseArrays.nonzeros(A))
            @test all(>=(0), SparseArrays.nonzeros(A))
            @test A[nm[1,1], nm[1,2]] == A[nm[1,1], nm[2,1]] == 1
            @test SparseArrays.nnz(A) == 2 * Circuitscape.count_graph_edges(g, nm, four)
        end
    end
end

@testset "create_new_polymap: focal regions over polygons, randomised" begin
    # Reference: Circuitscape 4's get_poly_map_temp / get_overlap_polymap,
    # relabelling the running map once per point id.
    function cs4_poly_map_temp(poly_map, point_map)
        m = copy(poly_map)
        k = maximum(poly_map)
        for p in sort!(unique(point_map[point_map .!= 0]))
            k += 1
            cells = findall(==(p), point_map)
            vals = unique(m[c] for c in cells if m[c] != 0)
            for i in eachindex(m)
                m[i] in vals && (m[i] = k)
            end
            m[cells] .= k
        end
        m
    end
    rng = Random.MersenneTwister(523)
    gmap = ones(12, 12)
    for trial in 1:300
        poly = rand(rng, 0:6, 12, 12) .* (rand(rng, 12, 12) .< 0.3)
        pm = rand(rng, 0:8, 12, 12) .* (rand(rng, 12, 12) .< 0.25)
        I = findall(!iszero, pm)
        points_rc = (getindex.(I, 1), getindex.(I, 2), pm[I])
        length(unique(points_rc[3])) == length(points_rc[3]) && continue
        r = @inferred Circuitscape.create_new_polymap(gmap, poly, points_rc, 0, 0, pm)
        @test r == cs4_poly_map_temp(poly, pm)
    end
end

import Circuitscape: resolve_conflicts

# The finite grounds are those left after the conflicts are resolved (a
# removed ground does not stay on the solver diagonal), and rmvall removes
# both sides; -9999 marks "no finite grounds".
@test resolve_conflicts([1.,0.,0.], [1.,0.,0.], :rmvgnd) == ([1, 0, 0], [0, 0, 0], [-9999])
@test resolve_conflicts([1.,0.,0.], [1.,0.,0.], :rmvsrc) == ([0, 0, 0], [1, 0, 0], [1, 0, 0])
@test resolve_conflicts([1.,0.,0.], [1.,0.,0.], :keepall) == ([1, 0, 0], [1, 0, 0], [1, 0, 0])
@test resolve_conflicts([1.,0.,0.], [1.,0.,0.], :rmvall) == ([0, 0, 0], [0, 0, 0], [-9999])
@test resolve_conflicts([1.,1.,0.], [1.,0.,2.], :rmvgnd) == ([1, 1, 0], [0, 0, 2], [0, 0, 2])
@test resolve_conflicts([1.,1.,0.], [1.,0.,2.], :rmvall) == ([0, 1, 0], [0, 0, 2], [0, 0, 2])
# A direct ground under a source is removed whatever the source's sign
@test resolve_conflicts([-1.,1.,0.], [Inf,0.,2.], :keepall) == ([-1, 1, 0], [0, 0, 2], [0, 0, 2])

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

        # laplacian! inserts the diagonal, or fills it in place when stored;
        # the isolated node 5 gets no entry either way.
        L = Circuitscape.laplacian!(copy(A))
        @test L ≈ SparseArrays.spdiagm(0 => vec(sum(A, dims = 1))) - A
        @test SparseArrays.nnz(L) == SparseArrays.nnz(A) + 4
        B = A + LinearAlgebra.I
        @test Circuitscape.laplacian!(B) === B
        @test B == L && B.rowval == L.rowval
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
    @test all(Circuitscape._is_unset(v) for (k, v) in Circuitscape.INI_PLACEHOLDERS
              if Symbol(k) in fieldnames(Circuitscape.CSConfig) &&
                 fieldtype(Circuitscape.CSConfig, Symbol(k)) === String)
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

# INI files are read as Circuitscape 4 (Python's ConfigParser) read them.
# Before, a comment or blank line became a key "", `precision: single` was
# dropped (running double) and `Solver = cholmod` was ignored (running cg+amg).
@testset "INI syntax" begin
    ini = joinpath(mktempdir(), "job.ini")
    write(ini, "[Circuitscape Mode]\n# a comment\n; another\n\n" *
               "Scenario = pairwise\nSOLVER = cholmod\nprecision: single\n" *
               "habitat_file = C:\\data\\cell map.asc\r\n")
    Logging.disable_logging(Logging.Debug)
    cfg = try
        # no warning: comments no longer become a key ""
        @test_logs Circuitscape.parse_config(ini)
    finally
        Logging.disable_logging(Logging.Warn)
    end
    @test cfg.scenario == Circuitscape.sc_pairwise
    @test cfg.solver == Circuitscape.st_cholmod
    @test cfg.precision == Circuitscape.pr_single
    @test cfg.habitat_file == "C:\\data\\cell map.asc"

    # A line that is not `key = value`, a comment or a section is an error
    write(ini, "scenario = pairwise\nprecision single\n")
    @test_throws "job.ini, line 2: expected `key = value`" Circuitscape.parse_config(ini)

    # write_config lays out every option once, in Circuitscape 4's sections
    keys = reduce(vcat, last.(Circuitscape.INI_SECTIONS))
    @test allunique(keys) && sort(keys) == sort(collect(fieldnames(Circuitscape.CSConfig)))
end

# The options reference restates each key's default and accepted values;
# keep it in step with the code.
@testset "options reference matches the code" begin
    C = Circuitscape
    doc = read(joinpath(@__DIR__, "..", "docs", "src", "options.md"), String)
    rows = Dict(m[1] => (default = strip(m[2], ['`', ' ']), text = m[3])
                for m in eachmatch(r"^\| `(\w+)` \| \w+ \| ([^|]*) \|([^\n]*)$"m, doc))
    fields = String.(fieldnames(C.CSConfig))
    @test sort!(collect(keys(rows))) == sort!(setdiff(fields, ["version"]))
    # defaults; "—" is a path with no default
    @test isempty([k for (k, r) in rows if !(C._same(r.default, C.INI_DEFAULTS[k]) ||
                                              r.default == "—" && C._is_unset(C.INI_DEFAULTS[k]))])
    # every accepted spelling of every enumerated option
    missing_names = String[]
    for (k, r) in rows
        T = fieldtype(C.CSConfig, Symbol(k))
        (T <: C.Named && T !== Bool) || continue
        for (_, names) in C._names(T), n in names
            occursin(lowercase("`$n`"), lowercase(r.text)) || push!(missing_names, "$k: $n")
        end
    end
    @test isempty(missing_names)
end

# Keys and option values in any case, in a dict as in an INI file (#534);
# paths keep their case, and typos are still errors.
@testset "Case-insensitive configuration" begin
    C = Circuitscape
    d = C.init_config()
    C.update!(d, Dict("Scenario" => "pairwise", "solver" => "CHOLMOD",
                      "write_cur_maps" => "TRUE", "habitat_file" => "Data/Cells.ASC"))
    @test !haskey(d, "Scenario")
    cfg = C.CSConfig(d)
    @test cfg.scenario == C.sc_pairwise && cfg.solver == C.st_cholmod && cfg.write_cur_maps
    @test cfg.habitat_file == "Data/Cells.ASC"

    # CSConfig(dict) folds key case itself, and refuses a conflict
    @test C.CSConfig(Dict("scenario" => "pairwise", "Solver" => "cholmod")).solver == C.st_cholmod
    @test_throws "given twice, in different case" C.CSConfig(
        Dict("scenario" => "pairwise", "Solver" => "cholmod", "solver" => "cg+amg"))

    @test_throws "solver = \"CHOLMODD\" is not recognised" C.CSConfig(
        Dict("scenario" => "pairwise", "solver" => "CHOLMODD"))
    @test_throws "write_cur_maps = \"TURE\" is not recognised" C.CSConfig(
        Dict("scenario" => "pairwise", "write_cur_maps" => "TURE"))
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

# Advanced mode resolves source/ground conflicts the way Python does: on the
# raster cell maps as they are read (rmvall drops both sides), then on the
# nodes; a component is solved as long as it has a nonzero source and ground.
@testset "advanced conflict policies" begin
    dir = mktempdir()
    hdr = "ncols 3\nnrows 3\nxllcorner 0\nyllcorner 0\ncellsize 1\nNODATA_value -9999\n"
    asc(rows) = hdr * join((join(r, " ") for r in rows), "\n") * "\n"
    write(joinpath(dir, "cell.asc"), asc(fill([1, 1, 1], 3)))
    # (2,2) is both a source and a finite ground; (1,1) a plain source and
    # (3,3) a plain ground.
    write(joinpath(dir, "src.asc"), asc([[1, 0, 0], [0, 1, 0], [0, 0, 0]]))
    write(joinpath(dir, "gnd.asc"), asc([[0, 0, 0], [0, 0.5, 0], [0, 0, 2]]))
    function config(; kw...)
        d = Circuitscape.init_config()
        d["data_type"] = "raster"; d["scenario"] = "advanced"
        d["habitat_file"] = joinpath(dir, "cell.asc")
        d["source_file"] = joinpath(dir, "src.asc")
        d["ground_file"] = joinpath(dir, "gnd.asc")
        d["ground_file_is_resistances"] = "False"
        d["connect_four_neighbors_only"] = "True"
        d["output_file"] = joinpath(dir, "out.out")
        for (k, v) in kw; d[string(k)] = v; end
        d
    end
    function vectors(policy)
        prob = Circuitscape.build_problem(Float64, Circuitscape.CSConfig(config(remove_src_or_gnd = policy)))
        nm = prob.geometry.nodemap
        at(v) = Dict(cell => v[nm[cell...]] for cell in ((1,1), (2,2), (3,3)))
        at(prob.sources), at(prob.grounds), at(prob.finitegrounds)
    end
    s, g, f = vectors("keepall")
    @test s == Dict((1,1) => 1, (2,2) => 1, (3,3) => 0)
    @test g == Dict((1,1) => 0, (2,2) => 0.5, (3,3) => 2)
    @test f == g
    s, g, f = vectors("rmvsrc")
    @test s == Dict((1,1) => 1, (2,2) => 0, (3,3) => 0)
    @test g == Dict((1,1) => 0, (2,2) => 0.5, (3,3) => 2)
    @test f == g
    # rmvgnd removes the finite ground everywhere, the solver diagonal included
    s, g, f = vectors("rmvgnd")
    @test s == Dict((1,1) => 1, (2,2) => 1, (3,3) => 0)
    @test g == Dict((1,1) => 0, (2,2) => 0, (3,3) => 2)
    @test f == g
    # rmvall removes both sides of the conflict
    s, g, f = vectors("rmvall")
    @test s == Dict((1,1) => 1, (2,2) => 0, (3,3) => 0)
    @test g == Dict((1,1) => 0, (2,2) => 0, (3,3) => 2)
    @test f == g

    # rmvgnd is the same run as one without the conflicting ground
    v_rmvgnd = compute(config(remove_src_or_gnd = "rmvgnd"))
    write(joinpath(dir, "gnd_only33.asc"), asc([[0, 0, 0], [0, 0, 0], [0, 0, 2]]))
    v_ref = compute(config(ground_file = joinpath(dir, "gnd_only33.asc")))
    @test v_rmvgnd ≈ v_ref

    # A map emptied by the read-time removal is an error, as in Python
    write(joinpath(dir, "src22.asc"), asc([[0, 0, 0], [0, 1, 0], [0, 0, 0]]))
    write(joinpath(dir, "gnd22.asc"), asc([[0, 0, 0], [0, 0.5, 0], [0, 0, 0]]))
    @test_throws ErrorException("No valid sources detected. Please check source file") compute(
        config(source_file = joinpath(dir, "src22.asc"), remove_src_or_gnd = "rmvsrc"))
    @test_throws ErrorException("No valid grounds detected. Please check ground file") compute(
        config(ground_file = joinpath(dir, "gnd22.asc"), remove_src_or_gnd = "rmvgnd"))
    @test_throws ErrorException("No valid sources detected. Please check source file") compute(
        config(source_file = joinpath(dir, "src22.asc"), remove_src_or_gnd = "rmvall"))
    # keepall keeps both, and a source on a direct ground (the only ground of
    # its component) leaves the node-level pass with nothing to solve
    @test compute(config(source_file = joinpath(dir, "src22.asc"), remove_src_or_gnd = "keepall")) isa Matrix
    @test_throws ErrorException compute(
        config(source_file = joinpath(dir, "src22.asc"), ground_file = joinpath(dir, "gnd22.asc"),
               use_direct_grounds = "True"))

    # Sources of +1 and -1 sum to zero but are still sources
    write(joinpath(dir, "srcpm.asc"), asc([[1, 0, -1], [0, 0, 0], [0, 0, 0]]))
    v = compute(config(source_file = joinpath(dir, "srcpm.asc"), ground_file = joinpath(dir, "gnd_only33.asc")))
    @test size(v) == (3, 3) && v[1, 1] > 0 && v[1, 3] < 0

    # Network: building the problem twice from one NetworkData must not
    # invert the ground resistances twice
    write(joinpath(dir, "g.txt"), "1\t2\t1.0\n2\t3\t1.0\n")
    write(joinpath(dir, "nsrc.txt"), "1\t1.0\n")
    write(joinpath(dir, "ngnd.txt"), "3\t4.0\n")
    d = Circuitscape.init_config()
    d["data_type"] = "network"; d["scenario"] = "advanced"
    d["habitat_file"] = joinpath(dir, "g.txt"); d["habitat_map_is_resistances"] = "True"
    d["source_file"] = joinpath(dir, "nsrc.txt"); d["ground_file"] = joinpath(dir, "ngnd.txt")
    d["ground_file_is_resistances"] = "True"; d["output_file"] = joinpath(dir, "n.out")
    cfg = Circuitscape.CSConfig(d)
    data = Circuitscape.load_data(Float64, cfg)
    p1 = Circuitscape.build_problem(data, cfg)
    p2 = Circuitscape.build_problem(data, cfg)
    @test data.ground_map[1, 2] == 4.0
    @test p1.grounds == p2.grounds && maximum(p1.grounds) ≈ 0.25

    # Network conflicts are resolved on the nodes only. rmvgnd takes the
    # ground off the solver diagonal too, and rmvall removes both sides (Python
    # keeps the finite ground under rmvgnd and the ground under rmvall here).
    write(joinpath(dir, "nsrc12.txt"), "1\t1.0\n2\t1.0\n")
    write(joinpath(dir, "ngnd13.txt"), "1\t4.0\n3\t1.0\n")
    write(joinpath(dir, "nsrc2.txt"), "2\t1.0\n")
    write(joinpath(dir, "ngnd3.txt"), "3\t1.0\n")
    function netrun(src, gnd, pol)
        dd = copy(d)
        dd["source_file"] = joinpath(dir, src); dd["ground_file"] = joinpath(dir, gnd)
        dd["remove_src_or_gnd"] = pol
        compute(dd)
    end
    @test netrun("nsrc12.txt", "ngnd13.txt", "rmvgnd") ≈ netrun("nsrc12.txt", "ngnd3.txt", "keepall")
    @test netrun("nsrc12.txt", "ngnd13.txt", "rmvall") ≈ netrun("nsrc2.txt", "ngnd3.txt", "keepall")
    @test netrun("nsrc12.txt", "ngnd13.txt", "rmvsrc") ≈ netrun("nsrc2.txt", "ngnd13.txt", "keepall")
    # A component left with no source is an error in advanced mode
    @test_throws ErrorException("All sources conflicted with grounds and were removed. There is nothing to solve.") netrun(
        "nsrc.txt", "ngnd13.txt", "rmvsrc")

    # One-to-all / all-to-one: a focal node short-circuited to every other
    # focal node (here by a polygon) has nothing to solve and reports -1, as
    # a lone focal node does, instead of ending the run.
    write(joinpath(dir, "pts.asc"), asc([[1, 2, -9999], [-9999, -9999, -9999], [-9999, -9999, -9999]]))
    write(joinpath(dir, "poly.asc"), asc([[5, 5, -9999], [-9999, -9999, -9999], [-9999, -9999, -9999]]))
    for sc in ("one-to-all", "all-to-one")
        dd = config(scenario = sc, point_file = joinpath(dir, "pts.asc"),
                    use_polygons = "True", polygon_file = joinpath(dir, "poly.asc"))
        @test compute(dd) == [1 -1; 2 -1]
    end
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
    # One pass: chained rows must not compound (1 -> 2 -> 3), whatever their
    # order in the file. Circuitscape 4 compounded them; that is not kept.
    m = [1.0 2.0; 2.0 3.0]
    Circuitscape.reclassify!(m, [1.0 2.0; 2.0 3.0])
    @test m == [2.0 3.0; 3.0 3.0]
    m = [1.0 2.0; 2.0 3.0]
    Circuitscape.reclassify!(m, [2.0 3.0; 1.0 2.0])
    @test m == [2.0 3.0; 3.0 3.0]
    # Swapping two values works as a lookup
    m = [1.0 2.0; 2.0 3.0]
    Circuitscape.reclassify!(m, [1.0 2.0; 2.0 1.0])
    @test m == [2.0 1.0; 1.0 3.0]
    # Unrelated rows are a plain lookup
    m = [1.0 2.0; 2.0 3.0]
    Circuitscape.reclassify!(m, [1.0 5.0; 3.0 7.0])
    @test m == [5.0 2.0; 2.0 7.0]

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

# With no map output requested, postprocessing is skipped entirely: no
# per-pair files in either mode, and resistances identical to a run with maps.
@testset "postprocess skipped without map outputs" begin
    dir = mktempdir()
    hdr = "ncols 5\nnrows 5\nxllcorner 0\nyllcorner 0\ncellsize 1\nNODATA_value -9999\n"
    write(joinpath(dir, "cell.asc"), hdr * join(fill(join(fill("1", 5), " "), 5), "\n") * "\n")
    write(joinpath(dir, "pts.asc"), hdr * "1 0 0 0 2\n0 0 0 0 0\n0 0 0 0 0\n0 0 0 0 0\n3 0 0 0 0\n")
    write(joinpath(dir, "excl.txt"), "mode\texclude\n1\t2\n")   # defeat the resistance shortcut
    d = Circuitscape.init_config()
    d["data_type"] = "raster"; d["scenario"] = "pairwise"
    d["habitat_file"] = joinpath(dir, "cell.asc"); d["point_file"] = joinpath(dir, "pts.asc")
    d["use_included_pairs"] = "True"; d["included_pairs_file"] = joinpath(dir, "excl.txt")
    d["output_file"] = joinpath(dir, "maps.out"); d["write_cur_maps"] = "True"
    r_maps = compute(d)
    @test any(f -> occursin("_curmap_", f), readdir(dir))
    d["output_file"] = joinpath(dir, "nomaps.out"); d["write_cur_maps"] = "False"
    r_none = compute(d)
    @test r_none == r_maps
    @test !any(f -> startswith(f, "nomaps") && (occursin("curmap", f) || occursin("voltmap", f)), readdir(dir))

    # Network mode used to write per-pair current files regardless of the flag
    write(joinpath(dir, "g.txt"), "1\t2\t1.0\n2\t3\t1.0\n3\t4\t1.0\n")
    write(joinpath(dir, "fp.txt"), "1\n4\n")
    n = Circuitscape.init_config()
    n["data_type"] = "network"; n["scenario"] = "pairwise"
    n["habitat_file"] = joinpath(dir, "g.txt"); n["point_file"] = joinpath(dir, "fp.txt")
    n["output_file"] = joinpath(dir, "net.out"); n["write_cur_maps"] = "False"
    rn = compute(n)
    @test rn[2, 3] ≈ 3.0
    @test !any(f -> startswith(f, "net") && occursin("currents", f), readdir(dir))
end

# 32-bit indices: chosen automatically when the problem fits, forced to
# 64-bit by use_64bit_indexing, and kept through the whole pipeline. The
# regularization used to promote a 32-bit matrix (and so the CHOLMOD factor)
# back to 64-bit indices via `sparse(c * I, n, n)`.
@testset "index type" begin
    SparseMatrixCSC, sparse, nnz = SparseArrays.SparseMatrixCSC, SparseArrays.sparse, SparseArrays.nnz
    I = LinearAlgebra.I
    cfg =Circuitscape.parse_config("input/raster/one_to_all/11/oneToAllVerify11.ini")
    cfg64 = Circuitscape.CSConfig(cfg; use_64bit_indexing = true)
    @test Circuitscape.index_type(cfg, 10, 10) == Int32
    @test Circuitscape.index_type(cfg, Circuitscape.INT32_MAX_NODES_ITERATIVE - 1, typemax(Int32) - 1) == Int32
    @test Circuitscape.index_type(cfg, Circuitscape.INT32_MAX_NODES_ITERATIVE, 10) == Int64
    cfgch = Circuitscape.CSConfig(cfg; solver = Circuitscape.st_cholmod)
    @test Circuitscape.index_type(cfgch, Circuitscape.INT32_MAX_NODES_DIRECT - 1, 10) == Int32
    @test Circuitscape.index_type(cfgch, Circuitscape.INT32_MAX_NODES_DIRECT, 10) == Int64
    @test Circuitscape.index_type(cfgch, Circuitscape.INT32_MAX_NODES_ITERATIVE - 1, 10) == Int64
    @test Circuitscape.index_type(cfg, 10, typemax(Int32)) == Int64
    @test Circuitscape.index_type(cfg64, 10, 10) == Int64
    @test Circuitscape.raster_nnz_bound(1_000_000) > nnz(model_problem(1000)) 

    @test Circuitscape.load_raster_data(Float64, cfg) isa Circuitscape.RasterData{Float64,Int32}
    @test Circuitscape.load_raster_data(Float64, cfg64) isa Circuitscape.RasterData{Float64,Int64}
    @test Circuitscape.load_raster_data(Float32, Int64, cfg) isa Circuitscape.RasterData{Float32,Int64}
    net = Circuitscape.parse_config("input/network/sgNetworkVerify1.ini")
    @test Circuitscape.get_network_data(Float64, net) isa Circuitscape.NetworkData{Float64,Int32}
    @test Circuitscape.get_network_data(Float64, Int64, net) isa Circuitscape.NetworkData{Float64,Int64}
    @test Circuitscape.get_network_data(Float64, net).coords ==
          Circuitscape.get_network_data(Float64, Int64, net).coords

    for V in (Int32, Int64)
        prob = Circuitscape.build_problem(Circuitscape.load_data(Float64, V, cfg), cfg)
        @test prob.G isa SparseMatrixCSC{Float64,V}
        @test prob.cc isa Vector{Vector{V}}
        n = size(prob.G, 1)
        R = Circuitscape.regularize(prob.G)
        @test R isa SparseMatrixCSC{Float64,V}
        @test R == prob.G + sparse(10eps() * I, n, n)
        F = Circuitscape.construct_cholesky_factor(prob.G, Circuitscape.CholmodSolver(10))
        @test F isa SparseArrays.CHOLMOD.Factor{Float64,V}
        # Advanced mode's finite grounds go on the diagonal with the same index type
        a = deepcopy(prob.G)
        sources = zeros(n); sources[1] = 1
        grounds = zeros(n); grounds[end] = 2.0
        finitegrounds = copy(grounds)
        v32 = Circuitscape.multiple_solver(cfg, Circuitscape.CholmodSolver(10), a, sources, grounds, finitegrounds)
        @test length(v32) == n
        # With a direct ground too, the finite one goes on the reduced copy
        # in place; `a` itself is never written to.
        grounds[2] = Inf
        v = Circuitscape.multiple_solver(cfg, Circuitscape.CholmodSolver(10), a, copy(sources), grounds, finitegrounds)
        @test a == prob.G
        keep = [1; 3:n]
        @test v[2] == 0
        @test v[keep] ≈ (Matrix(a) + LinearAlgebra.Diagonal(finitegrounds))[keep, keep] \ sources[keep]
    end
    # Same answers whichever index type is used
    r32 = Circuitscape.run_onetoall(Circuitscape.load_data(Float64, Int32, cfg), cfg)
    r64 = Circuitscape.run_onetoall(Circuitscape.load_data(Float64, Int64, cfg), cfg)
    @test r32 == r64
end

# Variable source strengths are matched to focal points by id, as in Python
# (`get_strengths_rc`): the order of the rows in the file does not matter and
# an id the file leaves out gets strength 1. Focal ids are taken as written,
# never shifted as network node ids are.
@testset "Variable source strengths by id" begin
    fixture = joinpath(@__DIR__, "input", "raster", "one_to_all", "12")
    dir = mktempdir()
    function job(scenario, strengths_file; write_cur_maps = false)
        d = Dict{String,String}(parse_config(joinpath(fixture, "oneToAllVerify12.ini")))
        for k in ("habitat_file", "point_file", "polygon_file", "mask_file", "included_pairs_file")
            d[k] = joinpath(@__DIR__, d[k])
        end
        d["scenario"] = scenario
        d["variable_source_file"] = strengths_file
        d["output_file"] = joinpath(dir, scenario, "out.out")
        d["write_cur_maps"] = write_cur_maps ? "True" : "False"
        d["solver"] = "cholmod"
        d["parallelize"] = "false"
        mkpath(joinpath(dir, scenario))
        compute(d)
    end
    ref = readdlm(joinpath(@__DIR__, "output_verify", "oneToAllVerify12_resistances.out"))
    rows = readdlm(joinpath(fixture, "variable_source_list.txt"))

    # Same rows in another order: identical result (id 3 has strength 0 and
    # must stay -1; by row position that zero would land on id 1)
    shuffled = joinpath(dir, "shuffled.txt")
    writedlm(shuffled, rows[[3, 1, 2, 5, 4, 6, 7, 9, 8, 10, 11], :])
    r = job("one-to-all", shuffled)
    @test check_resistances(ref, r, 1e-6, label = "shuffled strengths")

    # An id missing from the file has strength 1: one-to-all resistances
    # (voltage over strength) are unchanged rather than an error
    missing5 = joinpath(dir, "missing5.txt")
    writedlm(missing5, rows[rows[:, 1] .!= 5, :])
    @test check_resistances(ref, job("one-to-all", missing5), 1e-6, label = "id 5 missing")

    # All-to-one: the current maps depend on which strength each source gets.
    # (allToOneVerify12's maps are CS5's own, not Circuitscape 4's: see the
    # "All-to-one excluded points carry no source" testset.)
    cum_map() = readdlm(joinpath(dir, "all-to-one", "out_cum_curmap.asc"), skipstart = 6)
    job("all-to-one", shuffled; write_cur_maps = true)
    cum_ref = readdlm(joinpath(@__DIR__, "output_verify", "allToOneVerify12_cum_curmap.asc"), skipstart = 6)
    @test sum(abs2, cum_map() - cum_ref) < 1e-6
    # ... and a missing id sources 1, the same as listing it with strength 1
    job("all-to-one", missing5; write_cur_maps = true)
    cum_missing = cum_map()
    @test sum(abs2, cum_missing - cum_ref) > 1e-3
    unit5 = joinpath(dir, "unit5.txt")
    writedlm(unit5, vcat(rows[rows[:, 1] .!= 5, :], [5 1.0]))
    job("all-to-one", unit5; write_cur_maps = true)
    @test sum(abs2, cum_map() - cum_missing) < 1e-10

    # A focal id of 0 in the file must not shift the other ids (only network
    # node lists are 0-based)
    with0 = joinpath(dir, "with0.txt")
    writedlm(with0, vcat([0 7.0], rows))
    @test Circuitscape.read_variable_strengths(Float64, with0)[:, 1] == vcat(0, rows[:, 1])
    @test check_resistances(ref, job("one-to-all", with0), 1e-6, label = "id 0 row")
    write(joinpath(dir, "bad.txt"), "1 2 3\n")
    @test_throws ArgumentError Circuitscape.read_variable_strengths(Float64, joinpath(dir, "bad.txt"))
end

# All-to-one with an include file: a focal point the file excludes for the
# current ground is not a source. With unit strengths this is what Python
# does (it zeroes the point out of `unique_point_map`). With variable
# strengths Python 4.0.5 still injects 1 A at every excluded point (the
# zeroed id falls back to strength 1 in `get_strengths_rc`); CS5 does not,
# so the allToOneVerify12 reference maps were regenerated.
@testset "All-to-one excluded points carry no source" begin
    n = 7
    hdr = "ncols $n\nnrows $n\nxllcorner 0\nyllcorner 0\ncellsize 1\nNODATA_value -9999\n"
    grid(f) = join([join([string(f(i, j)) for j in 1:n], " ") for i in 1:n], "\n") * "\n"
    function run(dir, pts, include; strengths = nothing)
        write(joinpath(dir, "cell.asc"), hdr * grid((i, j) -> 1))
        write(joinpath(dir, "pts.asc"), hdr * grid((i, j) -> get(pts, (i, j), 0)))
        d = Circuitscape.init_config()
        d["data_type"] = "raster"; d["scenario"] = "all-to-one"
        d["habitat_file"] = joinpath(dir, "cell.asc")
        d["habitat_map_is_resistances"] = "True"
        d["point_file"] = joinpath(dir, "pts.asc")
        d["connect_four_neighbors_only"] = "True"
        d["output_file"] = joinpath(dir, "out.out")
        d["write_cur_maps"] = "True"
        d["solver"] = "cholmod"
        if include !== nothing
            write(joinpath(dir, "inc.txt"), include)
            d["use_included_pairs"] = "True"
            d["included_pairs_file"] = joinpath(dir, "inc.txt")
        end
        if strengths !== nothing
            write(joinpath(dir, "str.txt"), strengths)
            d["use_variable_source_strengths"] = "True"
            d["variable_source_file"] = joinpath(dir, "str.txt")
        end
        compute(d)
        readdlm(joinpath(dir, "out_curmap_1.asc"), skipstart = 6)
    end
    pts = Dict((1,1) => 1, (7,7) => 2, (1,7) => 3)
    # Pairs (1,2) and (2,3) included: with 1 as the ground, 3 is excluded
    c_inc = run(mktempdir(), pts, "mode\tinclude\n1\t2\n2\t3\n")
    # Reference: point 3 absent from the point file altogether
    c_ref = run(mktempdir(), filter(kv -> kv[2] != 3, pts), nothing)
    @test c_inc ≈ c_ref atol = 1e-8
    # Whereas with pair (1,3) included as well, point 3 does source current
    c_all = run(mktempdir(), pts, "mode\tinclude\n1\t2\n2\t3\n1\t3\n")
    @test c_all[1, 7] > c_ref[1, 7] + 0.5
    @test c_all[1, 1] ≈ 2 atol = 1e-8       # two unit sources into ground 1
    @test c_ref[1, 1] ≈ 1 atol = 1e-8

    # Variable strengths: the excluded point 3 injects nothing either (not
    # its own 4 A, and not the 1 A Python 4.0.5 gives it)
    str = "1\t3.0\n2\t2.5\n3\t4.0\n"
    v_inc = run(mktempdir(), pts, "mode\tinclude\n1\t2\n2\t3\n"; strengths = str)
    v_ref = run(mktempdir(), filter(kv -> kv[2] != 3, pts), nothing; strengths = str)
    @test v_inc ≈ v_ref atol = 1e-8
    @test v_inc[1, 1] ≈ 2.5 atol = 1e-8     # only point 2's 2.5 A reaches ground 1
    v_all = run(mktempdir(), pts, "mode\tinclude\n1\t2\n2\t3\n1\t3\n"; strengths = str)
    @test v_all[1, 1] ≈ 6.5 atol = 1e-8

    # The same on the allToOneVerify12 fixture: with ground 1 the include
    # file excludes points 4 and 9, so its current map is that of a run
    # whose point file lacks them (and the ids the include file never names)
    fixture = joinpath(@__DIR__, "input", "raster", "all_to_one", "12")
    d = Dict{String,String}(parse_config(joinpath(fixture, "allToOneVerify12.ini")))
    for k in ("habitat_file", "point_file", "polygon_file", "mask_file", "variable_source_file")
        d[k] = joinpath(@__DIR__, d[k])
    end
    dir = mktempdir()
    hdr12 = readlines(d["point_file"])[1:6]
    pts12 = readdlm(d["point_file"], skipstart = 6)
    open(joinpath(dir, "pts.asc"), "w") do io
        foreach(l -> println(io, l), hdr12)
        writedlm(io, map(x -> x in (1, 2, 3, 5, 6, 8) ? Int(x) : -9999, pts12), ' ')
    end
    d["point_file"] = joinpath(dir, "pts.asc")
    d["use_included_pairs"] = "False"
    d["output_file"] = joinpath(dir, "out.out")
    d["solver"] = "cholmod"; d["parallelize"] = "false"
    compute(d)
    @test sum(abs2, readdlm(joinpath(dir, "out_curmap_1.asc"), skipstart = 6) -
              readdlm(joinpath(@__DIR__, "output_verify", "allToOneVerify12_curmap_1.asc"), skipstart = 6)) < 1e-10

    # A multi-cell focal region left on its own by the include file: nothing
    # to solve, reported as -1 rather than attempted
    pts_region = Dict((1,1) => 1, (1,2) => 1, (7,7) => 2)
    dir = mktempdir()
    write(joinpath(dir, "cell.asc"), hdr * grid((i, j) -> 1))
    write(joinpath(dir, "pts.asc"), hdr * grid((i, j) -> get(pts_region, (i, j), 0)))
    write(joinpath(dir, "inc.txt"), "mode\texclude\n1\t2\n")
    d = Circuitscape.init_config()
    d["data_type"] = "raster"; d["scenario"] = "one-to-all"
    d["habitat_file"] = joinpath(dir, "cell.asc")
    d["habitat_map_is_resistances"] = "True"
    d["point_file"] = joinpath(dir, "pts.asc")
    d["connect_four_neighbors_only"] = "True"
    d["output_file"] = joinpath(dir, "out.out")
    d["use_included_pairs"] = "True"
    d["included_pairs_file"] = joinpath(dir, "inc.txt")
    d["solver"] = "cholmod"
    @test compute(d) == [1 -1; 2 -1]
end

# Network input read the way Circuitscape 4 reads it (compute.py read_graph /
# read_focal_nodes, csio.py load_graph).
@testset "Network I/O matches Circuitscape 4" begin
    dir = mktempdir()
    function network(name, graph, focal; extra...)
        d = Circuitscape.init_config()
        d["data_type"] = "network"; d["scenario"] = "pairwise"
        d["habitat_file"] = joinpath(dir, graph); d["habitat_map_is_resistances"] = "True"
        d["point_file"] = joinpath(dir, focal); d["output_file"] = joinpath(dir, name * ".out")
        for (k, v) in extra; d[string(k)] = v; end
        compute(d)
    end

    # Focal ids name graph nodes: on a 0-based graph they shift with the graph
    # (0-1-2-3 with 1, 2, 3 ohm; nodes 1 and 3 are 2 + 3 apart), a 1-based
    # graph never shifts them
    write(joinpath(dir, "g0.txt"), "0 1 1\n1 2 2\n2 3 3\n")
    write(joinpath(dir, "fp13.txt"), "1\n3\n")
    r = network("shift", "g0.txt", "fp13.txt")
    @test r[1, 2:end] == [2, 4] && r[2, 3] ≈ 5.0
    write(joinpath(dir, "fp03.txt"), "0\n3\n")
    r = network("shift0", "g0.txt", "fp03.txt")
    @test r[1, 2:end] == [1, 4] && r[2, 3] ≈ 6.0
    write(joinpath(dir, "g1.txt"), "1 2 1\n2 3 2\n3 4 3\n")
    r = network("noshift", "g1.txt", "fp13.txt")
    @test r[1, 2:end] == [1, 3] && r[2, 3] ≈ 3.0
    @test Circuitscape.read_focal_points(Int, joinpath(dir, "fp13.txt"), false) == [1, 3]
    @test Circuitscape.read_focal_points(Int, joinpath(dir, "fp13.txt"), true) == [2, 4]

    # The edge list can double as the focal node list: its value column is
    # not read as node ids (values 5, 7, 9 are not nodes of the graph)
    write(joinpath(dir, "gv.txt"), "1 2 5\n2 3 7\n3 4 9\n")
    @test Circuitscape.read_focal_points(Int, joinpath(dir, "gv.txt"), false; is_graph = true) == [1, 2, 3, 4]
    r = network("graphfp", "gv.txt", "gv.txt")
    @test r[1, 2:end] == [1, 2, 3, 4] && r[2, 5] ≈ 21.0
    write(joinpath(dir, "gv0.txt"), "0 1 5\n1 2 7\n2 3 9\n")
    r = network("graphfp0", "gv0.txt", "gv0.txt")
    @test r[1, 2:end] == [1, 2, 3, 4] && r[2, 5] ≈ 21.0

    # The focal list is sorted and made unique (np.unique)
    write(joinpath(dir, "fpdup.txt"), "3\n1\n3\n")
    @test Circuitscape.read_focal_points(Int, joinpath(dir, "fpdup.txt"), false) == [1, 3]
    r = network("dup", "g1.txt", "fpdup.txt")
    @test size(r) == (3, 3) && r[1, 2:end] == [1, 3] && r[2, 3] ≈ 3.0

    # Comment lines and comma-delimited edge lists are accepted
    write(joinpath(dir, "gc.txt"), "# a comment\n1,2,1\n2,3,2\n# trailing comment\n3,4,3\n")
    i, j, v, from_zero = Circuitscape.load_graph(Int, joinpath(dir, "gc.txt"), Float64)
    @test (i, j, v) == ([1, 2, 3], [2, 3, 4], [1.0, 2.0, 3.0]) && !from_zero
    write(joinpath(dir, "gcw.txt"), "# a comment\n1 2 1\n2 3 2\n3 4 3\n")
    @test Circuitscape.load_graph(Int, joinpath(dir, "gcw.txt"), Float64) ==
          Circuitscape.load_graph(Int, joinpath(dir, "g1.txt"), Float64)
    @test network("comma", "gc.txt", "fp13.txt") ≈ network("ws", "g1.txt", "fp13.txt")
    write(joinpath(dir, "gbad.txt"), "1 2 x\n")
    @test_throws ArgumentError Circuitscape.load_graph(Int, joinpath(dir, "gbad.txt"), Float64)

    # The reclass table applies to edge values, before the conversion to
    # conductances; chained rows do not compound
    write(joinpath(dir, "rc.txt"), "1 10\n")
    r = network("reclass", "g1.txt", "fp13.txt";
                use_reclass_table = "True", reclass_file = joinpath(dir, "rc.txt"))
    @test r[2, 3] ≈ 12.0    # 10 + 2
    write(joinpath(dir, "rc2.txt"), "2 3\n1 2\n")
    r = network("reclass2", "g1.txt", "fp13.txt";
                use_reclass_table = "True", reclass_file = joinpath(dir, "rc2.txt"))
    @test r[2, 3] ≈ 5.0     # 1 -> 2 and 2 -> 3 in one pass

    # A zero resistance edge is refused at load rather than becoming Inf
    write(joinpath(dir, "gz.txt"), "1 2 0\n2 3 1\n")
    @test_throws ArgumentError network("zero", "gz.txt", "fp13.txt")
    d = Circuitscape.init_config()
    d["data_type"] = "network"; d["scenario"] = "pairwise"
    d["habitat_file"] = joinpath(dir, "gz.txt"); d["habitat_map_is_resistances"] = "True"
    d["point_file"] = joinpath(dir, "fp13.txt")
    @test_throws ArgumentError Circuitscape.get_network_data(Float64, Circuitscape.CSConfig(d))
    d["habitat_map_is_resistances"] = "False"   # zero conductance is fine
    @test Circuitscape.get_network_data(Float64, Circuitscape.CSConfig(d)).coords[3] == [0.0, 1.0]

    # Advanced mode: a source or ground at an id that is not in the graph is
    # ignored with a warning. Unit source at 1, direct ground at 3 of the
    # 1 -(1)- 2 -(2)- 3 -(3)- 4 chain: 3 volts across 1 and 3.
    write(joinpath(dir, "src.txt"), "1 1\n9 1\n")
    write(joinpath(dir, "gnd.txt"), "3 0\n7 0\n")
    a = Circuitscape.init_config()
    a["data_type"] = "network"; a["scenario"] = "advanced"
    a["habitat_file"] = joinpath(dir, "g1.txt"); a["habitat_map_is_resistances"] = "True"
    a["source_file"] = joinpath(dir, "src.txt"); a["ground_file"] = joinpath(dir, "gnd.txt")
    a["ground_file_is_resistances"] = "True"; a["output_file"] = joinpath(dir, "adv.out")
    # runtests.jl silences warnings globally; capture them for this call
    Logging.disable_logging(Logging.Debug)
    r = try
        @test_logs (:warn, r"Ignoring source at node 9") (:warn, r"Ignoring ground at node 7") match_mode = :any compute(a)
    finally
        Logging.disable_logging(Logging.Warn)
    end
    @test r[:, 1] == [1, 2, 3, 4]
    @test r[1, 2] - r[3, 2] ≈ 3.0
    @test abs(r[3, 2]) < 1e-10

    # Source and ground ids shift with the graph only: on a 1-based graph a
    # source listed at 0 is not a node and does not shift the other ids
    write(joinpath(dir, "src0.txt"), "0 1\n1 1\n")
    @test Circuitscape.read_point_strengths(Float64, joinpath(dir, "src0.txt"), false)[:, 1] == [0, 1]
    @test Circuitscape.read_point_strengths(Float64, joinpath(dir, "src0.txt"), true)[:, 1] == [1, 2]
    a["source_file"] = joinpath(dir, "src0.txt"); a["output_file"] = joinpath(dir, "adv0.out")
    Logging.disable_logging(Logging.Debug)
    r = try
        @test_logs (:warn, r"Ignoring source at node 0") match_mode = :any compute(a)
    finally
        Logging.disable_logging(Logging.Warn)
    end
    @test r[1, 2] - r[3, 2] ≈ 3.0
    @test abs(r[3, 2]) < 1e-10
end

# Output files, as Circuitscape 4 writes them: the prefix is output_file
# without its extension, the INI copy goes to <prefix>.ini, one-to-all writes
# a resistances file, and the map output options are applied to every grid
# as it is written rather than before accumulation.
@testset "output files" begin
    # The prefix helper: splitext, so a directory containing ".out" is safe
    @test Circuitscape.output_prefix("a/b.out") == "a/b"
    @test Circuitscape.output_prefix("a.out/b.out") == "a.out/b"
    @test Circuitscape.output_prefix("a.out/b") == "a.out/b"
    @test Circuitscape.output_prefix("run") == "run"
    @test Circuitscape.output_prefix("out/res.txt") == "out/res"

    root = mktempdir()
    dir = joinpath(root, "job.out")          # a directory named with ".out"
    mkpath(dir)
    hdr = "ncols 5\nnrows 5\nxllcorner 0\nyllcorner 0\ncellsize 1\nNODATA_value -9999\n"
    grid(f) = hdr * join([join([string(f(i, j)) for j in 1:5], " ") for i in 1:5], "\n") * "\n"
    write(joinpath(dir, "cell.asc"), grid((i, j) -> (i, j) == (3, 3) ? -9999 : 1))
    write(joinpath(dir, "pts.asc"), grid((i, j) -> get(Dict((1,1) => 1, (5,5) => 2, (1,5) => 3), (i, j), 0)))
    read_map(name) = readdlm(joinpath(dir, name), skipstart = 6)
    files(prefix) = filter(f -> startswith(f, prefix), readdir(dir))

    d = Circuitscape.init_config()
    d["data_type"] = "raster"; d["scenario"] = "pairwise"
    d["habitat_file"] = joinpath(dir, "cell.asc"); d["point_file"] = joinpath(dir, "pts.asc")
    d["connect_four_neighbors_only"] = "True"
    d["write_cur_maps"] = "True"; d["write_max_cur_maps"] = "True"; d["write_volt_maps"] = "True"

    # Raw maps: every output lands in the ".out" directory under the prefix
    d["output_file"] = joinpath(dir, "raw.out")
    r_raw = compute(d)
    @test Set(files("raw")) == Set(["raw.ini", "raw_resistances.out", "raw_resistances_3columns.out",
        "raw_cum_curmap.asc", "raw_max_curmap.asc",
        "raw_curmap_1_2.asc", "raw_curmap_1_3.asc", "raw_curmap_2_3.asc",
        "raw_voltmap_1_2.asc", "raw_voltmap_1_3.asc", "raw_voltmap_2_3.asc"])
    @test !isfile(joinpath(dir, "raw.out"))
    ini = read(joinpath(dir, "raw.ini"), String)
    @test occursin("[Circuitscape Mode]", ini) && occursin("scenario = pairwise", ini)
    @test occursin("output_file = $(joinpath(dir, "raw.out"))", ini)
    @test readdlm(joinpath(dir, "raw_resistances.out")) ≈ r_raw
    raw_cum = read_map("raw_cum_curmap.asc")
    raw_max = read_map("raw_max_curmap.asc")
    raw_pair = read_map("raw_curmap_1_2.asc")
    raw_volt = read_map("raw_voltmap_1_2.asc")
    @test raw_cum ≈ raw_pair + read_map("raw_curmap_1_3.asc") + read_map("raw_curmap_2_3.asc")
    @test raw_max ≈ max.(raw_pair, read_map("raw_curmap_1_3.asc"), read_map("raw_curmap_2_3.asc"))
    @test raw_cum[3, 3] == 0 && raw_volt[3, 3] == 0

    # log_transform_maps: log10 of the sum, not the sum of the logs, and
    # cells with no current in a single pair are nodata only in that pair
    logt(m) = map(x -> x > 0 ? log10(x) : -9999.0, m)
    d["output_file"] = joinpath(dir, "log.out"); d["log_transform_maps"] = "True"
    @test compute(d) ≈ r_raw
    @test read_map("log_cum_curmap.asc") ≈ logt(raw_cum)
    @test read_map("log_max_curmap.asc") ≈ logt(raw_max)
    @test read_map("log_curmap_1_2.asc") ≈ logt(raw_pair)
    @test read_map("log_voltmap_1_2.asc") ≈ raw_volt      # voltages are never transformed
    d["log_transform_maps"] = "False"

    # set_null_*_to_nodata: only the NODATA cell (no node) becomes -9999,
    # in every written grid; a dry cell with a node stays 0
    nulls(m) = (n = copy(m); n[3, 3] = -9999.0; n)
    d["output_file"] = joinpath(dir, "null.out")
    d["set_null_currents_to_nodata"] = "True"; d["set_null_voltages_to_nodata"] = "True"
    @test compute(d) ≈ r_raw
    @test read_map("null_cum_curmap.asc") == nulls(raw_cum)
    @test read_map("null_max_curmap.asc") == nulls(raw_max)
    @test read_map("null_curmap_1_2.asc") == nulls(raw_pair)
    @test read_map("null_voltmap_1_2.asc") == nulls(raw_volt)

    # write_cum_cur_map_only: no per-pair current files, cumulative still there
    d["output_file"] = joinpath(dir, "cumonly.out"); d["write_cum_cur_map_only"] = "True"
    compute(d)
    @test !any(f -> occursin("_curmap_", f), files("cumonly"))
    @test read_map("cumonly_cum_curmap.asc") == nulls(raw_cum)
    d["write_cum_cur_map_only"] = "False"

    # One-to-all: a two-column resistances file, the same options applied
    # to the per-node and cumulative grids, and no per-node current maps
    # under write_cum_cur_map_only
    d["scenario"] = "one-to-all"
    d["set_null_currents_to_nodata"] = "False"; d["set_null_voltages_to_nodata"] = "False"
    d["output_file"] = joinpath(dir, "ota.out")
    r_ota = compute(d)
    @test size(r_ota) == (3, 2) && r_ota[:, 1] == [1, 2, 3]
    @test readdlm(joinpath(dir, "ota_resistances.out")) ≈ r_ota
    @test isfile(joinpath(dir, "ota.ini"))
    ota_cum = read_map("ota_cum_curmap.asc")
    ota_1 = read_map("ota_curmap_1.asc")
    @test ota_cum ≈ ota_1 + read_map("ota_curmap_2.asc") + read_map("ota_curmap_3.asc")
    d["output_file"] = joinpath(dir, "otalog.out")
    d["log_transform_maps"] = "True"; d["set_null_currents_to_nodata"] = "True"
    d["set_null_voltages_to_nodata"] = "True"
    @test compute(d) ≈ r_ota
    @test read_map("otalog_cum_curmap.asc") ≈ nulls(logt(ota_cum))
    @test read_map("otalog_curmap_1.asc") ≈ nulls(logt(ota_1))
    @test read_map("otalog_voltmap_1.asc") == nulls(read_map("ota_voltmap_1.asc"))
    d["output_file"] = joinpath(dir, "otacum.out"); d["write_cum_cur_map_only"] = "True"
    @test compute(d) ≈ r_ota
    @test !any(f -> occursin("_curmap_", f), files("otacum"))
    @test read_map("otacum_cum_curmap.asc") ≈ nulls(logt(ota_cum))
    # ... and the cumulative map is still accumulated without write_cur_maps
    d["output_file"] = joinpath(dir, "otacum2.out"); d["write_cur_maps"] = "False"
    compute(d)
    @test read_map("otacum2_cum_curmap.asc") ≈ nulls(logt(ota_cum))

    # Advanced mode: the options apply to its one voltage and current grid,
    # and write_cum_cur_map_only alone still writes (and computes) that grid
    write(joinpath(dir, "src.asc"), grid((i, j) -> (i, j) == (1, 1) ? 1 : 0))
    write(joinpath(dir, "gnd.asc"), grid((i, j) -> (i, j) == (5, 5) ? 1 : -9999))
    a = Circuitscape.init_config()
    a["data_type"] = "raster"; a["scenario"] = "advanced"
    a["habitat_file"] = joinpath(dir, "cell.asc"); a["connect_four_neighbors_only"] = "True"
    a["source_file"] = joinpath(dir, "src.asc"); a["ground_file"] = joinpath(dir, "gnd.asc")
    a["ground_file_is_resistances"] = "False"
    a["write_cur_maps"] = "True"; a["write_volt_maps"] = "True"
    a["output_file"] = joinpath(dir, "adv.out")
    v_adv = compute(a)
    adv_cur = read_map("adv_curmap.asc"); adv_volt = read_map("adv_voltmap.asc")
    @test adv_volt ≈ v_adv
    @test isapprox(adv_cur[1, 1], 1, atol = 1e-6) && adv_cur[3, 3] == 0
    a["output_file"] = joinpath(dir, "advlog.out")
    a["log_transform_maps"] = "True"; a["set_null_currents_to_nodata"] = "True"
    a["set_null_voltages_to_nodata"] = "True"
    @test compute(a) ≈ v_adv
    @test read_map("advlog_curmap.asc") ≈ nulls(logt(adv_cur))
    @test read_map("advlog_voltmap.asc") == nulls(adv_volt)
    a["output_file"] = joinpath(dir, "advcum.out")
    a["write_cur_maps"] = "False"; a["write_cum_cur_map_only"] = "True"
    compute(a)
    @test read_map("advcum_curmap.asc") ≈ nulls(logt(adv_cur))
    a["output_file"] = joinpath(dir, "advnone.out"); a["write_cum_cur_map_only"] = "False"
    compute(a)
    @test !any(f -> occursin("curmap", f), files("advnone"))

    # Network pairwise: per-pair current files only with write_cur_maps and
    # without write_cum_cur_map_only; the prefix helper is used there too
    write(joinpath(dir, "g.txt"), "1\t2\t1.0\n2\t3\t1.0\n3\t4\t1.0\n")
    write(joinpath(dir, "fp.txt"), "1\n4\n")
    n = Circuitscape.init_config()
    n["data_type"] = "network"; n["scenario"] = "pairwise"
    n["habitat_file"] = joinpath(dir, "g.txt"); n["point_file"] = joinpath(dir, "fp.txt")
    n["write_cur_maps"] = "True"; n["write_volt_maps"] = "True"
    n["output_file"] = joinpath(dir, "net.out")
    compute(n)
    @test Set(files("net")) == Set(["net.ini", "net_resistances.out", "net_resistances_3columns.out",
        "net_node_currents_1_4.txt", "net_branch_currents_1_4.txt",
        "net_node_currents_cum.txt", "net_branch_currents_cum.txt", "net_voltages_1_4.txt"])
    n["output_file"] = joinpath(dir, "netcum.out"); n["write_cum_cur_map_only"] = "True"
    compute(n)
    @test Set(files("netcum")) == Set(["netcum.ini", "netcum_resistances.out",
        "netcum_resistances_3columns.out", "netcum_node_currents_cum.txt",
        "netcum_branch_currents_cum.txt", "netcum_voltages_1_4.txt"])
end

# A text-list focal node exactly on the lower-left corner of the habitat map
# maps to column 0, which is outside a 1-based grid: it must be reported as
# such, not fail later with a BoundsError.
@testset "text-list focal node bounds" begin
    dir = mktempdir()
    meta = Circuitscape.RasterMeta(5, 5, 0.0, 0.0, 1.0, -9999.0, [0.0], "")
    write(joinpath(dir, "corner.txt"), "1 0.0 0.0\n2 4.5 4.5\n")
    @test_throws "At least one focal node location falls outside of habitat map" begin
        Circuitscape.read_point_map(Int, joinpath(dir, "corner.txt"), meta)
    end
    write(joinpath(dir, "inside.txt"), "1 0.5 0.5\n2 4.5 4.5\n")
    i, j, v = Circuitscape.read_point_map(Int, joinpath(dir, "inside.txt"), meta)
    @test v == [1, 2] && i == [5, 1] && j == [1, 5]
end

# The INI copy of a run: complete, so it reads back as the configuration that
# was run; written to <prefix>.ini, but never over the INI file being run.
@testset "configuration copy" begin
    dir = mktempdir()
    hdr = "ncols 3\nnrows 3\nxllcorner 0\nyllcorner 0\ncellsize 1\nNODATA_value -9999\n"
    write(joinpath(dir, "cell.asc"), hdr * "1 1 1\n1 1 1\n1 1 1\n")
    write(joinpath(dir, "pts.asc"), hdr * "1 0 0\n0 0 0\n0 0 2\n")
    d = Circuitscape.init_config()
    d["data_type"] = "raster"; d["scenario"] = "pairwise"
    d["habitat_file"] = joinpath(dir, "cell.asc"); d["point_file"] = joinpath(dir, "pts.asc")
    d["output_file"] = joinpath(dir, "job.out")
    d["write_cur_maps"] = "True"; d["set_focal_node_currents_to_zero"] = "True"
    d["precision"] = "single"; d["cholmod_batch_size"] = "7"; d["write_as_tif"] = "False"
    d["log_level"] = "WARNING"
    cfg = Circuitscape.CSConfig(d)
    # `version` records the writer ([Version] section), not an option
    reads_as(path, c) = Circuitscape.CSConfig(Circuitscape.parse_config(path); version = c.version) == c

    # write_config: every option round-trips; the default path is output_file
    Circuitscape.write_config(cfg)
    @test reads_as(joinpath(dir, "job.out"), cfg)
    Circuitscape.write_config(cfg, joinpath(dir, "copy.ini"))
    @test reads_as(joinpath(dir, "copy.ini"), cfg)
    rm(joinpath(dir, "job.out")); rm(joinpath(dir, "copy.ini"))
    # ... from a Dict too, by default to output_file
    Circuitscape.write_config(d)
    @test reads_as(joinpath(dir, "job.out"), cfg)
    Circuitscape.write_config(d, joinpath(dir, "copy.ini"))
    @test reads_as(joinpath(dir, "copy.ini"), cfg)
    rm(joinpath(dir, "job.out")); rm(joinpath(dir, "copy.ini"))

    # compute(dict) writes <prefix>.ini and not output_file
    d["log_level"] = "INFO"
    compute(d)
    @test isfile(joinpath(dir, "job.ini")) && !isfile(joinpath(dir, "job.out"))
    @test reads_as(joinpath(dir, "job.ini"), Circuitscape.CSConfig(d))

    # compute(path) with the INI at <prefix>.ini leaves the user's file alone
    ini = joinpath(dir, "job.ini")
    text = "# my job\n" * join(["$k = $v" for (k, v) in d], "\n") * "\n"
    write(ini, text)
    compute(ini)
    @test read(ini, String) == text
    # ... but an INI elsewhere still gets its copy at <prefix>.ini
    other = joinpath(dir, "other.ini")
    write(other, text)
    rm(ini)
    compute(other)
    @test read(other, String) == text
    @test reads_as(ini, Circuitscape.CSConfig(d))

    # an output_file naming a directory is a configuration error
    d["output_file"] = dir * "/"
    @test_throws "names a directory" Circuitscape.validate(Circuitscape.CSConfig(d))
end
