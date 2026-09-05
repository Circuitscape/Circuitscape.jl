# read_point_map must drop every negative focal-node value, not shift indices.
# Deleting one index at a time from a precomputed index list makes every later
# index off by one, so with two or more negatives the wrong nodes were removed.

import Circuitscape: read_point_map, RasterMeta
using Circuitscape, Test

let
    dir = mktempdir()
    meta = RasterMeta(5, 5, 0.0, 0.0, 1.0, -9999.0, [0.0], "")

    # ASCII grid. findall walks column-major, so the nonzero values are met in
    # the order 1, -1, 2, -2, 3 -> negatives at positions 2 and 4.
    pts = joinpath(dir, "pts.asc")
    write(pts, """ncols         5
nrows         5
xllcorner     0
yllcorner     0
cellsize      1
NODATA_value  -9999
0 -1 0 0 0
1 0 0 0 0
0 0 0 -2 0
0 0 2 0 0
0 0 0 0 3
""")
    i, j, v = read_point_map(Int, pts, meta)
    @test v == [1, 2, 3]
    @test i == [2, 4, 5]
    @test j == [1, 3, 5]

    # Text list (value x y) in the same order and at the same cells. With
    # xll = yll = 0 and cellsize 1, row r / col c is at x = c - 0.5, y = 5 - r + 0.5.
    txt = joinpath(dir, "pts.txt")
    write(txt, """1 0.5 3.5
-1 1.5 4.5
2 2.5 1.5
-2 3.5 2.5
3 4.5 0.5
""")
    i, j, v = read_point_map(Int, txt, meta)
    @test v == [1, 2, 3]
    @test i == [2, 4, 5]
    @test j == [1, 3, 5]

    # Negative entry last: the shifted index used to run past the end (BoundsError).
    txt2 = joinpath(dir, "pts2.txt")
    write(txt2, """1 0.5 3.5
2 2.5 1.5
-1 1.5 4.5
3 4.5 0.5
-2 3.5 2.5
""")
    i, j, v = read_point_map(Int, txt2, meta)
    @test v == [1, 2, 3]
    @test i == [2, 4, 5]
    @test j == [1, 3, 5]
end
