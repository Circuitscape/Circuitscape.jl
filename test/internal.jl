gmap = [1 0 1
        0 1 0
        1 0 1]

polymap = [1 0 1
           0 2 0
           2 0 0]

r = Circuitscape.construct_node_map(gmap, polymap)

@test r == [1 0 1
            0 2 0
            2 0 3]

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

    nodemap = Circuitscape.construct_node_map(gmap, polymap)

    @test nodemap == [ 3.0  18.0  0.0  10.0  14.0
                       0.0   0.0  0.0  11.0  15.0
                       1.0   4.0  7.0  12.0  16.0
                       2.0   5.0  8.0  13.0  17.0
                       3.0   6.0  9.0   0.0  18.0]
end

let 

    cfg = Circuitscape.parse_config("input/raster/one_to_all/11/oneToAllVerify11.ini")
    r, hbmeta = Circuitscape.load_maps(cfg)

    cellmap = r.cellmap
    polymap = r.polymap
    points_rc = r.points_rc
    point_map = [ 1.0  2.0  0.0  0.0  0.0
                  0.0  0.0  0.0  0.0  0.0
                  3.0  0.0  0.0  7.0  0.0
                  4.0  0.0  0.0  0.0  0.0
                  1.0  0.0  0.0  0.0  2.0 ]

    r = Circuitscape.create_new_polymap(cellmap, polymap, points_rc, point_map = point_map)

    @test r == [ 1.0  2.0  0.0  0.0  0.0
                 0.0  0.0  0.0  0.0  0.0
                 12.0  0.0  0.0  2.0  0.0
                 1.0  0.0  0.0  0.0  0.0
                 1.0  0.0  0.0  0.0  2.0 ]
end

import Circuitscape: resolve_conflicts

@test resolve_conflicts([1,0,0], [1,0,0], :rmvgnd) == ([1, 0, 0], [0, 0, 0], [1, 0, 0])

@test resolve_conflicts([1,0,0], [1,0,0], :rmvsrc) == ([0, 0, 0], [1, 0, 0], [1, 0, 0])

@test resolve_conflicts([1,0,0], [1,0,0], :keepall) == ([1, 0, 0], [1, 0, 0], [1, 0, 0])

@test resolve_conflicts([1,0,0], [1,0,0], :rmvall) == ([0, 0, 0], [1, 0, 0], [1, 0, 0])