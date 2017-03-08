gmap = [1 0 1
        0 1 0
        1 0 1]

polymap = [1 0 1
           0 2 0
           2 0 0]

r = CircuitScape.construct_node_map(gmap, polymap)

@test r == [1 0 1
            0 2 0
            2 0 3]
