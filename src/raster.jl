function compute_raster(cfg::Inifile)

    # Read inputs
    gmap, polymap, points_rc = load_maps(cfg)
    c = count(x -> x > 0, gmap)
    info("Resistance/Conductance map has $c nodes")
    four_neighbors = get(cfg, "Connection scheme for raster habitat data",
                                "connect_four_neighbors_only") == "True"
    average_resistances = get(cfg, "Connection scheme for raster habitat data",
                                "connect_using_avg_resistances") == "True"

    resistances = pairwise_module(gmap, polymap, points_rc, four_neighbors, 
                                    average_resistances)

    #gmap, polymap, points_rc
end

function load_maps(cfg::Inifile)

    # Read raster map
    info("Reading Maps")
    habitat_file = get(cfg, "Habitat raster or graph", "habitat_file")
    is_res = get(cfg, "Habitat raster or graph", "habitat_map_is_resistances") == "True"

    cellmap, habitatmeta = read_cell_map(habitat_file, is_res)

    # Read polygon map
    use_polygons = get(cfg, "Short circuit regions (aka polygons)", "use_polygons") == "True"
    polymap_file = get(cfg, "Short circuit regions (aka polygons)", "polygon_file")
    polymap = use_polygons ? read_polymap(polymap_file, habitatmeta) : Array{Float64,2}()

    scenario = get(cfg, "Circuitscape Mode", "scenario")
    point_file = get(cfg, "Options for pairwise and one-to-all and all-to-one modes",
                            "point_file")

    points_rc = (Vector{Int}(), Vector{Int}(), Vector{Float64}())
    if scenario == "advanced"
    else
        points_rc = read_point_map(point_file, habitatmeta)
    end

    cellmap, polymap, points_rc
end

function pairwise_module(gmap, polymap, points_rc, four_neighbors, average_resistances)

    point_file_contains_polygons = length(points_rc[1]) != length(unique(points_rc[3]))
    f = average_resistances ? res_avg : cond_avg

    #if !point_file_contains_polygons
        nodemap = construct_node_map(gmap, polymap)
        I = Int64[]
        J = Int64[]
        V = Float64[]
        for j = 1:size(gmap, 2)
            for i = 1:size(gmap, 1)
                if nodemap[i,j] == 0
                    continue
                else
                    # Horizontal neighbour
                    if j != size(gmap, 2) && nodemap[i,j+1] != 0
                        push!(I, nodemap[i,j])
                        push!(J, nodemap[i,j+1])
                        push!(V, f(gmap[i,j], gmap[i,j+1]))
                    end

                    # Vertical neighbour
                    if i != size(gmap, 1) && nodemap[i+1, j] != 0
                        push!(I, nodemap[i,j])
                        push!(J, nodemap[i+1,j])
                        push!(V, f(gmap[i,j], gmap[i+1,j]))
                    end

                    if !four_neighbors
                        # Diagonal neighbour
                        if i != size(gmap, 1) && j != size(gmap, 2) && nodemap[i+1, j+1] != 0
                            push!(I, nodemap[i,j])
                            push!(J, nodemap[i+1,j+1])
                            push!(V, weird_avg(gmap[i,j], gmap[i+1,j+1]))
                        end
                        
                        if i != 1 && j != size(gmap, 2) && nodemap[i-1, j+1] != 0
                            push!(I, nodemap[i,j])
                            push!(J, nodemap[i-1,j+1])
                            push!(V, weird_avg(gmap[i,j], gmap[i-1,j+1]))
                        end
                    end
                end
            end
        end
        m = maximum(nodemap)
        a = sparse(I,J,V, m, m)
        a = a + a'
        g = Graph(a)

        c = zeros(Int, length(points_rc[3]))
        for (i,v) in enumerate(zip(points_rc[1], points_rc[2]))
            c[i] = nodemap[v...]
        end

        resistances = single_ground_all_pair_resistances(a, g, c)
        return resistances
    #end
    return nothing
end

res_avg(x, y) = 1 / ((1/x + 1/y) / 2)
cond_avg(x, y) = (x + y) / 2
weird_avg(x,y) = (x + y) / (2*√2)

function construct_node_map(gmap, polymap)
    nodemap = zeros(size(gmap)) 
    if isempty(polymap)
         ind = find(x -> x > 0, gmap)
         nodemap[ind] = 1:length(ind)
    else
        d = Dict{Int, Vector{Int}}()
        for i in unique(polymap)
            d[i] = find(x -> x == i, polymap)
        end
        k = 1
        for i in find(gmap)
            if i in d[0]
                nodemap[i] = k
                k += 1
            else
                for key in keys(d)
                    if i in d[key]
                        if i == first(d[key])
                            nodemap[i] += k
                            k += 1
                        else
                            nodemap[i] = nodemap[first(d[key])]
                        end
                    end
                end
            end
        end
    end

    nodemap 
end
