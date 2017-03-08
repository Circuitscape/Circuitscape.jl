function compute_raster(cfg::Inifile)

    # Read inputs
    gmap, polymap, points_rc = load_maps(cfg)
    c = count(x -> x > 0, gmap)
    info("Resistance/Conductance map has $c nodes")

    resistances = pairwise_module(gmap, polymap, points_rc)

    #gmap, polymap, points_rc
end

function load_maps(cfg::Inifile)

    # Read raster map
    info("Reading Maps")
    habitat_file = get(cfg, "Habitat raster or graph", "habitat_file")
    is_res = get(cfg, "Habitat raster or graph", "habitat_file_is_resistances") == "True"
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

function pairwise_module(gmap, polymap, points_rc)

    point_file_contains_polygons = length(points_rc[1]) != length(unique(points_rc[3]))

    if !point_file_contains_polygons
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
                        push!(V, res_avg(gmap[i,j], gmap[i,j+1]))
                    end

                    # Vertical neighbour
                    if i != size(gmap, 1) && nodemap[i+1, j] != 0
                        push!(I, nodemap[i,j])
                        push!(J, nodemap[i+1,j])
                        push!(V, res_avg(gmap[i,j], gmap[i+1,j]))
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
    end
    return nothing
end

res_avg(x, y) = 1 / ((1/x + 1/y) / 2)
function construct_node_map(gmap, polymap)
   nodemap = zeros(size(gmap)) 
   ind = find(gmap) # Replace by more sophisticated logic later
   nodemap[ind] = 1:length(ind)
   ind_sc = find(polymap)
   for i  = 2:length(ind_sc)
       nodemap[ind_sc[i]] = nodemap[ind_sc[i-1]]
   end
   ind_notsc = setdiff(ind, ind_sc)  
   nodemap[ind_notsc] -= 1
   nodemap 
end
