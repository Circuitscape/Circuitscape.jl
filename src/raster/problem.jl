# The raster side of problem construction: how a loaded RasterData becomes
# the graph, focal nodes and cumulative maps of a GraphProblem or
# AdvancedProblem (see problem.jl), and the pairwise loop over focal regions,
# which needs a graph per pair. Grid-to-graph construction itself is in
# graph.jl.

"""
    build_graph(data::RasterData, cfg, polymap = data.polymap) -> (G, cc, geometry)

Laplacian, connected components and geometry of the landscape graph. A
`polymap` other than the loaded one lets callers merge extra cells into
nodes, as the focal-region paths do.
"""
function build_graph(data::RasterData{T,V}, cfg, polymap = data.polymap) where {T,V}
    cellmap = data.cellmap
    nodemap = construct_node_map(cellmap, polymap)
    A = construct_graph(cellmap, nodemap, cfg.connect_using_avg_resistances,
                        cfg.connect_four_neighbors_only)
    cc = connected_components(A)
    G = laplacian!(A)
    G, cc, RasterGeometry(nodemap, polymap, data.hbmeta, cellmap)
end

# Focal nodes are the nodes under the focal cells, paired with their ids.
function focal_nodes(data::RasterData{T,V}, geometry::RasterGeometry) where {T,V}
    points_rc = data.points_rc
    included_pairs = data.included_pairs
    exclude_pairs = isempty(included_pairs) ? Tuple{V,V}[] :
                        generate_exclude_pairs(points_rc, included_pairs)
    nodemap = geometry.nodemap
    points = V[nodemap[i, j] for (i, j) in zip(points_rc[1], points_rc[2])]
    points, points_rc[3], exclude_pairs
end

initialize_cum(data::RasterData, cfg, num_nodes) =
    initialize_cum_maps(data.cellmap, cfg.write_max_cur_maps)

# A point file that names one id in several cells describes focal regions.
has_focal_regions(data::RasterData) =
    length(data.points_rc[1]) != length(unique(data.points_rc[3]))

"""
    pairwise_regions(rasterdata, cfg)

Pairwise mode over focal regions: every pair gets its own graph, with the
two regions being solved collapsed into nodes and the others left as
ordinary habitat.
"""
function pairwise_regions(rasterdata::RasterData{T,V}, cfg)::Matrix{T} where {T,V}

    # get unique list of points
    # for every point pair do
        # construct new polymap
        # construct new nodemap
        # construct new graph
        # solve for two points
    # end
	

    # Data
    gmap = rasterdata.cellmap
    polymap = rasterdata.polymap
    points_rc = rasterdata.points_rc
	included_pairs = rasterdata.included_pairs
	exclude_pairs = isempty(included_pairs) ? Vector{Tuple{V,V}}() : 
						generate_exclude_pairs(points_rc, included_pairs)

    # Cumulative maps
    cum = initialize_cum_maps(gmap, cfg.write_max_cur_maps)

    pts = unique(points_rc[3])
    resistances = -1 * ones(length(pts), length(pts))

    pairs = pairs_to_solve(pts, exclude_pairs)
    n = length(pairs)
    @info("Total number of pair solves = $n")

    for (k, (i, j)) in enumerate(pairs)
        pt1, pt2 = pts[i], pts[j]
        @info("Solving pair $k of $n")
        graphdata = region_pair_problem(rasterdata, pt1, pt2, cum, cfg)
        pairwise_resistance = single_ground_all_pairs(graphdata, cfg, false)
        resistances[i,j] = resistances[j,i] = pairwise_resistance[2,3]
    end
    for i = 1:size(pts, 1)
        resistances[i,i] = 0
    end
    P = [0, pts...]
    r = hcat(P, vcat(pts', resistances))

    write_cum_maps(cum, rasterdata.hbmeta, cfg)

    # save resistances
    save_resistances(r, cfg)

    r
end

"""
    pairs_to_solve(pts, exclude_pairs)

Index pairs `(i, j)` with `i < j` into `pts` whose node IDs are not listed in
`exclude_pairs` (in either order). This is exactly the set of pairs the
polygon path solves, so counting it keeps the reported number of pair solves
consistent with the include/exclude file (issue #341).
"""
function pairs_to_solve(pts, exclude_pairs)
    excluded = Set(exclude_pairs)
    [(i, j) for i in 1:length(pts) for j in i+1:length(pts)
        if (pts[i], pts[j]) ∉ excluded && (pts[j], pts[i]) ∉ excluded]
end


function region_pair_problem(rasterdata::RasterData{T,V},
                            pt1, pt2, cum, cfg)::GraphProblem{T,V} where {T,V}

    points_rc = rasterdata.points_rc

    # Construct new polymap and the graph over it
    newpoly = create_new_polymap(rasterdata.cellmap, rasterdata.polymap, points_rc, pt1, pt2)
    G, cc, geometry = build_graph(rasterdata, cfg, newpoly)

    # Construct points vector
    nodemap = geometry.nodemap
    x = something(findfirst(isequal(pt1), points_rc[3]), 0)
    y = something(findfirst(isequal(pt2), points_rc[3]), 0)
    c1 = nodemap[points_rc[1][x], points_rc[2][x]]
    c2 = nodemap[points_rc[1][y], points_rc[2][y]]
    points = V[c1, c2]

    # Exclude pairs array
    exclude_pairs = Tuple{V,V}[]

    GraphProblem(G, cc, points, [pt1, pt2], exclude_pairs, geometry, cum, get_solver(cfg))
end

function generate_exclude_pairs(points_rc, included_pairs::IncludeExcludePairs{V}) where V

    # Include mode prunes focal nodes not mentioned in the include file;
    # exclude mode keeps every node and only drops the listed pairs.
    if included_pairs.mode == :include
        prune_points!(points_rc, included_pairs.point_ids)
    end

    exclude_pairs_from(included_pairs)
end
