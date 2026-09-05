"""
    load_data(T, V, cfg)

Read the inputs named in `cfg`: a `RasterData` or a `NetworkData`. This and
the geometry construction in `build_graph` are the only steps that know
which of the two it is; everything downstream works on one representation.
"""
load_data(T, V, cfg) = is_raster(cfg) ? load_raster_data(T, V, cfg) :
                                        get_network_data(T, V, cfg)

"""
    build_problem(data, cfg)
    build_problem(T, V, cfg)

The problem `cfg` describes over the loaded (or, given `T, V`, freshly read)
data: a `GraphProblem` in pairwise mode, an `AdvancedProblem` in advanced
mode. Both are built the same way for rasters and networks; `build_graph`,
`focal_nodes` and `initialize_cum` dispatch on the data type to construct the
graph and its output geometry.
"""
build_problem(T, V, cfg) = build_problem(load_data(T, V, cfg), cfg)
build_problem(data::Data, cfg) = is_advanced(cfg) ? advanced_problem(data, cfg) :
                                                    graph_problem(data, cfg)

function graph_problem(data::Data, cfg)
    G, cc, geometry = @timeit CSTIMER[] "construct graph" build_graph(data, cfg)
    points, user_points, exclude_pairs = focal_nodes(data, geometry)
    cum = initialize_cum(data, cfg, size(G, 1))
    GraphProblem(G, cc, points, user_points, exclude_pairs, geometry, cum, get_solver(cfg))
end

function advanced_problem(data::Data, cfg)
    G, cc, geometry = @timeit CSTIMER[] "construct graph" build_graph(data, cfg)
    sources, grounds, finitegrounds =
        sources_and_grounds(geometry, data.source_map, data.ground_map, G, cfg)
    V = eltype(rowvals(G))
    AdvancedProblem(G, cc, geometry, sources, grounds, data.source_map,
                    finitegrounds, V(-1), V(0), get_solver(cfg))
end

"""
    run_pairwise(data, cfg)
    run_pairwise(problem::GraphProblem, cfg)

Pairwise mode: solve every pair of focal nodes and write the cumulative
maps. A raster whose focal nodes are regions rebuilds the graph per pair
(`pairwise_regions`); everything else solves one shared problem.
"""
run_pairwise(data::Data, cfg) =
    has_focal_regions(data) ? pairwise_regions(data, cfg) :
                              run_pairwise(build_problem(data, cfg), cfg)

function run_pairwise(prob::GraphProblem{T}, cfg)::Matrix{T} where T
    r = @timeit CSTIMER[] "solve pairwise resistances" single_ground_all_pairs(prob, cfg)
    @timeit CSTIMER[] "write cumulative current maps" write_cum_maps(prob.cum, prob.geometry, cfg)
    r
end

"""
    run_advanced(problem::AdvancedProblem, cfg)

Advanced mode: one solve with the user's sources and grounds.
"""
function run_advanced(prob::AdvancedProblem{T}, cfg)::Matrix{T} where T
    v, _ = advanced_kernel(prob, cfg)
    v
end
