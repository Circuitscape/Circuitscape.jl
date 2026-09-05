"""
    Geometry

How graph nodes map back onto the user's coordinates when maps are written.
Once loaded, a raster landscape and a network are the same problem: a
symmetric Laplacian over graph nodes, its connected components and focal
nodes. The only thing that differs downstream of loading is where a node's
voltage or current goes in the output, and that is what a `Geometry` knows.
Output code dispatches on it rather than on the data type in the config.
"""
abstract type Geometry end

"""
    RasterGeometry(nodemap, polymap, hbmeta, cellmap)

Nodes are cells of a grid: `nodemap[i,j]` is the node of cell `(i,j)`, zero
for cells outside the habitat, and cells of one short-circuit polygon share a
node. `hbmeta` gives the grid's extent and georeferencing, `cellmap` the
conductance per cell (needed to mark null cells as NODATA).
"""
struct RasterGeometry{T,V} <: Geometry
    nodemap::Matrix{V}
    polymap::Matrix{V}
    hbmeta::RasterMeta
    cellmap::Matrix{T}
end

"""
    NetworkGeometry(nodes)

Nodes are the user's own node ids: row `k` of the matrix being solved is node
`nodes[k]`. For the whole graph this is the identity; restricted to a
connected component it is the component's node list. Branch coordinates live
on `Cumulative`.
"""
struct NetworkGeometry{V} <: Geometry
    nodes::Vector{V}
end

is_raster(::RasterGeometry) = true
is_raster(::NetworkGeometry) = false

"""
    restrict(geometry, component)

The geometry of one connected component, whose matrix rows are numbered
`1:length(component)`: a raster nodemap relabelled to those local indices, or
the component's node ids.
"""
restrict(g::RasterGeometry, component) =
    RasterGeometry(construct_local_node_map(g.nodemap, component, g.polymap),
                   g.polymap, g.hbmeta, g.cellmap)
restrict(::NetworkGeometry, component) = NetworkGeometry(component)

# A map to accumulate per-cell values into; networks have no map.
alloc_map(g::RasterGeometry) = zeros(g.hbmeta.nrows, g.hbmeta.ncols)
alloc_map(::NetworkGeometry) = zeros(0, 0)
