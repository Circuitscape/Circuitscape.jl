const AAGRID = 2
const TXTLIST = 3
const PAIRS_AAGRID = 4
const PAIRS_LIST = 5
const truelist = ["True", "true"]

abstract type Data end
abstract type Flags end
abstract type Polygon end
abstract type VarSource end
abstract type IncludedPairs end
abstract type Mask end
abstract type InputFlags <: Flags end
abstract type ComputeFlags <: Flags end
abstract type PointFilePolygons end

struct PointFileContainsPolygons <: PointFilePolygons
end
struct PointFileNoPolygons <: PointFilePolygons
end

struct UsePoly <: Polygon
    file::String
end
struct NoPoly <: Polygon
end
struct UseVarSrc <: VarSource
    file::String
end
struct NoVarSrc <: VarSource
end
struct UseIncPairs <: IncludedPairs
    file::String
end
struct NoIncPairs <: IncludedPairs
end
struct UseMask <: Mask
    file::String
end
struct NoMask <: Mask
end

struct IncludeExcludePairs
    mode::Symbol
    point_ids::Vector{Int64}
    include_pairs::Matrix{Int64}
end
function IncludeExcludePairs()
    IncludeExcludePairs(:undef, Int64[], Matrix{Int64}(0, 0))
end

struct NetPairFlags{T} <: InputFlags
    precision::T
    hab_is_res::Bool
    hab_file::String
    fp_file::String
end
struct NetAdvFlags{T} <: InputFlags
    precision::T
    hab_is_res::Bool
    hab_file::String
    source_file::String
    ground_file::String
end

struct RasInputFlags{T,P,M,V,I} <: InputFlags
    precision::T
    hab_is_res::Bool
    hab_file::String
    poly::P
    mask::M
    point_file::String
    var_source::V
    included_pairs::I
    source_file::String
    ground_file::String
    ground_is_res::Bool
end

struct RasCompFlags{P} <: ComputeFlags
    four_neighbors::Bool
    avg_res::Bool
    ptpoly::P
end

struct NetPairData{Ti,Tv} <: Data
    A::SparseMatrixCSC{Ti,Tv}
    fp::Vector{Tv}
end

struct NetAdvData{Ti,Tv} <: Data
    A::SparseMatrixCSC{Ti,Tv}
    source_map::Matrix{Ti}
    ground_map::Matrix{Ti}
end

struct RasData{T,V} <: Data
    cellmap::Matrix{T}
    polymap::Matrix{V}
    source_map::Matrix{V}
    ground_map::Matrix{V}
    points_rc::Tuple{Vector{V},Vector{V},Vector{V}}
    strengths::Matrix{T}
    included_pairs::IncludeExcludePairs
end

struct RasterMeta
    ncols::Int64
    nrows::Int64
    xllcorner::Float64
    yllcorner::Float64
    cellsize::Float64
    nodata::Float64
    file_type::Int64
end
function RasterMeta()
    RasterMeta(0,0,0,0,0,0,0)
end

function read_graph(a, gpath::String)
    i,j,v = load_graph(gpath, T)
    idx = findfirst(x -> x < 1, i)
    idx != 0 && throw("Indices no good")
    idx = findfirst(x -> x < 1, j)
    idx != 0 && throw("Indices no good")
    is_res = a["habitat_map_is_resistances"]
    if is_res == "True"
        v = 1./v
    end
    m = max(i[end], j[end])
    A = sparse(i,j,v,m,m)
    A + A'
end

function read_graph{T}(is_res::Bool, gpath::String, ::Type{T})
    i,j,v = load_graph(gpath, T)
    idx = findfirst(x -> x < 1, i)
    idx != 0 && throw("Indices no good")
    idx = findfirst(x -> x < 1, j)
    idx != 0 && throw("Indices no good")
    if is_res
        v = 1./v
    end
    m = max(i[end], j[end])
    A = sparse(i,j,v,m,m)
    A + A'
end

function load_graph{T}(gpath::String, ::Type{T})
    g = readdlm(gpath, T)
    i = zeros(Int, size(g, 1))
    j = zeros(Int, size(g, 1))
    v = zeros(T, size(g, 1))
    for iter = 1:size(g, 1)
        i[iter] = g[iter,1] + 1
        j[iter] = g[iter,2] + 1
        v[iter] = g[iter,3]
    end
    i,j,v
end

read_focal_points(path::String) = Int.(vec(readcsv(path)) + 1)

function read_point_strengths{T}(path::String, ::Type{T}, inc = true)
    a = readdlm(path, T)
    if inc
        a[:,1] = a[:,1] + 1
    end
    a
end

function read_cellmap{T}(habitat_file::String, is_res::Bool, ::Type{T})

    cell_map, rastermeta = _ascii_grid_reader(habitat_file, T)

    gmap = similar(cell_map)
    ind = find(x -> x == -9999, cell_map)
    if is_res
        if count(x -> x == 0, cell_map) > 0
            throw("Error: zero resistance values are not currently supported for habitat maps. Use a short-circuit region file instead.")
        else
            for i in eachindex(cell_map)
                gmap[i] = 1./cell_map[i]
            end
            gmap[ind] = 0
        end
    else
        copy!(gmap, cell_map)
        gmap[ind] = 0
    end
    gmap, rastermeta
end

function _ascii_grid_reader{T}(file, ::Type{T})
    f = endswith(file, ".gz") ? GZip.open(file, "r") : open(file, "r")
    rastermeta = _ascii_grid_read_header(file, f)
    c = Matrix{T}(0,0)
    ss = 6
    if rastermeta.nodata == -Inf
        ss = 5
    end
    try
        c = readdlm(f, T; skipstart = ss)
    catch
        seek(f, 0)
        d = readdlm(f; skipstart = ss)
        d = d[:, 1:end-1]
        c = map(T, d)
    end
    map!(x -> x == rastermeta.nodata ? -9999. : x , c, c)
    c, rastermeta
end

function _ascii_grid_read_header(habitat_file, f)
    file_type = _guess_file_type(habitat_file, f)
    ncols = parse(Int, split(readline(f))[2])
    nrows = parse(Int, split(readline(f))[2])
    xllcorner = float(split(readline(f))[2])
    yllcorner = float(split(readline(f))[2])
    cellsize = float(split(readline(f))[2])
    nodata = -Inf
    s = split(readline(f))
    if contains(s[1], "NODATA") || contains(s[1], "nodata")
        nodata = float(s[2])
    end
    seek(f, 0)
    RasterMeta(ncols, nrows, xllcorner, yllcorner, cellsize, nodata, file_type)
end

function _guess_file_type(filename, f)
    s = readline(f)
    seek(f, 0)

    if startswith(s, "min")
        return PAIRS_AAGRID
    elseif startswith(s, "mode")
        return PAIRS_LIST
    elseif contains(filename, ".asc")
        return AAGRID
    elseif endswith(filename, ".txt")
        return TXTLIST
    else
        throw("Check file format")
    end

end

function read_polymap{T}(file::String, habitatmeta, ::Type{T};
                            nodata_as = 0, resample = true)
    polymap, rastermeta = _ascii_grid_reader(file, T)

    ind = find(x -> x == rastermeta.nodata, polymap)
    if nodata_as != -1
        polymap[ind] = nodata_as
    end

    if rastermeta.cellsize != habitatmeta.cellsize
        warn("cellsize is not the same")
    elseif rastermeta.ncols != habitatmeta.ncols
        warn("ncols is not the same")
    elseif rastermeta.nrows != habitatmeta.nrows
        warn("nrows is not the same")
    elseif rastermeta.yllcorner != habitatmeta.yllcorner
        warn("yllcorner is not the same")
    elseif rastermeta.xllcorner != habitatmeta.xllcorner
        warn("xllcorner is not the same")
    end

    polymap
end

function read_point_map(::Raster{Pairwise}, file, habitatmeta)
    f = endswith(file, ".gz") ? GZip.open(file, "r") : open(file, "r")
    filetype = _guess_file_type(file, f)
    _points_rc = filetype == TXTLIST ? readdlm(file, Int64) :
                        read_polymap(file, habitatmeta, Int64)

    i = Int64[]
    j = Int64[]
    v = Int64[]
    if filetype == TXTLIST
        I = _points_rc[:,2]
        J = _points_rc[:,3]
        v = _points_rc[:,1]
        i  = ceil.(Int, habitatmeta.nrows - (J - habitatmeta.yllcorner) / habitatmeta.cellsize)
        j = ceil.(Int, (I - habitatmeta.xllcorner) / habitatmeta.cellsize)
    else
        (i,j,v) = findnz(_points_rc)
    end

    ind = find(x -> x < 0, v)

    # Get rid of negative resistances
    for index in ind
        deleteat!(i, index)
        deleteat!(j, index)
        deleteat!(v, index)
    end

    # Sort them
    idx = sortperm(v)
    i = i[idx]
    j = j[idx]
    v = v[idx]

    i, j, v
end

function read_source_and_ground_maps{T}(source_file, ground_file, habitatmeta,
                                        is_res, ::Type{T})

    ground_map = Matrix{T}(0,0)
    source_map = Matrix{T}(0,0)

    f = endswith(ground_file, "gz") ? Gzip.open(ground_file, "r") : open(ground_file, "r")
    filetype = _guess_file_type(ground_file, f)

    if filetype == AAGRID
        ground_map = read_polymap(ground_file, habitatmeta, T; nodata_as = -1)
        ground_map = map(Float64, ground_map)
    else
        rc = readdlm(ground_file, Int)
        ground_map = -9999 * ones(T, habitatmeta.nrows, habitatmeta.ncols)
        ground_map[rc[:,2], rc[:,3]] = rc[:,1]
    end

    f = endswith(source_file, "gz") ? Gzip.open(source_file, "r") : open(source_file, "r")
    filetype = _guess_file_type(source_file, f)

    if filetype == AAGRID
        source_map = read_polymap(source_file, habitatmeta, T)
        source_map = map(T, source_map)
    else
        rc = readdlm(source_file, Int)
        source_map = -9999 * ones(T, habitatmeta.nrows, habitatmeta.ncols)
        source_map[rc[:,2], rc[:,3]] = rc[:,1]
    end

    if is_res
        ind = find(x -> x == -9999, ground_map)
        ground_map = 1 ./ ground_map
        ground_map[ind] = 0
    else
        ind = find(x -> x == -9999, ground_map)
        ground_map[ind] = 0
    end

    source_map, ground_map
end

function read_included_pairs(file)

    f = endswith(file, "gz") ? Gzip.open(file, "r") : open(file, "r")
    filetype = _guess_file_type(file, f)
    minval = 0
    maxval = 0
    mode = :undef

    if filetype == PAIRS_AAGRID
        open(file, "r") do f
            minval = float(split(readline(f))[2])
            maxval = float(split(readline(f))[2])
        end
        included_pairs = readdlm(file, skipstart=2)
        point_ids = Int.(included_pairs[:,1])
        deleteat!(point_ids, 1)
        included_pairs = included_pairs[2:end, 2:end]
        map!(x -> x > maxval ? 0 : x, included_pairs, included_pairs)
        idx = find(x -> x >= minval, included_pairs)
        mode = :include
        bin = map(x -> x >= minval ? 1 : 0, included_pairs)
        IncludeExcludePairs(mode, point_ids, bin)
    else
        open(file, "r") do f
            mode = Symbol(split(readline(f))[2])
        end
        included_pairs = readdlm(file, skipstart = 1)
        point_ids = Int.(sort!(unique(included_pairs)))
        if point_ids[1] == 0
            deleteat!(point_ids, 1)
        end

        mat = zeros(size(point_ids, 1), size(point_ids, 1))

        for i = 1:size(included_pairs, 1)
            idx1 = findfirst(x -> x == included_pairs[i, 1], point_ids)
            idx2 = findfirst(x -> x == included_pairs[i, 2], point_ids)
            if idx1 != 0 && idx2 != 0
                mat[idx1,idx2] = 1
                mat[idx2,idx1] = 1
            end
        end
        IncludeExcludePairs(mode, point_ids, mat)
    end
end

function inputflags{S}(obj::Raster{S}, cfg)

    # Precision for computation
    p = cfg["precision"] == "Single" ? Float32 : Float64

    # Habitat file
    hab_file = cfg["habitat_file"]
    hab_is_res = cfg["habitat_map_is_resistances"] in truelist

    # Polygons
    poly = cfg["use_polygons"] in truelist ?
                    UsePoly(cfg["polygon_file"]) : NoPoly()

    # Mask file
    mask = cfg["use_mask"] in truelist ? UseMask(cfg["mask_file"]) : NoMask()

    # Point file
    point_file = cfg["point_file"]

    # Variable source strengths
    var_source = cfg["use_variable_source_strengths"] in truelist ?
                    UseVarSrc(cfg["variable_source_file"]) : NoVarSrc()

    included_pairs = cfg["use_included_pairs"] in truelist ?
                        UseIncPairs(cfg["included_pairs_file"]) : NoIncPairs()

    source_file = cfg["source_file"]
    ground_file = cfg["ground_file"]
    ground_is_res = cfg["ground_file_is_resistances"] in truelist

    RasInputFlags(p, hab_is_res, hab_file, poly, mask, point_file, var_source,
                included_pairs, source_file, ground_file, ground_is_res)
end

function inputflags{S}(obj::Network{S}, cfg)
    p = cfg["precision"] == "Single" ? Float32 : Float64
    hab_is_res = cfg["habitat_map_is_resistances"] in truelist
    hab_file = cfg["habitat_file"]
    _ipflags(obj, p, hab_is_res, hab_file, cfg)
end
function _ipflags(::Network{Pairwise}, p, hab_is_res, hab_file, cfg)
    fp_file = cfg["point_file"]
    NetPairFlags(p, hab_is_res, hab_file, fp_file)
end
function _ipflags(::Network{Advanced}, p, hab_is_res, hab_file, cfg)
    source_file = cfg["source_file"]
    ground_file = cfg["ground_file"]
    NetAdvFlags(p, hab_is_res, hab_file, source_file, ground_file)
end
function grab_input{S}(obj::Network{S}, flags)
    A = read_graph(flags.hab_is_res, flags.hab_file, flags.precision)
    _grab_input(obj, A, flags)
end
function _grab_input(::Network{Pairwise}, A, flags)
    fp = read_focal_points(flags.fp_file)
    NetPairData(A, fp)
end
function _grab_input(::Network{Advanced}, A, flags)
    source_map = read_point_strengths(flags.source_file, flags.precision)
    ground_map = read_point_strengths(flags.ground_file, flags.precision)
    NetAdvData(A, source_map, ground_map)
end
function grab_input{S}(obj::Raster{S}, flags)

    # Precision
    p = flags.precision

    info("Reading maps")

    # Read cell map
    cellmap, hbmeta = read_cellmap(flags.hab_file, flags.hab_is_res, p)
    c = count(x -> x > 0, cellmap)
    info("Resistance/Conductance map has $c nodes")

    # Read polymap
    polymap = read_polymap(flags.poly, hbmeta, p)

    # Read and update cellmap with mask file
    update!(cellmap, flags.mask, hbmeta)
    sum(cellmap) == 0 && throw("Mask file deleted everything!")

    # Read point file
    points_rc = read_point_map(obj, flags.point_file, hbmeta)

    # Advanced mode reading
    source_map, ground_map = advanced_read(obj, p, hbmeta, flags)

    # Included Pairs
    included_pairs = read_included_pairs(obj, flags)

    # Variable source strengths
    strengths = read_point_strengths(obj, flags)

    RasterData(cellmap, polymap, source_map, ground_map, points_rc, strengths,
                    included_pairs), hbmeta
end
read_polymap{T}(::NoPoly, hbmeta, ::Type{T}) = Matrix{T}(0,0)
read_polymap{T}(p::UsePoly, hbmeta, ::Type{T}) = read_polymap(p.file, hbmeta, T)
update!(cellmap, ::NoMask, hbmeta) = cellmap
function update!{T}(cellmap::Matrix{T}, m::Mask, hbmeta)
    mask = read_polymap(m.file, hbmeta, T)
    map!(x -> x > 0 ? 1 : 0, mask, mask)
    cellmap .= cellmap .* mask
end
read_point_map(::Raster{Advanced}, x, y) = (Int[], Int[], Int[])
read_point_map(::Raster{OneToAll}, args...) = read_point_map(Raster{Pairwise}(), args...)
read_point_map(::Raster{AllToOne}, args...) = read_point_map(Raster{Pairwise}(), args...)
function advanced_read{T<:Union{Pairwise,AllToOne,OneToAll},V}(::Raster{T}, ::Type{V}, args...)
    Matrix{V}(0,0), Matrix{V}(0,0)
end
function advanced_read{T}(::Raster{Advanced}, ::Type{T}, hbmeta, flags)
    read_source_and_ground_maps(flags.source_file, flags.ground_file,
                                hbmeta, flags.ground_is_res, T)
end
read_included_pairs(::Raster{Advanced}, args...) = IncludeExcludePairs()
read_included_pairs{T<:Union{Pairwise,AllToOne,OneToAll}}(::Raster{T}, flags) =
    read_included_pairs(flags.included_pairs)
read_included_pairs(p::UseIncPairs) = read_included_pairs(p.file)
read_included_pairs(::NoIncPairs) = IncludeExcludePairs()
read_point_strengths{T<:Union{OneToAll,AllToOne}}(::Raster{T}, flags) =
    read_point_strengths(flags.var_source, flags.precision)
read_point_strengths{T<:Union{Pairwise,Advanced}}(::Raster{T}, flags) =
    read_point_strengths(NoVarSrc(), flags.precision)
read_point_strengths{T}(::NoVarSrc, ::Type{T}) = Matrix{T}(0,0)
read_point_strengths{T}(v::UseVarSrc, ::Type{T}) = read_point_strengths(v.file, T, false)

function computeflags{S}(::Raster{S}, cfg, points_rc)
    four_neighbors = cfg["connect_four_neighbors_only"] == "True"
    avg_res = cfg["connect_using_avg_resistances"] == "True"
    ptpoly = length(points_rc[1]) != length(unique(points_rc[3])) ?
                PointFileContainsPolygons() : PointFileNoPolygons()
    RasCompFlags(four_neighbors, avg_res, ptpoly)
end
