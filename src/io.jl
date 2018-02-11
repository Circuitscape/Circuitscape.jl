abstract type Data end 

struct IncludeExcludePairs{T}
    mode::Symbol
    point_ids::Vector{T}
    include_pairs::Matrix{T}
end
function IncludeExcludePairs()
    IncludeExcludePairs(:undef, Int[], Matrix{Int}(0, 0))
end

struct RasInputFlags{T}
    is_advanced::Bool
    precision::T
    hab_is_res::Bool
    hab_file::String
    use_polygon::Bool
    polygon::String
    use_mask::Bool
    mask::String
    point_file::String
    use_var_src_str::Bool
    var_source::String
    use_included_pairs::Bool
    included_pairs::String
    source_file::String
    ground_file::String
    ground_is_res::Bool
end

struct RasCompFlags
    four_neighbors::Bool
    avg_res::Bool
    ptpoly::Bool
end

struct NetworkData{T,V} <: Data
    coords::Tuple{Vector{T},Vector{T},Vector{V}}
    fp::Vector{T}
    source_map::Matrix{V}
    ground_map::Matrix{V}
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

struct RasData{T,V} <: Data
    cellmap::Matrix{T}
    polymap::Matrix{V}
    source_map::Matrix{T}
    ground_map::Matrix{T}
    points_rc::Tuple{Vector{V},Vector{V},Vector{V}}
    strengths::Matrix{T}
    included_pairs::IncludeExcludePairs{V}
    hbmeta::RasterMeta
end

function read_graph(is_res::Bool, gpath::String, ::Type{T}) where {T}
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

function load_graph(gpath::String, ::Type{T}) where {T}
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

function read_point_strengths(path::String, ::Type{T}, inc = true) where {T}
    a = readdlm(path, T)
    if inc
        a[:,1] = a[:,1] + 1
    end
    a
end

function read_cellmap(habitat_file::String, is_res::Bool, ::Type{T}) where {T}

    cell_map, rastermeta = _ascii_grid_reader(T, habitat_file)

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

function _ascii_grid_reader(T, file)
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

function read_polymap(file::String, habitatmeta;
                            nodata_as = 0, resample = true)

    polymap, rastermeta = _ascii_grid_reader(Int64, file)

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

function read_point_map(file, habitatmeta)

    # Advanced mode 
    if file == "none" 
        return (Int[], Int[], Int[])
    end

    f = endswith(file, ".gz") ? GZip.open(file, "r") : open(file, "r")
    filetype = _guess_file_type(file, f)
    _points_rc = filetype == TXTLIST ? readdlm(file, Int64) :
                        read_polymap(file, habitatmeta)

    i = Int[]
    j = Int[]
    v = Int[]
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

    Int.(i), Int.(j), Int.(v)
end

function read_source_and_ground_maps(source_file, ground_file, habitatmeta,
                                        is_res, ::Type{T}) where {T}

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

function get_network_data(T, cfg)
    
    hab_is_res = cfg["habitat_map_is_resistances"] in truelist
    hab_file = cfg["habitat_file"]
    fp_file = cfg["point_file"]
    source_file = cfg["source_file"]
    ground_file = cfg["ground_file"]

    is_pairwise = cfg["scenario"] in PAIRWISE

    I,J,V = load_graph(hab_file, T)
    if hab_is_res
        V = 1./V
    end

    fp = read_focal_points(fp_file)

    if !is_pairwise
        source_map = read_point_strengths(source_file, T)
        ground_map = read_point_strengths(ground_file, T)
    else
        source_map = Matrix{T}(0,0)
        ground_map = Matrix{T}(0,0)
    end

    NetworkData((I,J,V), fp, source_map, ground_map)
end

function load_raster_data(T, cfg)

    # Habitat file
    hab_file = cfg["habitat_file"]
    hab_is_res = cfg["habitat_map_is_resistances"] in truelist

    # Polygons
    use_polygons = cfg["use_polygons"] in truelist
    polygon_file = cfg["polygon_file"]

    # Mask file
    use_mask = cfg["use_mask"] in truelist 
    mask_file = cfg["mask_file"]

    # Point file
    point_file = cfg["point_file"]

    # Variable source strengths
    use_var_source = cfg["use_variable_source_strengths"] in truelist 
    var_source_file = cfg["variable_source_file"]

    # Included Pairs
    use_inc_pairs = cfg["use_included_pairs"] in truelist 
    inc_pairs_file = ["included_pairs_file"]

    # Advanced mode
    is_pairwise = cfg["scenario"] in PAIRWISE
    source_file = cfg["source_file"]
    ground_file = cfg["ground_file"]
    ground_is_res = cfg["ground_file_is_resistances"] in truelist
    
    csinfo("Reading maps")

    # Read cell map
    cellmap, hbmeta = read_cellmap(hab_file, hab_is_res, T)
    c = count(x -> x > 0, cellmap)
    csinfo("Resistance/Conductance map has $c nodes")

    # Read polymap
    if use_polygons
        polymap = read_polymap(polygon_file, hbmeta)
    else
        polymap = Matrix{Int}(0,0)
    end

    # Read and update cellmap with mask file
    if use_mask
        update!(cellmap, mask_file, hbmeta)
        sum(cellmap) == 0 && throw("Mask file deleted everything!")
    end

    # Read point file
    if is_pairwise
        points_rc = read_point_map(point_file, hbmeta)
    else
        points_rc = (Int64[], Int64[], Int64[]) 
    end

    # Advanced mode reading
    if !is_pairwise
        source_map, ground_map = advanced_read(T, hbmeta, flags)
    else
        source_map, ground_map = Matrix{T}(0,0), Matrix{T}(0,0)
    end

    # Included Pairs
    if use_inc_pairs
        included_pairs = read_included_pairs(flags)
    else
        included_pairs = IncludeExcludePairs()
    end

    # Variable source strengths
    if use_var_source
        strengths = read_point_strengths(flags)
    else
        strengths = Matrix{T}(0,0)
    end
    
    RasData(cellmap, polymap, source_map, ground_map, points_rc, strengths,
                    included_pairs, hbmeta)
end

function inputflags_raster(cfg)

    # Advanced mode
    is_advanced = cfg["scenario"] in ADVANCED

    # Precision for computation
    p = cfg["precision"] == "Single" ? Float32 : Float64

    # Habitat file
    hab_file = cfg["habitat_file"]
    hab_is_res = cfg["habitat_map_is_resistances"] in truelist

    # Polygons
    use_poly = cfg["use_polygons"] in truelist
    poly = use_poly ? cfg["polygon_file"] : "none"

    # Mask file
    use_mask = cfg["use_mask"] in truelist
    mask = use_mask ? cfg["mask_file"] : "none"

    # Point file
    point_file = cfg["scenario"] == "advanced" ? "none" : cfg["point_file"]

    # Variable source strengths
    use_var_src_str = cfg["use_variable_source_strengths"] in truelist
    var_source = use_var_src_str ? cfg["variable_source_file"] : "none"

    # Included Pairs
    use_included_pairs = cfg["use_included_pairs"] in truelist && 
                            cfg["scenario"] != "advanced"
    included_pairs = use_included_pairs ?
                            cfg["included_pairs_file"] : "none"

    source_file = cfg["source_file"]
    ground_file = cfg["ground_file"]
    ground_is_res = cfg["ground_file_is_resistances"] in truelist

    RasInputFlags(is_advanced,
                p, hab_is_res, hab_file, 
                use_poly, poly, # Polygons 
                use_mask, mask, # Mask
                point_file, 
                use_var_src_str, var_source, # Variable Source Strengths
                use_included_pairs, included_pairs, # Included Pairs
                source_file, ground_file, ground_is_res)
end

function grab_input_raster(flags)

    # Precision
    p = flags.precision

    info("Reading maps")

    # Read cell map
    cellmap, hbmeta = read_cellmap(flags.hab_file, flags.hab_is_res, p)
    c = count(x -> x > 0, cellmap)
    info("Resistance/Conductance map has $c nodes")

    # Read polymap
    if flags.use_polygon
        polymap = read_polymap(flags.polygon, hbmeta, Int32)
    else
        polymap = empty_matrix(Int32)
    end

    # Read and update cellmap with mask file
    if flags.use_mask
        update!(cellmap, flags.mask, hbmeta)
    end
    sum(cellmap) == 0 && throw("Mask file deleted everything!")

    # Read point file
    if flags.is_advanced
        points_rc = (Int32[], Int32[], Int32[])
    else
        points_rc = read_point_map(flags.point_file, hbmeta)
    end

    # Advanced mode reading
    if flags.is_advanced
        source_map, ground_map = advanced_read(p, hbmeta, flags)
    else
        source_map, ground_map = empty_matrix(p), empty_matrix(p)
    end

    # Included Pairs
    if flags.use_included_pairs
        included_pairs = read_included_pairs(obj, flags)
    else
        included_pairs = IncludeExcludePairs()
    end

    # Variable source strengths
    if flags.use_var_src_str
        strengths = read_point_strengths(obj, flags)
    else
        strengths = empty_matrix(p)
    end

    RasData(cellmap, polymap, source_map, ground_map, 
            points_rc, strengths, included_pairs), hbmeta
end

function update!(cellmap::Matrix{T}, m::String, hbmeta) where {T}
    mask = read_polymap(m, hbmeta, T)
    map!(x -> x > 0 ? 1 : 0, mask, mask)
    cellmap .= cellmap .* mask
end