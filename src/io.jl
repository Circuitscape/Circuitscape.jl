abstract type Data end

struct IncludeExcludePairs{V}
    mode::Symbol
    point_ids::Vector{V}
    include_pairs::Matrix{V}
end
function IncludeExcludePairs(V)
    IncludeExcludePairs(:undef, V[], Matrix{V}(undef,0,0))
end

struct NetworkData{T,V} <: Data
    coords::Tuple{Vector{V},Vector{V},Vector{T}}
    fp::Vector{V}
    source_map::Matrix{T}
    ground_map::Matrix{T}
end

struct RasterMeta
    ncols::Int
    nrows::Int
    xllcorner::Float64
    yllcorner::Float64
    cellsize::Float64
    nodata::Float64
    file_type::Int
    crs::String
    affine_map::AffineMap
end
function RasterMeta()
    RasterMeta(0,0,0,0,0,0,0,"",AffineMap([0 0; 0 0], [0, 0]))
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

function load_graph(V, gpath::String, ::Type{T}) where {T}
    g = readdlm(gpath, T)
    i = zeros(V, size(g, 1))
    j = zeros(V, size(g, 1))
    v = zeros(T, size(g, 1))
    for iter = 1:size(g, 1)
        i[iter] = g[iter,1] + 1
        j[iter] = g[iter,2] + 1
        v[iter] = g[iter,3]
    end
    i,j,v
end

read_focal_points(V, path::String) = V.(vec(readdlm(path)) .+ 1)

function read_point_strengths(T, path::String, inc = true)
    a = readdlm(path, T)
    if inc
        @. a[:,1] = a[:,1] + 1
    end
    a
end

function read_cellmap(habitat_file::String, is_res::Bool, ::Type{T}) where {T}

    if is_geotiff(habitat_file)
        cell_map, rastermeta = read_geotiff(T, habitat_file)
    else
        cell_map, rastermeta = _ascii_grid_reader(T, habitat_file)
    end

    gmap = similar(cell_map)
    ind = findall(x -> x == -9999, cell_map)
    if is_res
        if count(x -> x == 0, cell_map) > 0
            throw("Error: zero resistance values are not currently supported for habitat maps. Use a short-circuit region file instead.")
        else
            for i in eachindex(cell_map)
                gmap[i] = 1 ./ cell_map[i]
            end
            gmap[ind] .= 0
        end
    else
        copyto!(gmap, cell_map)
        gmap[ind] .= 0
    end
    gmap, rastermeta
end

function _ascii_grid_reader(T, file)
    f = endswith(file, ".gz") ? GZip.open(file, "r") : open(file, "r")
    rastermeta = _ascii_grid_read_header(file, f)
    c = Matrix{T}(undef,0,0)
    ss = 6
    if rastermeta.nodata == -Inf
        ss = 5
    end
    try
        c = readdlm(f, T; skipstart = ss)
    catch
        seek(f, 0)
        try
            d = readdlm(f; skipstart = ss)
            d = d[:, 1:end-1]
            c = map(T, d)
        catch
            error("Failed to read habitat map. There may be errors in your file.")
        end
    end
    map!(x -> x == rastermeta.nodata ? -9999. : x , c, c)
    close(f)
    c, rastermeta
end

function _ascii_grid_read_header(habitat_file, f)
    file_type = _guess_file_type(habitat_file, f)
    ncols = parse(Int, split(readline(f))[2])
    nrows = parse(Int, split(readline(f))[2])
    xllcorner = parse(Float64, split(readline(f))[2])
    yllcorner = parse(Float64, split(readline(f))[2])
    cellsize = parse(Float64, split(readline(f))[2])
    nodata = -Inf
    s = split(readline(f))
    if occursin("NODATA", s[1]) || occursin("nodata", s[1])
        nodata = parse(Float64, s[2])
    end
    seek(f, 0)
    RasterMeta(ncols, nrows, xllcorner, yllcorner, cellsize, nodata, file_type ,"", AffineMap([0 0; 0 0], [0, 0]))
end

#=function _guess_file_type(filename, f)
    s = readline(f)
    seek(f, 0)

    if startswith(s, "min")
        return PAIRS_AAGRID
    elseif startswith(s, "mode")
        return PAIRS_LIST
    elseif occursin(".asc", filename)
        return AAGRID
    elseif endswith(filename, ".txt")
        return TXTLIST
    else
        throw("Check file format")
    end

end=#

function is_geotiff(filename)
    endswith(filename,"tif") | endswith(filename,"tiff")
end

function read_geotiff(T, f)
    gt=GeoArrays.read(f)
    ncols = size(gt,1)
    nrows = size(gt,2)
    xll = gt.f.translation[1]
    yll= gt.f.translation[2]-gt.f.linear[1,1]-(gt.f.linear[1,1]*(size(gt,2)-1))
    cellsize = gt.f.linear[1,1]
    crs = gt.crs
    affine_map = gt.f
    file_type = FILE_TYPE_GEOTIFF
    nodata = -Inf
    Array{T}(permutedims(gt.A[:,:,1])), RasterMeta(ncols, nrows, xll, yll, cellsize, nodata, file_type,crs, affine_map)
end

function _guess_file_type(filename, f)

    if is_geotiff(filename)
        filetype = FILE_TYPE_GEOTIFF
    else
        hdr = readline(f)
        seek(f, 0)

        if startswith(hdr, FILE_HDR_NPY)
            filetype = FILE_TYPE_NPY
        elseif startswith(lowercase(hdr), FILE_HDR_AAGRID)
            filetype = FILE_TYPE_AAGRID
        elseif startswith(hdr, FILE_HDR_INCL_PAIRS_AAGRID)
            filetype = FILE_TYPE_INCL_PAIRS_AAGRID
        elseif startswith(hdr, FILE_HDR_INCL_PAIRS)
            filetype = FILE_TYPE_INCL_PAIRS
        else
            filetype = FILE_TYPE_TXTLIST
        end
    end

    return filetype
end

function read_polymap(T, file::String, habitatmeta;
                            nodata_as = 0, resample = true)

    if is_geotiff(file)
        polymap, rastermeta = read_geotiff(T, file)
    else
        polymap, rastermeta = _ascii_grid_reader(T, file)
    end

    ind = findall(x -> x == rastermeta.nodata, polymap)
    if nodata_as != -1
        polymap[ind] .= nodata_as
    end

    if(is_geotiff(file))
        if rastermeta.yllcorner != habitatmeta.yllcorner
            cswarn("yllcorner is not the same")
        elseif rastermeta.xllcorner != habitatmeta.xllcorner
            cswarn("xllcorner is not the same")
        elseif rastermeta.ncols != habitatmeta.ncols
            cswarn("ncols is not the same")
        elseif rastermeta.nrows != habitatmeta.nrows
            cswarn("nrows is not the same")
        elseif rastermeta.cellsize != habitatmeta.cellsize
            cswarn("cellsize is not the same")
        elseif rastermeta.crs != habitatmeta.crs
            cswarn("reference system is not the same")
        end
    else
        if rastermeta.cellsize != habitatmeta.cellsize
            cswarn("cellsize is not the same")
        elseif rastermeta.ncols != habitatmeta.ncols
            cswarn("ncols is not the same")
        elseif rastermeta.nrows != habitatmeta.nrows
            cswarn("nrows is not the same")
        elseif rastermeta.yllcorner != habitatmeta.yllcorner
            cswarn("yllcorner is not the same")
        elseif rastermeta.xllcorner != habitatmeta.xllcorner
            cswarn("xllcorner is not the same")
        end
    end

    polymap
end

function read_point_map(V, file, habitatmeta)

    # Advanced mode
    if file == "none"
        return (V[], V[], V[])
    end

    f = endswith(file, ".gz") ? GZip.open(file, "r") : open(file, "r")
    filetype = _guess_file_type(file, f)
    _points_rc = filetype == FILE_TYPE_TXTLIST ? readdlm(file) :
                        read_polymap(V, file, habitatmeta)

    i = V[]
    j = V[]
    v = V[]
    if filetype == FILE_TYPE_TXTLIST
        I = _points_rc[:,2]
        J = _points_rc[:,3]
        v = _points_rc[:,1]
        i  = ceil.(V, habitatmeta.nrows .- (J .- habitatmeta.yllcorner) ./ habitatmeta.cellsize)
        j = ceil.(V, (I .- habitatmeta.xllcorner) ./ habitatmeta.cellsize)
    else
        _I = findall(!iszero, _points_rc)
        (i,j,v) =  (getindex.(_I, 1), getindex.(_I, 2), _points_rc[_I])
    end

    ind = findall(x -> x < 0, v)

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

    if (minimum(i) < 0) || (minimum(j) < 0) ||
            (maximum(i) > (habitatmeta.nrows)) ||
            (maximum(j) > (habitatmeta.ncols))
         throw("At least one focal node location falls outside of habitat map")
    end

    if size(unique(v),1) < 2
        throw("Less than two valid focal nodes found. Please check focal node location file.")
    end

    close(f)
    V.(i), V.(j), V.(v)
end

function read_source_and_ground_maps(T, V, source_file, ground_file, habitatmeta,
                                        is_res)

    ground_map = Matrix{T}(undef,0,0)
    source_map = Matrix{T}(undef,0,0)

    f = endswith(ground_file, "gz") ? Gzip.open(ground_file, "r") : open(ground_file, "r")
    filetype = _guess_file_type(ground_file, f)

    if filetype == FILE_TYPE_AAGRID || filetype == FILE_TYPE_GEOTIFF
        ground_map = read_polymap(T, ground_file, habitatmeta; nodata_as = -1)
        ground_map = map(T, ground_map)
    else
        rc = readdlm(ground_file, V)
        ground_map = -9999 * ones(T, habitatmeta.nrows, habitatmeta.ncols)
        ground_map[rc[:,2], rc[:,3]] = rc[:,1]
    end
    close(f)

    f = endswith(source_file, "gz") ? Gzip.open(source_file, "r") : open(source_file, "r")
    filetype = _guess_file_type(source_file, f)

    if filetype == FILE_TYPE_AAGRID || filetype == FILE_TYPE_GEOTIFF
        source_map = read_polymap(T, source_file, habitatmeta)
        source_map = map(T, source_map)
    else
        rc = readdlm(source_file, V)
        source_map = -9999 * ones(T, habitatmeta.nrows, habitatmeta.ncols)
        source_map[rc[:,2], rc[:,3]] = rc[:,1]
    end
    close(f)

    if is_res
        ind = findall(x -> x == -9999, ground_map)
        ground_map = 1 ./ ground_map
        ground_map[ind] .= 0
    else
        ind = findall(x -> x == -9999, ground_map)
        ground_map[ind] .= 0
    end

    source_map, ground_map
end

function read_included_pairs(V, filename)

    f = endswith(filename, "gz") ? Gzip.open(filename, "r") : open(filename, "r")
    filetype = _guess_file_type(filename, f)
    minval = 0
    maxval = 0
    mode = :undef

    if filetype == FILE_TYPE_INCL_PAIRS_AAGRID
        open(filename, "r") do fV
            minval = parse(Float64, split(readline(f))[2])
            maxval = parse(Float64, split(readline(f))[2])
        end
        included_pairs = readdlm(filename, skipstart=2)
        point_ids = V.(included_pairs[:,1])
        deleteat!(point_ids, 1)
        included_pairs = included_pairs[2:end, 2:end]
        map!(x -> x > maxval ? 0 : x, included_pairs, included_pairs)
        idx = findall(x -> x >= minval, included_pairs)
        mode = :include
        bin = map(x -> x >= minval ? V(1) : V(0), included_pairs)
        return IncludeExcludePairs(mode, point_ids, bin)
    elseif filetype == FILE_TYPE_INCL_PAIRS
        open(filename, "r") do f
            mode = Symbol(split(readline(f))[2])
        end
        included_pairs = readdlm(filename, skipstart = 1)
        if size(included_pairs, 1) == 1
            pl = zeros(V, 1,2)
            pl[1,:] = included_pairs
            included_pairs = pl
        end
        point_ids = V.(sort!(unique(included_pairs)))
        idx = findall(x -> x == 0 , point_ids)
        if length(idx) > 0
            deleteat!(point_ids, idx)
            cswarn("Code to include pairs is activated, some entries did not match with focal node file. Some focal nodes may have been dropped")
        end

        mat = zeros(V, size(point_ids, 1), size(point_ids, 1))

        for i = 1:size(included_pairs, 1)
            idx1 = findfirst(x -> x == included_pairs[i, 1], point_ids)
            idx2 = findfirst(x -> x == included_pairs[i, 2], point_ids)
            if idx1 != nothing && idx2 != nothing
                mat[idx1,idx2] = 1
                mat[idx2,idx1] = 1
            end
        end

        return IncludeExcludePairs(mode, point_ids, mat)
    end

end

function get_network_data(T, V, cfg)::NetworkData{T,V}

    hab_is_res = cfg["habitat_map_is_resistances"] in TRUELIST
    hab_file = cfg["habitat_file"]
    fp_file = cfg["point_file"]
    source_file = cfg["source_file"]
    ground_file = cfg["ground_file"]

    is_pairwise = cfg["scenario"] in PAIRWISE

    i,j,v = load_graph(V, hab_file, T)
    if hab_is_res
        v = 1 ./ v
    end

    if is_pairwise
        fp = read_focal_points(V, fp_file)
    else
        fp = V[]
    end

    if !is_pairwise
        source_map = read_point_strengths(T, source_file)
        ground_map = read_point_strengths(T, ground_file)
    else
        source_map = Matrix{T}(undef,0,0)
        ground_map = Matrix{T}(undef,0,0)
    end

    NetworkData((i,j,v), fp, source_map, ground_map)
end

function load_raster_data(T, V, cfg)::RasData{T,V}

    # Habitat file
    hab_file = cfg["habitat_file"]
    hab_is_res = cfg["habitat_map_is_resistances"] in TRUELIST

    # Polygons
    use_polygons = cfg["use_polygons"] in TRUELIST
    polygon_file = cfg["polygon_file"]

    # Mask file
    use_mask = cfg["use_mask"] in TRUELIST
    mask_file = cfg["mask_file"]

    # Point file
    point_file = cfg["point_file"]

    # Variable source strengths
    use_var_source = cfg["use_variable_source_strengths"] in TRUELIST
    var_source_file = cfg["variable_source_file"]

    # Included Pairs
    use_inc_pairs = cfg["use_included_pairs"] in TRUELIST
    inc_pairs_file = cfg["included_pairs_file"]

    # Advanced mode
    is_pairwise = cfg["scenario"] in PAIRWISE
    is_advanced = cfg["scenario"] in ADVANCED
    source_file = cfg["source_file"]
    ground_file = cfg["ground_file"]
    ground_is_res = cfg["ground_file_is_resistances"] in TRUELIST

    csinfo("Reading maps")

    # Read cell map
    cellmap, hbmeta = read_cellmap(hab_file, hab_is_res, T)
    c = count(x -> x > 0, cellmap)
    csinfo("Resistance/Conductance map has $c nodes")

    # Read polymap
    if use_polygons
        polymap = read_polymap(V, polygon_file, hbmeta)
    else
        polymap = Matrix{V}(undef,0,0)
    end

    # Read and update cellmap with mask file
    if use_mask
        update!(cellmap, mask_file, hbmeta)
        sum(cellmap) == 0 && throw("Mask file deleted everything!")
    end

    # Read point file
    if !is_advanced
        points_rc = read_point_map(V, point_file, hbmeta)
    else
        points_rc = (V[], V[], V[])
    end

    # Advanced mode reading
    if is_advanced
        source_map, ground_map =
        read_source_and_ground_maps(T, V, source_file, ground_file,
                                    hbmeta, ground_is_res)
    else
        source_map, ground_map = Matrix{T}(undef,0,0), Matrix{T}(undef,0,0)
    end

    # Included Pairs
    if use_inc_pairs
        included_pairs = read_included_pairs(V, inc_pairs_file)
    else
        included_pairs = IncludeExcludePairs(V)
    end

    # Variable source strengths
    if use_var_source
        strengths = read_point_strengths(T, var_source_file)
    else
        strengths = Matrix{T}(undef, 0,0)
    end

    RasData(cellmap, polymap, source_map, ground_map, points_rc, strengths,
                    included_pairs, hbmeta)
end

function update!(cellmap::Matrix{T}, m::String, hbmeta) where {T}
    mask = read_polymap(T, m, hbmeta)
    map!(x -> x > 0 ? 1 : 0, mask, mask)
    cellmap .= cellmap .* mask
end
