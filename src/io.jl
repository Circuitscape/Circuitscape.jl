abstract type Data end 

struct IncludeExcludePairs{T}
    mode::Symbol
    point_ids::Vector{T}
    include_pairs::Matrix{T}
end
function IncludeExcludePairs()
    IncludeExcludePairs(:undef, INT[], Matrix{INT}(0, 0))
end

struct NetworkData{T,V} <: Data
    coords::Tuple{Vector{V},Vector{V},Vector{T}}
    fp::Vector{V}
    source_map::Matrix{T}
    ground_map::Matrix{T}
end

struct RasterMeta
    ncols::INT
    nrows::INT
    xllcorner::Float64
    yllcorner::Float64
    cellsize::Float64
    nodata::Float64
    file_type::INT
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
    i = zeros(INT, size(g, 1))
    j = zeros(INT, size(g, 1))
    v = zeros(T, size(g, 1))
    for iter = 1:size(g, 1)
        i[iter] = g[iter,1] + 1
        j[iter] = g[iter,2] + 1
        v[iter] = g[iter,3]
    end
    i,j,v
end

read_focal_points(path::String) = INT.(vec(readcsv(path)) + 1)

function read_point_strengths(T, path::String, inc = true)
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
    ncols = parse(INT, split(readline(f))[2])
    nrows = parse(INT, split(readline(f))[2])
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

function read_polymap(T, file::String, habitatmeta;
                            nodata_as = 0, resample = true)

    polymap, rastermeta = _ascii_grid_reader(T, file)

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
        return (INT[], INT[], INT[])
    end

    f = endswith(file, ".gz") ? GZip.open(file, "r") : open(file, "r")
    filetype = _guess_file_type(file, f)
    _points_rc = filetype == TXTLIST ? readdlm(file, INT) :
                        read_polymap(INT, file, habitatmeta)

    i = INT[]
    j = INT[]
    v = INT[]
    if filetype == TXTLIST
        I = _points_rc[:,2]
        J = _points_rc[:,3]
        v = _points_rc[:,1]
        i  = ceil.(INT, habitatmeta.nrows - (J - habitatmeta.yllcorner) / habitatmeta.cellsize)
        j = ceil.(INT, (I - habitatmeta.xllcorner) / habitatmeta.cellsize)
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

    INT.(i), INT.(j), INT.(v)
end

function read_source_and_ground_maps(T, source_file, ground_file, habitatmeta,
                                        is_res) 

    ground_map = Matrix{T}(0,0)
    source_map = Matrix{T}(0,0)

    f = endswith(ground_file, "gz") ? Gzip.open(ground_file, "r") : open(ground_file, "r")
    filetype = _guess_file_type(ground_file, f)

    if filetype == AAGRID
        ground_map = read_polymap(T, ground_file, habitatmeta; nodata_as = -1)
        ground_map = map(T, ground_map)
    else
        rc = readdlm(ground_file, INT)
        ground_map = -9999 * ones(T, habitatmeta.nrows, habitatmeta.ncols)
        ground_map[rc[:,2], rc[:,3]] = rc[:,1]
    end

    f = endswith(source_file, "gz") ? Gzip.open(source_file, "r") : open(source_file, "r")
    filetype = _guess_file_type(source_file, f)

    if filetype == AAGRID
        source_map = read_polymap(T, source_file, habitatmeta)
        source_map = map(T, source_map)
    else
        rc = readdlm(source_file, INT)
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
        point_ids = INT.(included_pairs[:,1])
        deleteat!(point_ids, 1)
        included_pairs = included_pairs[2:end, 2:end]
        map!(x -> x > maxval ? 0 : x, included_pairs, included_pairs)
        idx = find(x -> x >= minval, included_pairs)
        mode = :include
        bin = map(x -> x >= minval ? INT(1) : INT(0), included_pairs)
        IncludeExcludePairs(mode, point_ids, bin)
    else
        open(file, "r") do f
            mode = Symbol(split(readline(f))[2])
        end
        included_pairs = readdlm(file, skipstart = 1)
        point_ids = INT.(sort!(unique(included_pairs)))
        if point_ids[1] == 0
            deleteat!(point_ids, 1)
        end

        mat = zeros(INT, size(point_ids, 1), size(point_ids, 1))

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

function get_network_data(T, cfg)::NetworkData{T,INT}
    
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

    if is_pairwise
        fp = read_focal_points(fp_file)
    else
        fp = INT[]
    end

    if !is_pairwise
        source_map = read_point_strengths(T, source_file)
        ground_map = read_point_strengths(T, ground_file)
    else
        source_map = Matrix{T}(0,0)
        ground_map = Matrix{T}(0,0)
    end

    NetworkData((I,J,V), fp, source_map, ground_map)
end

function load_raster_data(T, cfg)::RasData{T,INT}

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
    inc_pairs_file = cfg["included_pairs_file"]

    # Advanced mode
    is_pairwise = cfg["scenario"] in PAIRWISE
    is_advanced = cfg["scenario"] in ADVANCED
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
        polymap = read_polymap(INT, polygon_file, hbmeta)
    else
        polymap = Matrix{INT}(0,0)
    end

    # Read and update cellmap with mask file
    if use_mask
        update!(cellmap, mask_file, hbmeta)
        sum(cellmap) == 0 && throw("Mask file deleted everything!")
    end

    # Read point file
    if !is_advanced
        points_rc = read_point_map(point_file, hbmeta)
    else
        points_rc = (INT[], INT[], INT[]) 
    end

    # Advanced mode reading
    if is_advanced
        source_map, ground_map = 
        read_source_and_ground_maps(T, source_file, ground_file,
                                    hbmeta, ground_is_res)
    else
        source_map, ground_map = Matrix{T}(0,0), Matrix{T}(0,0)
    end

    # Included Pairs
    if use_inc_pairs
        included_pairs = read_included_pairs(inc_pairs_file)
    else
        included_pairs = IncludeExcludePairs()
    end

    # Variable source strengths
    if use_var_source
        strengths = read_point_strengths(T, var_source_file)
    else
        strengths = Matrix{T}(0,0)
    end
    
    RasData(cellmap, polymap, source_map, ground_map, points_rc, strengths,
                    included_pairs, hbmeta)
end

function update!(cellmap::Matrix{T}, m::String, hbmeta) where {T}
    mask = read_polymap(T, m, hbmeta)
    map!(x -> x > 0 ? 1 : 0, mask, mask)
    cellmap .= cellmap .* mask
end
