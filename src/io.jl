const AAGRID = 2
const TXTLIST = 3
const PAIRS_AAGRID = 4
const PAIRS_LIST = 5

immutable RasterMeta
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

immutable IncludeExcludePairs
    mode::Symbol
    point_ids::Vector{Int64}
    include_pairs::Matrix{Int64}
end
function IncludeExcludePairs()
    IncludeExcludePairs(:undef, Int64[], Matrix{Int64}())
end

function read_graph(a, gpath::String)
    i,j,v = load_graph(gpath)
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

function load_graph(gpath::String)
    g = readdlm(gpath)
    i = zeros(Int, size(g, 1))
    j = zeros(Int, size(g, 1))
    v = zeros(size(g, 1))
    for iter = 1:size(g, 1)
        i[iter] = g[iter,1] + 1
        j[iter] = g[iter,2] + 1
        v[iter] = g[iter,3]
    end
    i,j,v
end

read_focal_points(path::String) = Int.(vec(readcsv(path)) + 1)

function read_point_strengths(path::String, inc = true)
    a = readdlm(path)
    if inc 
        a[:,1] = a[:,1] + 1
    end
    a
end

function read_cell_map(habitat_file, is_res)
    cell_map, rastermeta = _ascii_grid_reader(habitat_file)

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

function _ascii_grid_reader(file)
    rastermeta = _ascii_grid_read_header(file)
    c = Array{Float64,2}()
    ss = 6
    if rastermeta.nodata == -Inf
        ss = 5
    end
    try
        c = readdlm(file, Float64; skipstart = ss)
    catch
        c = readdlm(file; skipstart = ss)
        c = c[:, 1:end-1]
        map!(Float64, c)
    end
    map!(x -> x == rastermeta.nodata ? -9999. : x , c)
    c, rastermeta
end

function _ascii_grid_read_header(habitat_file)
    file_type = _guess_file_type(habitat_file)
    f = open(habitat_file, "r")
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
    RasterMeta(ncols, nrows, xllcorner, yllcorner, cellsize, nodata, file_type)
end

function _guess_file_type(filename) 
    f = open(filename, "r")
    s = readline(f)
    close(f)

    if startswith(s, "min")
        return PAIRS_AAGRID
    elseif startswith(s, "mode")
        return PAIRS_LIST
    elseif endswith(filename, ".asc")
        return AAGRID
    elseif endswith(filename, ".txt")
        return TXTLIST
    else
        throw("Check file format")
    end
end

function read_polymap(file, habitatmeta; nodata_as = 0, resample = true)
    polymap, rastermeta = _ascii_grid_reader(file)

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
    filetype = _guess_file_type(file)
    points_rc = filetype == TXTLIST ? readdlm(file) : read_polymap(file, habitatmeta)

    i = Int64[]
    j = Int64[]
    v = Int64[]
    if filetype == TXTLIST
        I = points_rc[:,2]
        J = points_rc[:,3]
        v = points_rc[:,1]
        i  = ceil(Int, habitatmeta.nrows - (J - habitatmeta.yllcorner) / habitatmeta.cellsize)
        j = ceil(Int, (I - habitatmeta.xllcorner) / habitatmeta.cellsize)
    else
        (i,j,v) = findnz(points_rc)
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

function read_source_and_ground_maps(source_file, ground_file, habitatmeta, is_res)

    ground_map = Array{Float64,2}()
    source_map = Array{Float64,2}()

    filetype = _guess_file_type(ground_file)

    if filetype == AAGRID
        ground_map = read_polymap(ground_file, habitatmeta; nodata_as = -1)
        ground_map = map(Float64, ground_map)
    else
        rc = readdlm(ground_file, Int)
        ground_map = -9999 * ones(habitatmeta.nrows, habitatmeta.ncols)
        ground_map[rc[:,2], rc[:,3]] = rc[:,1]
    end

    filetype = _guess_file_type(source_file)

    if filetype == AAGRID
        source_map = read_polymap(source_file, habitatmeta)
        source_map = map(Float64, source_map)
    else
        rc = readdlm(source_file, Int)
        source_map = -9999 * ones(habitatmeta.nrows, habitatmeta.ncols)
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

    filetype = _guess_file_type(file)
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
        map!(x -> x > maxval ? 0 : x, included_pairs)
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
