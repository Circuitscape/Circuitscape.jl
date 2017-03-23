const AAGRID = 2
const TXT = 3

immutable RasterMeta
    ncols::Int64
    nrows::Int64
    xllcorner::Float64
    yllcorner::Float64
    cellsize::Float64
    nodata::Int64
    file_type::Int64
end


function read_graph(a::Inifile, gpath::String)
    i,j,v = load_graph(gpath)
    idx = findfirst(x -> x < 1, i)
    idx != 0 && throw("Indices no good")
    idx = findfirst(x -> x < 1, j)
    idx != 0 && throw("Indices no good")
    is_res = get(a, "Habitat raster or graph", "habitat_map_is_resistances")
    if is_res == "True"
        v = 1./v
    end
    m = max(i[end], j[end])
    A = sparse(i,j,v,m,m)
    A + A'
end

function load_graph(gpath::String)
    g = readdlm(gpath)
    i = Int.(g[:,1]) + 1
    j = Int.(g[:,2]) + 1
    v = Float64.(g[:,3])
    i,j,v
end

read_focal_points(path::String) = Int.(vec(readcsv(path)) + 1)

function read_point_strengths(path::String)
    a = readdlm(path)
    a[:,1] = a[:,1] + 1
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
                gmap[i] =  1 ./ cell_map[i]
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
    g = readdlm(file; skipstart = 6)
    g, rastermeta
end

function _ascii_grid_read_header(habitat_file)
    file_type = _guess_file_type(habitat_file)
    f = open(habitat_file, "r")
    ncols = parse(Int, split(readline(f))[2])
    nrows = parse(Int, split(readline(f))[2])
    xllcorner = float(split(readline(f))[2])
    yllcorner = float(split(readline(f))[2])
    cellsize = float(split(readline(f))[2])
    nodata = parse(Int, split(readline(f))[2])
    RasterMeta(ncols, nrows, xllcorner, yllcorner, cellsize, nodata, file_type)
end

function _guess_file_type(filename) 
    if endswith(filename, ".asc")
        return AAGRID
    elseif endswith(filename, ".txt")
        return TXT
    else
        throw("Check file format")
    end
end

function read_polymap(file, habitatmeta; nodata_as = 0, resample = true)
    #rastermeta = _ascii_grid_read_header(file)
    polymap, rastermeta = _ascii_grid_reader(file)

    ind = find(x -> x == rastermeta.nodata, polymap)
    if rastermeta.nodata != -1
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
    points_rc = filetype == AAGRID ? read_polymap(file, habitatmeta) : read_txt_points(file)
    findnz(points_rc)
end

function read_txt_points(file)
    readdlm(file)
end

