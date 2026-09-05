module CircuitscapeArchGDALExt

using ArchGDAL
using Circuitscape
import Circuitscape: read_raster_gdal, write_raster_gdal, IO_LOCK

# GeoTIFF (and any other GDAL-readable) raster I/O. ESRI ASCII grids are
# handled natively in Circuitscape and never reach this extension.

# Inspired by GeoArrays.read()
function read_raster_gdal(path::String, T::Type)
    endswith(path, ".gz") && (path = "/vsigzip/$(path)")
    ArchGDAL.read(path) do raw
        # Extract 1st band (should only be one band anyway)
        # to get a 2D array instead of 3D
        band = ArchGDAL.getband(raw, 1)

        # Extract the array
        array_t = ArchGDAL.read(band)

        # This handles UInt tiff rasters that can still have negative NoData values
        # Need to convert the NoData value to Int64 in these cases
        if eltype(array_t) <: Integer
            ras_type = Int64
        else
            ras_type = eltype(array_t)
        end

        # Extract no data value, first converting it to the proper type (based on
        # the raster). Then, need to convert to T. Weird, yes,
        # but it's the only way I could get it to work for all raster types... -VL
        nodata_val = convert(T, convert(ras_type, ArchGDAL.getnodatavalue(band)))

        # Transpose the array -- ArchGDAL returns a x by y array, need y by x
        array = convert(Array{T}, permutedims(array_t, [2, 1]))

        array[array .== nodata_val] .= -9999.0

        # Line to handle NaNs in datasets read from tifs
        array[isnan.(array)] .= -9999.0

        transform = ArchGDAL.getgeotransform(raw)
        wkt = ArchGDAL.getproj(raw)
        array, wkt, transform # wkt and transform are needed later for write_raster
    end
end

# Inspired by GeoArrays.write()
function write_raster_gdal(fn_prefix::String, array::AbstractMatrix, wkt::String, transform::AbstractVector)
    # Prepare data outside the lock
    array_t = permutedims(array, [2, 1])
    width, height = size(array_t)
    fn = fn_prefix * ".tif"
    options = ["COMPRESS=LZW"]

    # Lock only the ArchGDAL calls (GDAL is not thread safe)
    lock(IO_LOCK) do
        ArchGDAL.create(fn_prefix,
                        driver = ArchGDAL.getdriver("MEM"),
                        width = width,
                        height = height,
                        nbands = 1,
                        dtype = eltype(array_t),
                        options = options) do dataset
            band = ArchGDAL.getband(dataset, 1)
            ArchGDAL.write!(band, array_t)
            ArchGDAL.setnodatavalue!(band, -9999.0)
            ArchGDAL.setgeotransform!(dataset, transform)
            ArchGDAL.setproj!(dataset, wkt)
            ArchGDAL.write(dataset, fn,
                           driver = ArchGDAL.getdriver("GTiff"),
                           options = options)
        end
    end
    nothing
end

end # module CircuitscapeArchGDALExt
