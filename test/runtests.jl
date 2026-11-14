using Circuitscape
using Test
using Logging

Logging.disable_logging(Logging.Warn)

# GDAL is optional: without ArchGDAL the GeoTIFF fixtures are skipped.
const GDAL_AVAILABLE = try
    @eval using ArchGDAL
    true
catch
    false
end
GDAL_AVAILABLE || println("ArchGDAL not available: skipping GeoTIFF fixtures")

include("test_utils.jl")

clean_output()

@testset "Unit tests" begin
    include("internal.jl")
    include("read_point_map.jl")
end

runtests(solver="cg+amg", parallel=true, precision="double")
runtests(solver="cg+amg", parallel=true, precision="single")
runtests(solver="cholmod", parallel=true, precision="double")
runtests(solver="cholmod", parallel=true, precision="single")

accelerate_available = try
    @eval using AppleAccelerate
    true
catch
    false
end
if accelerate_available
    runtests(solver="accelerate", parallel=true, precision="double")
    runtests(solver="accelerate", parallel=true, precision="single")
else
    println("Skipping Apple Accelerate tests (not available)")
end

pardiso_available = try
    @eval using Pardiso
    Pardiso.MKLPardisoSolver()
    true
catch
    false
end
if pardiso_available
    runtests(solver="pardiso", parallel=true, precision="double")
else
    println("Skipping Pardiso tests (MKL not available)")
end

# AMGX needs AMGX.jl in the environment and a working CUDA GPU.
amgx_available = try
    @eval using AMGX
    @eval AMGX.CUDA.functional()
catch
    false
end
if amgx_available
    runtests(solver="amgx", parallel=true, precision="double")
    runtests(solver="amgx", parallel=true, precision="single")
else
    println("Skipping AMGX tests (AMGX.jl or a CUDA GPU not available)")
end

@testset "Issue 341: included pairs" begin
    include("issue341.jl")
end

clean_output()
