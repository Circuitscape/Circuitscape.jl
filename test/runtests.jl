using Circuitscape
using Test
using Logging

Logging.disable_logging(Logging.Warn)

include("test_utils.jl")

clean_output()

@testset "Unit tests" begin
    include("internal.jl")
end

runtests(solver="cg+amg", parallel=true)
runtests(solver="cholmod", parallel=true)

accelerate_available = try
    @eval using AppleAccelerate
    true
catch
    false
end
if accelerate_available
    runtests(solver="accelerate", parallel=true)
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
    runtests(solver="pardiso", parallel=true)
else
    println("Skipping Pardiso tests (MKL not available)")
end

@testset "Issue 341: included pairs" begin
    include("issue341.jl")
end

clean_output()
