using Circuitscape
using Test
using Logging

include("test_utils.jl")

# Unit tests for internals
@testset "Unit tests" begin
    include("internal.jl")
end

for f in (compute, compute_cholmod, compute_parallel)
    runtests(f)
end

# Run Pardiso tests only if Pardiso.jl is available and MKL is present
pardiso_available = try
    @eval using Pardiso
    Pardiso.MKLPardisoSolver()
    true
catch
    false
end
if pardiso_available
    runtests(compute_pardiso)
end
