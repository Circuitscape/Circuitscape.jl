using Circuitscape
using Test
using Logging
Logging.disable_logging(Logging.Info)

include("test_utils.jl")

# Unit tests for internals
@testset "Unit tests" begin
    include("internal.jl")
end

for f in (compute,)
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
