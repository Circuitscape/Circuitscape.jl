using Circuitscape
using Test
using Pardiso
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

runtests(compute_pardiso)
