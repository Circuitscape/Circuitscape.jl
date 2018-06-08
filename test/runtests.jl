using Circuitscape
using Base.Test
import Circuitscape: compute_single, compute_cholmod, compute_parallel,
                     runtests

# Utility to compare output files
include("compare_output.jl")

# Unit tests for internals
@testset "Unit tests" begin 
    include("internal.jl")
end

for f in (:compute, :compute_single, :compute_cholmod, :compute_parallel)
    runtests(f)
end
