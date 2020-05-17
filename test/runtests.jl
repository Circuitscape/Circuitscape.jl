using Circuitscape
using Test
import Circuitscape: compute_single, compute_cholmod, compute_parallel,
                     runtests

filename = "input/raster/all_to_one/12/include_matrix.txt"

f = endswith(filename, "gz") ? GZip.open(filename, "r") : open(filename, "r")
filetype = Circuitscape._guess_file_type(filename, f)
close(f)

open(filename, "r") do f
    minval = parse(Float64, split(readline(f))[2])
    maxval = parse(Float64, split(readline(f))[2])
end

## Unit tests for internals
#@testset "Unit tests" begin
#    include("internal.jl")
#end

#for f in (compute, compute_single, compute_cholmod, compute_parallel)
#    runtests(f)
#end
