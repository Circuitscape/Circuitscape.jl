using Circuitscape
using Base.Test
import Circuitscape: compute_single, compute_cholmod

# Utility to compare output files
include("compare_output.jl")

# Unit tests for internals
@testset "Unit tests" begin 
include("internal.jl")
end


function f()
for f in (compute, compute_single, compute_cholmod)
str =  f == compute ? "Double" : "Single"

@testset "$str Precision Tests" begin 

@testset "Network Pairwise" begin 
# Network pairwise tests
for i = 1:3
    info("Testing sgNetworkVerify$i")
    r = f("input/network/sgNetworkVerify$(i).ini")
    calculate_cum_current_maps("output/sgNetworkVerify$(i).out")
    x = readdlm("output_verify/sgNetworkVerify$(i)_resistances.out")
    valx = x[2:end, 2:end]
    valr = r[2:end, 2:end]
    @test sum(abs2, valx - valr) < 1e-6
    pts_x = x[2:end,1]
    pts_r = r[2:end,1]
    @test pts_x + 1 == pts_r
    compare_all_output("sgNetworkVerify$(i)")
    info("Test sgNetworkVerify$i passed")
end
end

@testset "Network Advanced" begin
# Network advanced tests
for i = 1:3
    info("Testing mgNetworkVerify$i")
    r = f("input/network/mgNetworkVerify$(i).ini")
    calculate_cum_current_maps("output/mgNetworkVerify$(i).out")
    x = readdlm("output_verify/mgNetworkVerify$(i)_voltages.txt")
    x[:,1] = x[:,1] + 1
    @test sum(abs2, x - r) < 1e-6
    compare_all_output("mgNetworkVerify$(i)")
    info("Test mgNetworkVerify$i passed")
end
end


@testset "Raster Pairwise" begin 
# Raster pairwise tests
for i = 1:15
    info("Testing sgVerify$i")
    r = f("input/raster/pairwise/$i/sgVerify$(i).ini")
    calculate_cum_current_maps("output/sgVerify$(i).out")
    x = readdlm("output_verify/sgVerify$(i)_resistances.out")
    # x = x[2:end, 2:end]
    @test sum(abs2, x - r) < 1e-6
    compare_all_output("sgVerify$(i)")
    info("Test sgVerify$i passed")
end
end

@testset "Raster Advanced" begin 
# Raster advanced tests
for i in 1:5
    info("Testing mgVerify$i")
    r = f("input/raster/advanced/$i/mgVerify$(i).ini")
    calculate_cum_current_maps("output/mgVerify$(i).out")
    x = readdlm("output_verify/mgVerify$(i)_voltmap.asc"; skipstart = 6)
    @test sum(abs2, x - r) < 1e-5
    # compare_all_output("mgVerify$(i)")
    info("Test mgVerify$i passed")
end
end

@testset "Raster One to All" begin 
# Raster one to all test
for i in 1:13
    info("Testing oneToAllVerify$i")
    r = f("input/raster/one_to_all/$i/oneToAllVerify$(i).ini")
    calculate_cum_current_maps("output/oneToAllVerify$(i).out")
    x = readdlm("output_verify/oneToAllVerify$(i)_resistances.out")
    # x = x[:,2]
    @test sum(abs2, x - r) < 1e-6
    compare_all_output("oneToAllVerify$(i)")
    info("Test oneToAllVerify$i passed")
end
end

@testset "Raster ALl to One" begin
# Raster all to one test
for i in 1:12
    info("Testing allToOneVerify$i")
    r = f("input/raster/all_to_one/$i/allToOneVerify$(i).ini")
    calculate_cum_current_maps("output/allToOneVerify$(i).out")
    x = readdlm("output_verify/allToOneVerify$(i)_resistances.out")
    # x = x[:,2]

    @test sum(abs2, x - r) < 1e-6
    info("Test allToOneVerify$i passed")
end
end

end
end
end

f()

@testset "Parallel Tests" begin 
addprocs(2)

@everywhere using Circuitscape
@everywhere import Circuitscape: compute_cholmod, compute_single
f()
end
