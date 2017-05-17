using CircuitScape
using Base.Test

include("internal.jl")
include("compare_output.jl")

# Network pairwise tests
for i = 1:3
    info("Testing sgNetworkVerify$i")
    r = compute("input/network/sgNetworkVerify$(i).ini")
    x = readdlm("output_verify/sgNetworkVerify$(i)_resistances.out")
    x = x[2:end, 2:end]
    @test sumabs2(x - r) < 1e-6
    compare_all_output("sgNetworkVerify$(i)")
    info("Test sgNetworkVerify$i passed")
end

# Network advanced tests
for i = 1:3
    info("Testing mgNetworkVerify$i")
    r = compute("input/network/mgNetworkVerify$(i).ini")
    x = readdlm("output_verify/mgNetworkVerify$(i)_voltages.txt")
    x[:,1] = x[:,1] + 1
    @test sumabs2(x - r) < 1e-6
    compare_all_output("mgNetworkVerify$(i)")
    info("Test mgNetworkVerify$i passed")
end


# Raster pairwise tests
for i in deleteat!(collect(1:15), 12)
    info("Testing sgVerify$i")
    r = compute("input/raster/pairwise/$i/sgVerify$(i).ini")
    x = readdlm("output_verify/sgVerify$(i)_resistances.out")
    x = x[2:end, 2:end]
    @test sumabs2(x - r) < 1e-6
    compare_all_output("sgVerify$(i)")
    info("Test sgVerify$i passed")
end

# Raster advanced tests
for i in 1:5
    info("Testing mgVerify$i")
    r = compute("input/raster/advanced/$i/mgVerify$(i).ini")
    x = readdlm("output_verify/mgVerify$(i)_voltmap.asc"; skipstart = 6)
    @test sumabs2(x - r) < 1e-6
    compare_all_output("mgVerify$(i)")
    info("Test mgVerify$i passed")
end

# Raster one to all test
for i in 1:13
    info("Testing oneToAllVerify$i")
    r = compute("input/raster/one_to_all/$i/oneToAllVerify$(i).ini")
    x = readdlm("output_verify/oneToAllVerify$(i)_resistances.out")
    x = x[:,2]

    @test sumabs2(x - r) < 1e-6
    info("Test oneToAllVerify$i passed")
end

# Raster all to one test
for i in 1:12
    info("Testing allToOneVerify$i")
    r = compute("input/raster/all_to_one/$i/allToOneVerify$(i).ini")
    x = readdlm("output_verify/allToOneVerify$(i)_resistances.out")
    x = x[:,2]

    @test sumabs2(x - r) < 1e-6
    info("Test allToOneVerify$i passed")
end
