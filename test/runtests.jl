using CircuitScape
using Base.Test

include("internal.jl")
include("compare_output.jl")

# Network pairwise tests
for i = 1:3
    r = compute("input/network/sgNetworkVerify$(i).ini")
    x = readdlm("output_verify/sgNetworkVerify$(i)_resistances.out")
    x = x[2:end, 2:end]
    @test sumabs2(x - r) < 1e-6
    compare_all_output("sgNetworkVerify$(i)")
end

# Network advanced tests
for i = 1:3
    r = compute("input/network/mgNetworkVerify$(i).ini")
    x = readdlm("output_verify/mgNetworkVerify$(i)_voltages.txt")
    x = x[:,2]
    @test sumabs2(x - r) < 1e-6
    compare_all_output("mgNetworkVerify$(i)")
end


# Raster pairwise tests
for i in deleteat!(collect(1:15), 12)
    r = compute("input/raster/pairwise/$i/sgVerify$(i).ini")
    x = readdlm("output_verify/sgVerify$(i)_resistances.out")
    x = x[2:end, 2:end]
    @test sumabs2(x - r) < 1e-6
    compare_all_output("sgVerify$(i)")
end

# Raster advanced tests
for i in 1:5
    r = compute("input/raster/advanced/$i/mgVerify$(i).ini")
    x = readdlm("output_verify/mgVerify$(i)_voltmap.asc"; skipstart = 6)
    @test sumabs2(x - r) < 1e-6
end

# Raster one to all test
for i in 1:13
    r = compute("input/raster/one_to_all/$i/oneToAllVerify$(i).ini")
    x = readdlm("output_verify/oneToAllVerify$(i)_resistances.out")
    x = x[:,2]

    @test sumabs2(x - r) < 1e-6
end

# Raster all to one test
for i in 1:12
    r = compute("input/raster/all_to_one/$i/allToOneVerify$(i).ini")
    x = readdlm("output_verify/allToOneVerify$(i)_resistances.out")
    x = x[:,2]

    @test sumabs2(x - r) < 1e-6
end
