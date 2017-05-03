using CircuitScape
using Base.Test

include("internal.jl")

# Simple test with one connected component
r = compute("input/network/sgNetworkVerify2.ini")
x = readdlm("output_verify/sgNetworkVerify2_resistances.out")
x = x[2:end, 2:end]

@test sumabs2(x - r) < 1e-6

# Network test with multiple connected components
r = compute("input/network/sgNetworkVerify1.ini")
x = readdlm("output_verify/sgNetworkVerify1_resistances.out")
x = x[2:end, 2:end]

@test sumabs2(x - r) < 1e-6

# Network test with multiple cc and resistance specified
r = compute("input/network/sgNetworkVerify3.ini")
x = readdlm("output_verify/sgNetworkVerify3_resistances.out")
x = x[2:end, 2:end]

@test sumabs2(x - r) < 1e-6

# Simple Network test with advanced mode
r = compute("input/network/mgNetworkVerify1.ini")
x = readdlm("output_verify/mgNetworkVerify1_voltages.txt")
x = x[:,2]

@test sumabs2(x - r) < 1e-6 

# Network test with advanced mode and multiple grounds
r = compute("input/network/mgNetworkVerify2.ini")
x = readdlm("output_verify/mgNetworkVerify2_voltages.txt")
x = x[:,2]

@test sumabs2(x - r) < 1e-6 

# Network test with advanced mode and both multiple grounds and sources
r = compute("input/network/mgNetworkVerify3.ini")
x = readdlm("output_verify/mgNetworkVerify3_voltages.txt")
x = x[:,2]

@test sumabs2(x - r) < 1e-6 

# Simple Raster test
r = compute("input/raster/pairwise/2/sgVerify2.ini")
x = readdlm("output_verify/sgVerify2_resistances.out")
x = x[2:end, 2:end]

@test sumabs2(x - r) < 1e-6

# Another raster test, sort points map
r = compute("input/raster/pairwise/1/sgVerify1.ini")
x = readdlm("output_verify/sgVerify1_resistances.out")
x = x[2:end, 2:end]

@test sumabs2(x - r) < 1e-6

# Raster test, points in txt file
r = compute("input/raster/pairwise/7/sgVerify7.ini")
x = readdlm("output_verify/sgVerify7_resistances.out")
x = x[2:end, 2:end]

@test sumabs2(x - r) < 1e-6

# Raster test, no polygon
r = compute("input/raster/pairwise/4/sgVerify4.ini")
x = readdlm("output_verify/sgVerify4_resistances.out")
x = x[2:end, 2:end]

@test sumabs2(x - r) < 1e-6

# Raster test, no polygon
r = compute("input/raster/pairwise/15/sgVerify15.ini")
x = readdlm("output_verify/sgVerify15_resistances.out")
x = x[2:end, 2:end]

@test sumabs2(x - r) < 1e-6

# Another Raster test, no polygon
r = compute("input/raster/pairwise/14/sgVerify14.ini")
x = readdlm("output_verify/sgVerify14_resistances.out")
x = x[2:end, 2:end]

@test sumabs2(x - r) < 1e-6

# Raster test, with polygons in focal nodes
r = compute("input/raster/pairwise/3/sgVerify3.ini")
x = readdlm("output_verify/sgVerify3_resistances.out")
x = x[2:end, 2:end]

@test sumabs2(x - r) < 1e-6

# Raster test, with polygons in focal nodes
r = compute("input/raster/pairwise/5/sgVerify5.ini")
x = readdlm("output_verify/sgVerify5_resistances.out")
x = x[2:end, 2:end]

@test sumabs2(x - r) < 1e-6

# Raster test, with polygons in focal nodes
r = compute("input/raster/pairwise/8/sgVerify8.ini")
x = readdlm("output_verify/sgVerify8_resistances.out")
x = x[2:end, 2:end]

@test sumabs2(x - r) < 1e-6

# Raster test, with polygons in focal nodes
r = compute("input/raster/pairwise/9/sgVerify9.ini")
x = readdlm("output_verify/sgVerify9_resistances.out")
x = x[2:end, 2:end]

@test sumabs2(x - r) < 1e-6

# Raster test, with polygons in focal nodes, but no polygons file
r = compute("input/raster/pairwise/6/sgVerify6.ini")
x = readdlm("output_verify/sgVerify6_resistances.out")
x = x[2:end, 2:end]

@test sumabs2(x - r) < 1e-6

# Raster test, with polygons in focal nodes, no polygons file
r = compute("input/raster/pairwise/10/sgVerify10.ini")
x = readdlm("output_verify/sgVerify10_resistances.out")
x = x[2:end, 2:end]

@test sumabs2(x - r) < 1e-6

# Raster test, with polygons in focal nodes, no polygons file
r = compute("input/raster/pairwise/11/sgVerify11.ini")
x = readdlm("output_verify/sgVerify11_resistances.out")
x = x[2:end, 2:end]

@test sumabs2(x - r) < 1e-6

# Raster test, custom config file
r = compute("input/raster/pairwise/13/sgVerify13.ini")
x = readdlm("output_verify/sgVerify13_resistances.out")
x = x[2:end, 2:end]

@test sumabs2(x - r) < 1e-6

# Raster advanced test
r = compute("input/raster/advanced/1/mgVerify1.ini")
x = readdlm("output_verify/mgVerify1_voltmap.asc", skipstart = 6)

@test sumabs2(x - r) < 1e-6

# Another Raster advanced test
r = compute("input/raster/advanced/3/mgVerify3.ini")
x = readdlm("output_verify/mgVerify3_voltmap.asc", skipstart = 6)

@test sumabs2(x - r) < 1e-6

# Another Raster advanced test
r = compute("input/raster/advanced/4/mgVerify4.ini")
x = readdlm("output_verify/mgVerify4_voltmap.asc", skipstart = 6)

@test sumabs2(x - r) < 1e-6

# Another Raster advanced test
r = compute("input/raster/advanced/5/mgVerify5.ini")
x = readdlm("output_verify/mgVerify5_voltmap.asc", skipstart = 6)

@test sumabs2(x - r) < 1e-6

# Raster advanced test, multiple cc
r = compute("input/raster/advanced/2/mgVerify2.ini")
x = readdlm("output_verify/mgVerify2_voltmap.asc", skipstart = 6)

@test sumabs2(x - r) < 1e-6

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
