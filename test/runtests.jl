using CircuitScape
using Base.Test

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
r = compute("input/raster/2/sgVerify2.ini")
x = readdlm("output_verify/sgVerify2_resistances.out")
x = x[2:end, 2:end]

@test sumabs2(x - r) < 1e-6
