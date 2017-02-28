using CircuitScape
using Base.Test

# Simple test with one connected component
r = compute("sgNetworkVerify2.ini")
x = readdlm("sgNetworkVerify2_resistances.out")
x = x[2:end, 2:end]

@test sumabs2(x - r) < 1e-6

# Network test with multiple connected components
r = compute("sgNetworkVerify1.ini")
x = readdlm("sgNetworkVerify1_resistances.out")
x = x[2:end, 2:end]

@test sumabs2(x - r) < 1e-6

# Network test with multiple cc and resistance specified
r = compute("sgNetworkVerify3.ini")
x = readdlm("sgNetworkVerify3_resistances.out")
x = x[2:end, 2:end]

@test sumabs2(x - r) < 1e-6

# Network test with advanced mode
r = compute("mgNetworkVerify1.ini")
x = readdlm("mgNetworkVerify1_voltages.txt")
x = x[:,2]

@test sumabs2(x - r) < 1e-6 
