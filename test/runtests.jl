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
