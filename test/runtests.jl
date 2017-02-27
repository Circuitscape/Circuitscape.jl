using CircuitScape
using Base.Test

# First test does not work because of the solver
# TODO: Replace IterativeSolvers by PyAMG for now
#r = CircuitScape.network("sgNetworkVerify2_graph_conductances.txt", "sgNetworkVerify2_focal_nodes.txt")

#a = readdlm("sgNetworkVerify2_resistances.out")
#a = a[2:end, 2:end]

#@test a ≈ r

# Network test with multiple connected components
r = compute("sgNetworkVerify1.ini")

x = readdlm("sgNetworkVerify1_resistances.out")
x = x[2:end, 2:end]

@test r ≈ x
