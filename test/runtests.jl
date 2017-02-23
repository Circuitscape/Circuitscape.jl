using CircuitScape
using Base.Test

r = CircuitScape.network("sgNetworkVerify2_graph_conductances.txt", "sgNetworkVerify2_focal_nodes.txt")

a = readdlm("sgNetworkVerify2_resistances.out")
a = a[2:end, 2:end]

@test a ≈ r
