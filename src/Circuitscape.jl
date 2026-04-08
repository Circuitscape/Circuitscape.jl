module Circuitscape
using AlgebraicMultigrid
using ArchGDAL
using Graphs
using Krylov
using GZip

using LinearAlgebra
using SparseArrays
using DelimitedFiles
using Logging
using Dates
using SuiteSparse
using TimerOutputs

const CSTIMER = TimerOutput()

function __init__()
    LinearAlgebra.BLAS.set_num_threads(1)
end

include("consts.jl")
include("config.jl")
include("logging.jl")
include("utils.jl")
include("io.jl")
include("out.jl")
include("core.jl")
include("network/pairwise.jl")
include("raster/pairwise.jl")
include("raster/advanced.jl")
include("network/advanced.jl")
include("raster/onetoall.jl")
include("run.jl")
include("INIBuilder/INIBuilder.jl")

using .INIBuilder
export start

end
