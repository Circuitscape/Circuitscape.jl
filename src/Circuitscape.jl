module Circuitscape
using AlgebraicMultigrid
using Krylov
using GZip

using LinearAlgebra
using SparseArrays
using DelimitedFiles
using Logging
using Dates
using SuiteSparse
using TimerOutputs
using Base.ScopedValues

# One timer per compute() call. A scoped value rather than a global so that
# concurrent jobs (e.g. Omniscape solving windows in parallel) do not reset or
# merge into each other's timings; tasks spawned inside inherit it.
const CSTIMER = ScopedValue(TimerOutput())

function __init__()
    LinearAlgebra.BLAS.set_num_threads(1)
end

include("consts.jl")
include("config.jl")
include("logging.jl")
include("utils.jl")
include("io.jl")
include("geometry.jl")
include("core.jl")
include("out.jl")
include("network/pairwise.jl")
include("raster/pairwise.jl")
include("raster/advanced.jl")
include("network/advanced.jl")
include("raster/onetoall.jl")
include("run.jl")

"""
    Circuitscape.start()

Build an INI file interactively, step by step, and optionally run it. The
builder lives in an extension on the REPL stdlib; calling `start()` loads
REPL on demand, so it works from a script as well as from the prompt.
"""
function start(args...)
    isempty(args) || throw(MethodError(start, args))
    load_repl_extension()
    Base.invokelatest(start)
end

const REPL_PKGID = Base.PkgId(Base.UUID("3fa0cd96-eef1-5676-8a61-b3b8758bbffb"), "REPL")

# Load the REPL stdlib if it is not already, which activates
# CircuitscapeREPLExt and with it the zero-argument start() method.
function load_repl_extension()
    ext = Base.get_extension(@__MODULE__, :CircuitscapeREPLExt)
    ext === nothing || return ext
    Base.require(REPL_PKGID)
    ext = Base.get_extension(@__MODULE__, :CircuitscapeREPLExt)
    ext === nothing && error("Loading the REPL stdlib did not activate the INI builder extension")
    ext
end
end
