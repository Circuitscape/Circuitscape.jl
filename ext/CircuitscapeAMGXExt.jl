module CircuitscapeAMGXExt

using AMGX
using SparseArrays
using LinearAlgebra
import Circuitscape: AMGXSolver, construct_amgx_solver, release!, solve_linear_system,
                     solve_linear_system!, refine_columns!, regularize

# AMGX settings: CG preconditioned by one V-cycle of classical AMG (PMIS
# coarsening, distance-2 interpolation). On the BigTests grids it needs 7-8 CG
# iterations against 40-110 for AMGX's aggregation AMG, and solves 5-8x
# faster. CG needs a symmetric preconditioner, so the smoother is a symmetric
# Gauss-Seidel sweep with as many pre- as post-sweeps. AMGX stops on its own
# residual relative to the right-hand side; that is the true residual of the
# regularized matrix, and `refine_columns!` still checks every column against
# the user's tolerance on the unregularized one, as it does for the direct
# solvers.
amgx_tolerance(::Type{Float64}) = 1e-8
amgx_tolerance(::Type{Float32}) = 1e-6

amgx_config(::Type{T}) where T = join([
    "config_version=2",
    "solver(main)=PCG",
    "main:max_iters=1000",
    "main:convergence=RELATIVE_INI_CORE",
    "main:tolerance=$(amgx_tolerance(T))",
    "main:norm=L2",
    "main:monitor_residual=1",
    "main:print_solve_stats=0",
    "main:obtain_timings=0",
    "main:preconditioner(amg)=AMG",
    "amg:algorithm=CLASSICAL",
    "amg:selector=PMIS",
    "amg:interpolator=D2",
    "amg:cycle=V",
    "amg:max_iters=1",
    "amg:max_levels=50",
    "amg:presweeps=1",
    "amg:postsweeps=1",
    "amg:smoother(sm)=MULTICOLOR_GS",
    "sm:symmetric_GS=1",
    "sm:max_iters=1",
    "amg:coarse_solver(cs)=DENSE_LU_SOLVER",
    "amg:min_coarse_rows=2",
    "amg:dense_lu_num_rows=2048",
    "amg:print_grid_stats=0",
    "amg:print_solve_stats=0",
], ", ")

amgx_mode(::Type{Float64}) = AMGX.dDDI
amgx_mode(::Type{Float32}) = AMGX.dFFI

# AMGX is initialized once per process, on first use, and shares one set of
# resources (the CUDA device that is current then) between every solve. All
# calls into AMGX go through `LOCK`: one-to-all runs its points on several
# threads, and AMGX makes no promise about concurrent use of one resources
# object. The host-side work around the solves still runs in parallel.
const LOCK = ReentrantLock()
const RESOURCES = Ref{Any}(nothing)
const RESOURCES_CONFIG = Ref{Any}(nothing)

function amgx_resources()
    RESOURCES[] === nothing || return RESOURCES[]
    # AMGX prints its banner and statistics to stdout; route them to the
    # debug log instead.
    AMGX.register_print_callback(s -> (@debug("AMGX: " * rstrip(s)); nothing))
    AMGX.initialize()
    cfg = AMGX.Config("")
    RESOURCES_CONFIG[] = cfg
    RESOURCES[] = AMGX.Resources(cfg)
    atexit(amgx_shutdown)
    RESOURCES[]
end

function amgx_shutdown()
    lock(LOCK) do
        RESOURCES[] === nothing && return
        try
            close(RESOURCES[])
            close(RESOURCES_CONFIG[])
            AMGX.finalize()
        catch e
            @debug("AMGX shutdown failed", exception = e)
        end
        RESOURCES[] = nothing
        RESOURCES_CONFIG[] = nothing
    end
end

"""
    AMGXHandle

A component's matrix and multigrid hierarchy on the GPU, with the vectors
each solve uploads into and downloads from. Built by
`construct_amgx_solver`, freed by `release!`.
"""
mutable struct AMGXHandle{T}
    config::AMGX.Config
    matrix::AMGX.AMGXMatrix
    solver::AMGX.Solver
    rhs::AMGX.AMGXVector
    sol::AMGX.AMGXVector
    buf::Vector{T}
    open::Bool
end

function construct_amgx_solver(matrix::SparseMatrixCSC{T}, ::AMGXSolver) where T
    n = size(matrix, 1)
    A = regularize(matrix)
    nnz(A) < typemax(Cint) ||
        error("AMGX uses 32-bit indices and cannot hold a component with $(nnz(A)) " *
              "stored entries; use the cg+amg solver for this problem.")
    # AMGX takes CSR with zero-based `Cint` indices. The Laplacian is
    # symmetric, so its CSC arrays are also its CSR arrays.
    rowptr = Cint.(A.colptr) .- Cint(1)
    colind = Cint.(rowvals(A)) .- Cint(1)
    vals = nonzeros(A)
    mode = amgx_mode(T)
    lock(LOCK) do
        res = amgx_resources()
        config = AMGX.Config(amgx_config(T))
        created = Any[config]
        try
            m = AMGX.AMGXMatrix(res, mode); push!(created, m)
            AMGX.upload!(m, rowptr, colind, vals)
            s = AMGX.Solver(res, mode, config); push!(created, s)
            AMGX.setup!(s, m)
            rhs = AMGX.AMGXVector(res, mode); push!(created, rhs)
            sol = AMGX.AMGXVector(res, mode); push!(created, sol)
            AMGX.set_zero!(rhs, n)
            AMGX.set_zero!(sol, n)
            AMGXHandle{T}(config, m, s, rhs, sol, Vector{T}(undef, n), true)
        catch
            # Children before parents, as AMGX requires.
            foreach(close, reverse(created))
            rethrow()
        end
    end
end

function release!(h::AMGXHandle)
    lock(LOCK) do
        h.open || return
        foreach(close, (h.sol, h.rhs, h.solver, h.matrix, h.config))
        h.open = false
    end
    nothing
end

# One GPU solve of `h`'s regularized matrix: `x = A \ b`. A solve that stops
# short of AMGX's tolerance is not an error here; `refine_columns!` measures
# the true residual of the result and refines or rejects it.
function amgx_solve!(x::AbstractVector{T}, h::AMGXHandle{T}, b::AbstractVector) where T
    h.open || error("AMGX handle used after release!")
    copyto!(h.buf, b)
    lock(LOCK) do
        AMGX.upload!(h.rhs, h.buf)
        # Zero the guess explicitly: AMGX's solve with a zero initial guess
        # assumes, but does not ensure, that the solution vector holds zeros
        # (NVIDIA/AMGX#80), and on a reused vector returns its old contents
        # plus the new solution.
        AMGX.set_zero!(h.sol, length(h.buf))
        AMGX.solve!(h.sol, h.solver, h.rhs)
        status = AMGX.get_status(h.solver)
        status == AMGX.FAILED && error("AMGX solve failed")
        status == AMGX.SUCCESS ||
            @debug("AMGX solve stopped with status $status after " *
                   "$(AMGX.get_iterations_number(h.solver)) iterations")
        AMGX.download!(h.buf, h.sol)
    end
    copyto!(x, h.buf)
end

# `refine_columns!` corrects a column with `factor \ residual`.
Base.:\(h::AMGXHandle{T}, r::AbstractVector) where T =
    amgx_solve!(Vector{T}(undef, length(r)), h, r)

function solve_linear_system!(lhs, h::AMGXHandle, matrix, rhs;
                              tol = 1e-4, resid = similar(lhs))
    for i = 1:size(rhs, 2)
        amgx_solve!(view(lhs, :, i), h, view(rhs, :, i))
    end
    refine_columns!(lhs, h, matrix, rhs, tol, "AMGX"; resid)
end

solve_linear_system(h::AMGXHandle, matrix, rhs; tol = 1e-4) =
    solve_linear_system!(similar(rhs), h, matrix, rhs; tol)

end # module CircuitscapeAMGXExt
