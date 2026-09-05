module CircuitscapePardisoExt

using Pardiso
using SparseArrays
using LinearAlgebra
import Circuitscape: PardisoSolver, construct_cholesky_factor, solve_linear_system, refine_columns!

mutable struct PardisoFactorize
    const ps::Pardiso.MKLPardisoSolver
    const verbose::Bool
    firsttime::Bool
end
PardisoFactorize(;verbose=false) = PardisoFactorize(Pardiso.MKLPardisoSolver(), verbose, true)

function (p::PardisoFactorize)(x, A, b, update_matrix=false; kwargs...)
    if p.firsttime
        Pardiso.set_phase!(p.ps, Pardiso.ANALYSIS_NUM_FACT)
        Pardiso.pardiso(p.ps, x, A, b)
        p.firsttime = false
    end

    if update_matrix
        Pardiso.set_phase!(p.ps, Pardiso.NUM_FACT)
        Pardiso.pardiso(p.ps, x, A, b)
    end

    Pardiso.set_phase!(p.ps, Pardiso.SOLVE_ITERATIVE_REFINE)
    Pardiso.pardiso(p.ps, x, A, b)
end

construct_cholesky_factor(matrix, ::PardisoSolver) =
    PardisoFactorize()

# Lets refine_columns! use `\` on the Pardiso handle for a fixed matrix.
struct PardisoOp
    p::PardisoFactorize
    A::SparseMatrixCSC
end
function Base.:\(o::PardisoOp, r::AbstractVector)
    x = zeros(eltype(r), length(r))
    o.p(x, o.A, Vector(r))
    x
end

function solve_linear_system(factor::PardisoFactorize, matrix, rhs; tol = 1e-4)
    mat = sparse(10eps(eltype(matrix)) * I, size(matrix)...) + matrix
    op = PardisoOp(factor, mat)
    lhs = similar(rhs)
    for i = 1:size(lhs, 2)
        lhs[:, i] .= op \ rhs[:, i]
    end
    # Pardiso already runs SOLVE_ITERATIVE_REFINE; this is the same
    # residual-driven refinement the CHOLMOD and Accelerate backends apply,
    # so all three direct solvers are held to the same standard.
    refine_columns!(lhs, op, mat, rhs, tol, "Pardiso")
end

end # module CircuitscapePardisoExt
