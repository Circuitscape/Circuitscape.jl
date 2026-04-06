module CircuitscapePardisoExt

using Pardiso
using SparseArrays
using LinearAlgebra
import Circuitscape: PardisoSolver, construct_cholesky_factor, solve_linear_system

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

construct_cholesky_factor(matrix, ::PardisoSolver, suppress_info::Bool) =
    PardisoFactorize()

function solve_linear_system(factor::PardisoFactorize, matrix, rhs)
    lhs = similar(rhs)
    mat = sparse(10eps(eltype(matrix)) * I, size(matrix)...) + matrix
    x = zeros(eltype(matrix), size(matrix, 1))
    for i = 1:size(lhs, 2)
        factor(x, mat, rhs[:, i])
        @assert (norm(mat * x .- rhs[:, i]) / norm(rhs[:, i])) < 1e-4
        lhs[:, i] .= x
    end
    lhs
end

end # module CircuitscapePardisoExt
