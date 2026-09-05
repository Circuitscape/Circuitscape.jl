module CircuitscapeAppleAccelerateExt

using AppleAccelerate
using SparseArrays
using LinearAlgebra
import Circuitscape: AccelerateSolver, construct_cholesky_factor, solve_linear_system, refine_columns!

function construct_cholesky_factor(matrix, ::AccelerateSolver)
    T = eltype(matrix)
    regularized = matrix + sparse(T(10) * eps(T) * I, size(matrix)...)
    factor = AppleAccelerate.AAFactorization(regularized)
    factor
end

function solve_linear_system(factor::AppleAccelerate.AAFactorization, matrix, rhs; tol = 1e-4)
    refine_columns!(factor \ rhs, factor, matrix, rhs, tol, "Apple Accelerate")
end

end # module CircuitscapeAppleAccelerateExt
