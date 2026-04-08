module CircuitscapeAppleAccelerateExt

using AppleAccelerate
using SparseArrays
using LinearAlgebra
import Circuitscape: AccelerateSolver, construct_cholesky_factor, solve_linear_system

function construct_cholesky_factor(matrix, ::AccelerateSolver)
    T = eltype(matrix)
    regularized = matrix + sparse(T(10) * eps(T) * I, size(matrix)...)
    factor = AppleAccelerate.AAFactorization(regularized)
    factor
end

function solve_linear_system(factor::AppleAccelerate.AAFactorization, matrix, rhs)
    lhs = factor \ rhs
    for col = 1:size(rhs, 2)
        residual = norm(matrix * lhs[:, col] .- rhs[:, col]) / norm(rhs[:, col])
        residual < 1e-4 || error("Apple Accelerate solver residual $residual exceeds tolerance 1e-4 for column $col")
    end
    lhs
end

end # module CircuitscapeAppleAccelerateExt
