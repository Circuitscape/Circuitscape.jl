module CircuitscapeAppleAccelerateExt

using AppleAccelerate
using SparseArrays
using LinearAlgebra
import Circuitscape: AccelerateSolver, construct_cholesky_factor, solve_linear_system,
                     refine_columns!, regularize

# AppleAccelerate.jl only accepts `Int64`-indexed matrices, though it copies
# the indices into Accelerate's own 32-bit row index / 64-bit column start
# arrays either way, so a 32-bit matrix is converted on the way in; the copy
# is transient.
construct_cholesky_factor(matrix::SparseMatrixCSC{T}, ::AccelerateSolver) where T =
    AppleAccelerate.AAFactorization(convert(SparseMatrixCSC{T,Int64}, regularize(matrix)))

function solve_linear_system(factor::AppleAccelerate.AAFactorization, matrix, rhs; tol = 1e-4)
    refine_columns!(factor \ rhs, factor, matrix, rhs, tol, "Apple Accelerate")
end

end # module CircuitscapeAppleAccelerateExt
