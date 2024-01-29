module CircuitscapeMKLPardiso

using Pardiso
import Circuitscape: Solver, _check_eltype, MKLPardisoSolver, construct_cholesky_factor, solve_linear_system

mutable struct MKLPardisoFactorize
    A
    ps::Pardiso.MKLPardisoSolver
    verbose::Bool
    firsttime::Bool
end
MKLPardisoFactorize(;verbose=false) = MKLPardisoFactorize(nothing,Pardiso.MKLPardisoSolver(),verbose,true)
function (p::MKLPardisoFactorize)(x,A,b,update_matrix=false;kwargs...)
    if p.firsttime
        Pardiso.set_phase!(p.ps, Pardiso.ANALYSIS_NUM_FACT)
        Pardiso.pardiso(p.ps, x, A, b)
        p.firsttime = false
    end

    if update_matrix
        Pardiso.set_phase!(p.ps, Pardiso.NUM_FACT)
        Pardiso.pardiso(p.ps, x, A, b)
        p.A = A
    end

    Pardiso.set_phase!(p.ps, Pardiso.SOLVE_ITERATIVE_REFINE)
    Pardiso.pardiso(p.ps, x, A, b)
end

function compute_mklpardiso(str, batch_size = 5)
    cfg = parse_config(str)
    T = cfg["precision"] in SINGLE ? Float32 : Float64
    if T == Float32
        cswarn("Pardiso supports only double precision. Changing precision to double.")
        T = Float64
    end
    V = cfg["use_64bit_indexing"] in TRUELIST ? Int64 : Int32
    cfg["solver"] = "mklpardiso"
    _compute(T, V, cfg)
end


_check_eltype(a, solver::MKLPardisoSolver) = a
construct_cholesky_factor(matrix, ::MKLPardisoSolver, suppress_info::Bool) =
            MKLPardisoFactorize()

function solve_linear_system(factor::MKLPardisoFactorize, matrix, rhs)
    lhs = similar(rhs)
	mat = sparse(10eps(eltype(matrix))*I,size(matrix)...) + matrix
    x = zeros(eltype(matrix), size(matrix, 1))
    for i = 1:size(lhs, 2)
        factor(x, mat, rhs[:,i])
		@assert (norm(mat*x .- rhs[:,i]) / norm(rhs[:,i])) < 1e-4
        lhs[:,i] .= x
    end
    lhs
end

end # module CircuitscapeMKLPardiso
