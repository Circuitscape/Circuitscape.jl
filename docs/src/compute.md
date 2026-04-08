# Computational Limitations, Speed, and Landscape Size

We have tested this code on landscapes with up to 437 million cells. Increasing numbers of connections using diagonal (eight neighbor) connections will decrease the size of landscapes that can be analyzed. Also, increasing landscape size or numbers of focal nodes will increase computation time. Note that due to the matrix algebra involved with solving many pairs of focal nodes, Circuitscape will run much faster when focal points (each focal node falls within only one grid cell), rather than focal regions (at least one focal node occupies multiple grid cells), are used.

## Memory Limitations

There are several ways to increase the solvable grid size:

- Set impermeable areas of your resistance map to NODATA
- Use focal points instead of regions in pairwise mode
- Connect cells to their four neighbors only (`connect_four_neighbors_only = True`)
- Disable current and voltage maps (`write_cur_maps = False`, `write_volt_maps = False`)
- Use the one-to-all or all-to-one modes, which typically use less memory and run more quickly than pairwise mode
- Use the `cg+amg` solver instead of `cholmod`, `accelerate`, or `pardiso` for large grids (direct solvers use significantly more memory)
- Coarsen your grids (use larger cell sizes) -- this often produces qualitatively similar results (see McRae et al. 2008)

The all-to-one mode can be an alternative to pairwise mode when the goal is to produce a cumulative map of important connectivity areas among multiple source/target patches.

## Multi-threading

Circuitscape uses Julia's native threading for parallel computation. Start Julia
with multiple threads to take advantage of this:

```bash
julia -t 4    # use 4 threads
```

The `cg+amg` solver benefits most from threading — each focal point pair is solved
independently on a separate thread. For problems with many focal points, this can
provide significant speedups.

The `cholmod`, `accelerate`, and `pardiso` solvers perform batched direct solves.
Threading in these modes parallelizes the postprocessing step (current map
accumulation and output writing).

## Default Solvers: CG+AMG and CHOLMOD

Circuitscape ships with two solvers that work out of the box with no additional
packages:

- **CG+AMG** (`solver = cg+amg`, the default) — an iterative solver using
  conjugate gradient with an algebraic multigrid preconditioner. This is the
  best choice for large grids because memory usage scales well with problem
  size. It also parallelizes individual pair solves across threads.

- **CHOLMOD** (`solver = cholmod`) — a direct solver using Cholesky
  factorization from SuiteSparse. It can be significantly faster than CG+AMG
  for small to medium problems, but memory usage grows quickly due to fill-in,
  making it impractical for very large grids. In pairwise mode it performs
  batched solves controlled by the `cholmod_batch_size` parameter.

## Pardiso Solver

Circuitscape supports the [Pardiso](https://github.com/JuliaSparse/Pardiso.jl)
direct solver as a package extension. Pardiso uses Intel MKL's sparse direct
solver and can offer excellent performance on Intel hardware. To use it, first
install Pardiso.jl:

```julia
using Pkg
Pkg.add("Pardiso")
```

Then load it before (or alongside) Circuitscape:

```julia
using Pardiso
using Circuitscape
compute("config.ini")  # with solver = pardiso in the INI file
```

Pardiso requires double precision and will automatically switch if single
precision is requested. Like CHOLMOD, it is a direct solver best suited for
small to medium problem sizes, and uses the same batched solve strategy in
pairwise mode.

## Apple Accelerate Solver

On macOS (13.4 or later), Circuitscape can use Apple's
[Accelerate](https://developer.apple.com/documentation/accelerate/sparse_solvers)
framework for sparse Cholesky factorization. This is a direct solver that can
provide high performance on Apple Silicon hardware. To use it:

```julia
using AppleAccelerate
using Circuitscape
compute("config.ini")  # with solver = accelerate in the INI file
```

Or install `AppleAccelerate` first:

```julia
using Pkg
Pkg.add("AppleAccelerate")
```

The Accelerate solver supports both single and double precision, and uses the
same batched solve strategy as the CHOLMOD and Pardiso solvers.

!!! note
    Circuitscape sets BLAS to single-threaded at startup. Its workload is
    predominantly sparse matrix operations which do not benefit from
    multi-threaded BLAS.
