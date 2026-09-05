# Circuitscape

[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://docs.circuitscape.org/Circuitscape.jl/latest/)
[![Build Status](https://github.com/Circuitscape/Circuitscape.jl/workflows/CI/badge.svg)](https://github.com/Circuitscape/Circuitscape.jl/actions?query=workflow%3ACI)
[![codecov.io](https://codecov.io/github/Circuitscape/Circuitscape.jl/coverage.svg?branch=master)](https://codecov.io/github/Circuitscape/Circuitscape.jl?branch=master)

Circuitscape is an open-source program that uses circuit theory to model connectivity
in heterogeneous landscapes. Its most common applications include modeling movement and gene flow
of plants and animals, as well as identifying areas important for connectivity conservation.

Circuitscape is written in [Julia](https://julialang.org) for high performance and scalability.
More detail about the implementation can be found in the
[JuliaCon paper](https://proceedings.juliacon.org/papers/10.21105/jcon.00058).

> [!NOTE]
> Circuitscape.jl requires [Julia v1.11](https://julialang.org/downloads/) or later.

## Installation

1. [Install Julia](https://julialang.org/downloads/).

2. At the Julia prompt, install Circuitscape:

```julia
julia> using Pkg
julia> Pkg.add("Circuitscape")
```

To install the latest development version:

```julia
julia> Pkg.add(PackageSpec(name="Circuitscape", rev="master"))
```

Run the test suite with:

```julia
julia> Pkg.test("Circuitscape")
```

### Optional: GeoTIFF support

ESRI ASCII grids (`.asc`) work out of the box. To read or write GeoTIFF
rasters, add [ArchGDAL](https://github.com/yeesian/ArchGDAL.jl) and load it
before Circuitscape:

```julia
julia> Pkg.add("ArchGDAL")
julia> using ArchGDAL, Circuitscape
```

## Usage

Circuitscape jobs are configured via INI files. See the
[documentation](https://docs.circuitscape.org/Circuitscape.jl/latest/) for a full guide
on data types, calculation modes, and all available options.

```julia
julia> using Circuitscape
julia> Circuitscape.compute("path/to/config/file.ini")
```

Nothing is exported; every function is called as `Circuitscape.name`
(or use `using Circuitscape: compute` to bring the bare name in).

You can also build INI files interactively using the built-in terminal UI:

```julia
julia> using Circuitscape
julia> Circuitscape.start()
```

Or construct a configuration programmatically:

```julia
julia> using Circuitscape
julia> cfg = Circuitscape.init_config()
julia> cfg["habitat_file"] = "resistance_map.asc"
julia> cfg["point_file"] = "focal_nodes.asc"
julia> cfg["scenario"] = "pairwise"
julia> cfg["output_file"] = "output/results.out"
julia> Circuitscape.compute(cfg)
```

Example INI files can be found in the
[test folder](https://github.com/Circuitscape/Circuitscape.jl/tree/master/test/input).

## Features

### Solver Modes

Circuitscape supports four solver modes:

- **CG+AMG** (default) — an iterative solver using algebraic multigrid preconditioning. Scales well to large problems.
- **CHOLMOD** — a direct solver using [Cholesky decomposition](https://en.wikipedia.org/wiki/Cholesky_decomposition) via the [SuiteSparse](https://people.engr.tamu.edu/davis/suitesparse.html) library. Can be significantly faster for smaller problems, but memory use grows quickly with problem size due to [fill-in](https://algowiki-project.org/en/Cholesky_method#Reordering_to_reduce_the_number_of_fill-in_elements).
- **Apple Accelerate** — a direct solver using Apple's [Accelerate](https://developer.apple.com/documentation/accelerate/sparse_solvers) framework sparse Cholesky, available as a package extension on macOS (requires macOS 13.4+). Requires `using AppleAccelerate` before `using Circuitscape`.
- **Pardiso** — a direct solver via [Pardiso.jl](https://github.com/JuliaSparse/Pardiso.jl), available as a package extension. Requires `using Pardiso` before `using Circuitscape`.

To select a solver, add to your INI file:
```
solver = cholmod   # or accelerate, pardiso
```

### Parallel Computing

Circuitscape supports multi-threaded computation on Linux, macOS, and Windows.
Start Julia with `julia -t N` (or set `JULIA_NUM_THREADS`) to use N threads **and**
set `parallelize = True` in your INI file; both are needed.

### Single Precision

Circuitscape can run in single precision for reduced memory use at the cost of solution
accuracy. It is tested for the CG+AMG, CHOLMOD and Apple Accelerate solvers (Pardiso is
double-only). Add to your INI file:
```
precision = single
```

### Configuration checking

Since 6.0, misspelt values such as `solver = cholmodd` (for `cholmod`) or `write_cur_maps = Ture` (for `True`), missing
input files and a missing output directory are reported as errors before any data is
read, instead of silently running something else. See
[Upgrading to 6.0](https://docs.circuitscape.org/Circuitscape.jl/latest/migration/)
for what to change when upgrading.

## Performance

In a 2021 benchmark, Circuitscape.jl 5.x was up to **4x faster** than the legacy Python
version (Circuitscape v4.0.5) on 16 parallel processes, with the CHOLMOD solver providing
the best performance on suitable problem sizes. The benchmark predates the 6.0
performance work (threaded pair scheduling, native `.asc` I/O, lower memory use), which
is documented with measurements in the
[solvers and computation time](https://docs.circuitscape.org/Circuitscape.jl/latest/compute/)
page.

<img src="https://raw.githubusercontent.com/Circuitscape/Circuitscape.jl/master/docs/src/benchmark/benchmark.png" width=650 height=450>

Benchmarks were run on a Linux (Ubuntu) server with an Intel Xeon Silver 4114 CPU
(2.20 GHz, 20 cores, 384 GB RAM) using problems from the
[benchmark suite](https://github.com/Circuitscape/BigTests).

## Related Projects

- [Omniscape.jl](https://github.com/Circuitscape/Omniscape.jl) — Omnidirectional connectivity analysis built on Circuitscape
- [AlgebraicMultigrid.jl](https://github.com/JuliaLinearAlgebra/AlgebraicMultigrid.jl) — The default iterative solver used by Circuitscape

## Contributing

If you encounter any issues or would like to ask a question, please file
a report [here](https://github.com/Circuitscape/Circuitscape.jl/issues).
Contributions in the form of
[pull requests](https://github.com/Circuitscape/Circuitscape.jl/pulls) are welcome!
