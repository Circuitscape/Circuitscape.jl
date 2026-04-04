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

## Usage

Circuitscape jobs are configured via INI files. See the
[documentation](https://docs.circuitscape.org/Circuitscape.jl/latest/) for a full guide
on data types, calculation modes, and all available options.

```julia
julia> using Circuitscape
julia> compute("path/to/config/file.ini")
```

You can also build INI files interactively using the built-in terminal UI:

```julia
julia> using Circuitscape
julia> Circuitscape.INIBuilder.start()
```

Or construct a configuration programmatically:

```julia
julia> using Circuitscape
julia> cfg = Circuitscape.init_config()
julia> cfg["habitat_file"] = "resistance_map.asc"
julia> cfg["point_file"] = "focal_nodes.asc"
julia> cfg["scenario"] = "pairwise"
julia> cfg["output_file"] = "output/results.out"
julia> compute(cfg)
```

Example INI files can be found in the
[test folder](https://github.com/Circuitscape/Circuitscape.jl/tree/master/test/input).

## Features

### Solver Modes

Circuitscape supports two solver modes:

- **CG+AMG** (default) — an iterative solver using algebraic multigrid preconditioning. Scales well to large problems.
- **CHOLMOD** — a direct solver using [Cholesky decomposition](https://en.wikipedia.org/wiki/Cholesky_decomposition) via the [SuiteSparse](https://people.engr.tamu.edu/davis/suitesparse.html) library. Can be significantly faster for smaller problems, but memory use grows quickly with problem size due to [fill-in](https://algowiki-project.org/en/Cholesky_method#Reordering_to_reduce_the_number_of_fill-in_elements).

To use CHOLMOD, add to your INI file:
```
solver = cholmod
```

### Parallel Computing

Circuitscape supports parallel computation on Linux, macOS, and Windows using Julia's
built-in distributed computing.

### Single Precision

Circuitscape can run in single precision for reduced memory use at the cost of solution
accuracy. Add to your INI file:
```
precision = single
```

## Performance

Circuitscape.jl is up to **4x faster** than the legacy Python version (Circuitscape v4.0.5)
on 16 parallel processes, with the CHOLMOD solver providing the best performance on
suitable problem sizes.

<img src="docs/src/benchmark/benchmark.png" width=650 height=450>

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
