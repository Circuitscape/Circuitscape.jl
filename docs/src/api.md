# API Reference

```@meta
CurrentModule = Circuitscape
```

Circuitscape exports nothing: every function, `compute` included, is reached
as `Circuitscape.name` (or brought in explicitly with
`using Circuitscape: compute`). The names in the **Public API** section are covered by
semantic versioning: they keep working, with the documented meaning, within a
major version. Anything not listed there is internal and may change in any
release. [Omniscape.jl](https://github.com/Circuitscape/Omniscape.jl) depends
on exactly three of these, `compute_omniscape_current`, `init_config` and
`update!`, and they are stable.

## Public API

```@docs
Circuitscape.compute
Circuitscape.init_config
Circuitscape.update!
Circuitscape.compute_omniscape_current
Circuitscape.start
```

### Configuration

A configuration is a `Dict{String,String}` of INI keys and values while it is
being built (`init_config`, `update!`), and a `CSConfig` once parsed.
`compute` accepts the dictionary form and does the parsing and validation
itself; these functions expose the same steps for scripts that want to read,
check or write INI files without running them.

```@docs
Circuitscape.parse_config
Circuitscape.CSConfig
Circuitscape.validate
Circuitscape.write_config
```

## Internals

!!! warning "Internal"
    Nothing below is part of the public API. These functions and types carry
    no semantic-versioning guarantee and may be renamed, changed or removed
    in any release. They are documented for contributors and for readers who
    want to follow a run through the code.

### Problem geometry

After loading, a raster landscape and a network are the same problem: a
symmetric Laplacian over graph nodes, its connected components, and focal
nodes. The one thing that differs downstream is where a node's voltage or
current goes in the output, and that is what a `Geometry` knows.

```@docs
Circuitscape.Geometry
Circuitscape.RasterGeometry
Circuitscape.NetworkGeometry
Circuitscape.restrict
```

### Problem construction and mode drivers

`load_data` reads a `RasterData` or `NetworkData`; `build_problem` turns it
into a `GraphProblem` (pairwise) or an `AdvancedProblem` (advanced) with the
matching geometry; `run_pairwise` and `run_advanced` are the mode drivers, and
one-to-all/all-to-one build an `AdvancedProblem` per focal node and call
`advanced_kernel` on it.

```@docs
Circuitscape.build_problem
Circuitscape.run_pairwise
Circuitscape.run_advanced
Circuitscape.AdvancedProblem
Circuitscape.advanced_kernel
```

### Pairwise driver

The pairwise driver is shared by every solver. `solve` iterates the connected
components of a `GraphProblem`, calling `prepare!` once per component (an AMG
preconditioner or a factorization), `pair_jobs` to enumerate the `PairJob`s
of that component, and `solve_pairs!` to run them, which hands every solved
pair to `handle_pair` to store the resistance and to `postprocess` to
accumulate and write maps through the `Cumulative` accumulators. Only
`prepare!` and `solve_pairs!` differ between the iterative and direct paths.

```@docs
Circuitscape.GraphProblem
Circuitscape.Cumulative
Circuitscape.PairJob
Circuitscape.pair_jobs
Circuitscape.get_num_pairs
Circuitscape.solve
Circuitscape.prepare!
Circuitscape.solve_pairs!
```

### Solvers and convergence

```@docs
Circuitscape.Solver
Circuitscape.AMGSolver
Circuitscape.CholmodSolver
Circuitscape.PardisoSolver
Circuitscape.AccelerateSolver
Circuitscape.get_solver
Circuitscape.residual_tolerance
Circuitscape.solve_linear_system
Circuitscape.refine_columns!
```

### Input and output

```@docs
Circuitscape.read_raster
Circuitscape.write_raster
Circuitscape.connected_components
```
