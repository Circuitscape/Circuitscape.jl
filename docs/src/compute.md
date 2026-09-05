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

## Convergence Checks

Every linear solve, whichever solver produced it, is accepted only if its true
relative residual `‖Gv − b‖ / ‖b‖` is below `residual_tolerance` (by default
`1e-4` in double precision and `1e-3` in single). For `cg+amg` there is a subtlety: Krylov's stopping test is on the
*preconditioned* residual, and when the multigrid preconditioner is a poor fit
for a particular component the two norms can disagree — CG reports success
while the true residual is still too large. Circuitscape then tightens the
Krylov tolerance and retries (up to three attempts) before failing, so a single
awkward component does not abort a long run. The retry is logged as a warning;
if you see many of them, the grid probably has extreme resistance contrasts
worth examining.

Single precision (`precision = single`) halves memory at a cost in accuracy
that depends on the problem: on the largest raster in the test suite,
effective resistances differ from the double precision run by about 1%, and
the gap will be larger on bigger grids or with stronger resistance contrasts.
The difference does not change when `residual_tolerance` is tightened, so it
is a property of the Float32 arithmetic rather than of solver convergence;
verify against a double precision run before relying on single precision
results.

## What Makes a Run Fast

The facts below are measured; each cites the pull request that carries the
full table. Timings were taken on the development machines named in those
PRs and are for relative comparison only.

### Circuitscape 5.17 vs 6.0

Three end-to-end runs of v5.17.1 (from the General registry) against 6.0
(master at c867d5a), on an Apple M2 Max with Julia 1.12.7, 8 threads
(`JULIA_NUM_THREADS=8`, `parallelize = True`), solver `cg+amg`, double
precision, eight-neighbour connectivity, random integer resistances 1–10,
16 single-cell focal points, ESRI ASCII inputs. Each cell is a fresh Julia
process: one warm-up `compute` call, then one timed call (a single run, not a
best-of-N). "Allocated" is Julia's total allocation for that call
(`@allocated`); "peak RSS" is `Sys.maxrss()` at process end.

| Case | Metric | v5.17.1 | 6.0 |
|------|--------|---------|-----|
| pairwise, default settings (no maps), 2000×2000 (4M cells) | time | 72.3 s | 14.9 s |
| | allocated | 27.1 GiB | 22.1 GiB |
| | peak RSS | 9.7 GiB | 9.9 GiB |
| pairwise + current maps, 1000×1000, all 119 pairs (include file defeats the shortcut) | time | 61.5 s | 28.7 s |
| | allocated | 128.4 GiB | 78.1 GiB |
| | peak RSS | 5.7 GiB | 4.1 GiB |
| one-to-all, 1000×1000 | time | 6.1 s | 5.7 s |
| | allocated | 45.5 GiB | 42.4 GiB |
| | peak RSS | 11.7 GiB | 11.0 GiB |

The default path is 4.9× faster because the resistance-shortcut solves now
run in parallel (they were serial in 5.x) and the graph-construction and
postprocessing overheads were removed. With maps on, the gain (2.1×, 39% less
allocation) comes from computing node currents without sparse temporaries.
One-to-all is essentially unchanged because its cost is the per-focal-point
AMG solve, which 6.0 did not change. Peak RSS on the default path is flat
because it is dominated by the AMG hierarchy for a 4M-cell Laplacian, which
is identical in both versions. These are single runs on one machine and will
vary; the per-PR tables in the sections below are the component
measurements. Load time also dropped: `using Circuitscape` went from 2.71 s
to about 1.2 s, and the resolved dependency count from 195 to 104 packages
(ArchGDAL optional, Graphs removed).

### Ask only for what you need

- **The resistance shortcut.** In raster pairwise mode, when no current,
  voltage, cumulative or maximum map is requested and no include/exclude
  file is used, Circuitscape solves one linear system per focal node of a
  component instead of one per pair, and derives every pairwise resistance
  from those voltages: `n − 1` solves rather than `n(n − 1)/2`. The log says
  `Triggering resistance calculation shortcut` when it is active. Requesting
  any map, or an include/exclude file, switches to one solve per pair.
- **No maps, no postprocessing.** Since 6.0 ([#500](https://github.com/Circuitscape/Circuitscape.jl/pull/500)),
  when none of the four map options is set the per-pair postprocessing (node
  currents, a full-grid current map, the locked accumulation into a
  cumulative map) is skipped entirely. On a 400×400 grid, 16 focal points,
  119 pairs with the shortcut defeated, 8 threads: 4.96 s → 3.31 s and
  20.28 GiB → 3.49 GiB allocated per run. On a 1M-cell, 629-pair run this
  path had been 19% of profile samples and about 1 GB of allocation per pair.

### Threads

Set `parallelize = True` and start Julia with threads (`julia -t N`). Since
[#491](https://github.com/Circuitscape/Circuitscape.jl/pull/491) the `cg+amg`
pair solves are scheduled in equal-sized chunks, about four per thread,
rather than one task per source node. That matters most for the default
no-maps path: the shortcut anchors every solve at one source, so before 6.0
all of them ran on a single task whatever `-t` was. Measured on random
resistance grids, 8 threads, best of 3:

| Grid | Case | Before | After |
|------|------|--------|-------|
| 400×400 | shortcut (default), 16 focal points | 2.37 s | 0.64 s |
| 400×400 | all pairs, 16 focal points (119 pairs) | 5.20 s | 4.66 s |
| 2000×2000 (4M cells) | shortcut (default), 16 focal points (15 solves) | 73.9 s | 20.4 s |
| 2000×2000 | all pairs, 9 focal points (35 solves) | 59.6 s | 48.2 s |

Large all-pairs runs are memory-bandwidth bound, so extra threads help them
less. In one-to-all and all-to-one each focal node is an independent task
([#505](https://github.com/Circuitscape/Circuitscape.jl/pull/505) also stopped
rebuilding the graph per focal node when an include file is used with
single-cell focal points: 500×500 grid, 30 points, 5.33 GiB → 4.35 GiB
allocated).

### Memory

Load-time memory decides the largest grid that fits, and 6.0 removed several
copies of the graph from that path:

| Change | PR | Measured |
|--------|----|----------|
| ESRI ASCII grids read and written natively, no GDAL and no lock | [#492](https://github.com/Circuitscape/Circuitscape.jl/pull/492) | 100 MB `.asc`: write 2.20 s → 0.58 s, read at parity with GDAL; files ~17% smaller |
| GDAL optional (ArchGDAL is a package extension) | [#493](https://github.com/Circuitscape/Circuitscape.jl/pull/493) | default install resolves 110 packages instead of 195 (104 once Graphs was removed in #495) |
| Connected components computed on the CSC structure; no Graphs adjacency copy | [#495](https://github.com/Circuitscape/Circuitscape.jl/pull/495) | 2000×2000 grid: 1,574 MB → 140 MB allocated, 458 MB → 0 retained, 0.73 s → 0.32 s |
| `construct_graph` sizes its buffers exactly and emits both edge orientations; no `a + a'` | [#504](https://github.com/Circuitscape/Circuitscape.jl/pull/504) | 2000×2000 grid: 2,536 MiB → 1,311 MiB allocated (−48%) for the same 0.5 GB matrix |
| Node currents in one sweep over the matrix, no sparse temporaries | [#501](https://github.com/Circuitscape/Circuitscape.jl/pull/501) | 1000×1000 grid: 967 MB → 8 MB allocated per pair, about 10× faster |
| Network cumulative branch currents indexed once | [#497](https://github.com/Circuitscape/Circuitscape.jl/pull/497) | accumulation per pair is O(E) instead of O(E²) |

### Direct or iterative

- `cg+amg` (default) keeps memory proportional to the number of nonzeros in
  the Laplacian and parallelizes the pair solves. It is the only choice for
  the largest grids; a run with more than 5 million cells and
  `solver = cholmod` logs a warning to that effect.
- The direct solvers (`cholmod`, `accelerate`, `pardiso`) factorize each
  connected component once and then solve `cholmod_batch_size` right-hand
  sides at a time, so with many focal pairs on a small or medium grid they are
  usually faster. Their memory is dominated by fill-in in the factor, which
  grows faster than the grid. The batch size trades memory for fewer solver
  calls: the right-hand-side and solution arrays are `n × batch` dense
  matrices, so on a large grid a smaller batch may be necessary.

### Single precision

`precision = single` halves the memory of every vector and matrix. The cost
is accuracy: see the note under *Convergence Checks* above for what to
expect and how to check it.

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

## GeoTIFF Support

Circuitscape reads and writes ESRI ASCII grids (`.asc`, optionally gzipped)
natively and needs no extra packages for them. GeoTIFF input, and GeoTIFF
output via `write_as_tif = True`, use GDAL through the optional
[ArchGDAL.jl](https://github.com/yeesian/ArchGDAL.jl) package extension.
Install it once and load it before running:

```julia
using Pkg
Pkg.add("ArchGDAL")
```

```julia
using ArchGDAL
using Circuitscape
Circuitscape.compute("config.ini")   # may now name .tif rasters, or set write_as_tif = True
```

Keeping GDAL optional cuts the default install by a large set of binary
dependencies (GDAL, PROJ, GEOS, HDF5, NetCDF, …) and their precompile time. A
configuration that names a GeoTIFF without ArchGDAL loaded is refused before
any data is read, with a message saying exactly this.

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
Circuitscape.compute("config.ini")  # with solver = pardiso in the INI file
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
Circuitscape.compute("config.ini")  # with solver = accelerate in the INI file
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
