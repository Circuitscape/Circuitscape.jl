# Upgrading to 6.0

Circuitscape.jl 6.0 is the first release since 5.17.1 and collects everything
merged since then. Most INI files run unchanged. This page lists the changes
that can require action when upgrading.

## Breaking changes

- **`compute` is no longer exported.** Nothing is: call
  `Circuitscape.compute(...)`, or add `using Circuitscape: compute` to keep the
  bare name.
- **`start` is no longer exported.** Call the interactive INI builder as
  `Circuitscape.start()`; a bare `start()` after `using Circuitscape` is an
  `UndefVarError`.
- **Configuration values are validated.** Unrecognised strings for `solver`,
  `scenario`, `data_type`, `precision`, `log_level` and `remove_src_or_gnd`,
  non-boolean values for boolean options, a missing `scenario`, input files
  that do not exist and a missing output directory are now errors, reported
  together before any data is read, instead of silently falling back to a
  default. Run each INI file once and fix what the message names; see
  *Configuration Validation* in [Inputs, Outputs and Options](options.md).
- **GeoTIFF read and write require the optional ArchGDAL package.** Run
  `Pkg.add("ArchGDAL")` once, then `using ArchGDAL` before `using Circuitscape`
  whenever a job names a `.tif` or sets `write_as_tif = True`. ESRI ASCII grids
  need nothing extra.
- **`low_memory_mode` and `preemptive_memory_release` set to `True` are
  errors.** These Circuitscape 4 options were never implemented here; remove
  them or set them to `False`.

## Changed results

- **One-to-all and all-to-one variable source strengths are matched by focal
  node ID.** Previously the file's rows were applied in order, so a file not
  sorted by ID, or one that skipped an ID or listed an ID absent from the
  point file, gave strengths to the wrong nodes or failed. Sorted files that
  list exactly the focal nodes give the same results as before.
- **All-to-one with an include/exclude file:** a focal node excluded from a
  pair with the current ground no longer injects current into that solve.
  Previously it injected 1 A (and Circuitscape 4.0.5 does the same when
  variable source strengths are used), so current maps change wherever the
  file excludes pairs.

Everything else in 6.0 is additive or internal; see the
[release notes](https://github.com/Circuitscape/Circuitscape.jl/releases) for
the full list, and *What Makes a Run Fast* in
[On Solvers and Computation Time](compute.md) for the performance changes.
