# Calling Circuitscape from Other Programs

Circuitscape is a Julia package with one entry point, `Circuitscape.compute`,
which takes either the path of an INI file or a dictionary of INI keys and
string values. Nothing is exported, so every call is written fully qualified
(`Circuitscape.compute`, `Circuitscape.init_config`); `using Circuitscape:
compute` brings the bare name in if you prefer. Anything that can call Julia
can therefore run Circuitscape. Output files are
written next to `output_file`; `compute` also returns the result matrix
(pairwise resistances with the focal node IDs in the first row and column).

## From Julia

```julia
using Circuitscape
result = Circuitscape.compute("job.ini")
```

Or build the configuration in code, starting from the defaults:

```julia
using Circuitscape
cfg = Circuitscape.init_config()
cfg["habitat_file"] = "resistance_map.asc"
cfg["point_file"] = "focal_nodes.asc"
cfg["scenario"] = "pairwise"
cfg["output_file"] = "output/results.out"
cfg["write_cur_maps"] = "True"
result = Circuitscape.compute(cfg)
```

Values are the same strings you would write in an INI file (`"True"`, `"cholmod"`,
`"single"`, ...); unrecognised values are errors, see *Configuration Validation*
in [Inputs, Outputs and Options](options.md).

To run in parallel, start Julia with threads and turn parallelism on in the
configuration:

```bash
JULIA_NUM_THREADS=8 julia myscript.jl      # or: julia -t 8 myscript.jl
```

```julia
cfg["parallelize"] = "True"
```

GeoTIFF input or output (`write_as_tif = True`) needs the optional ArchGDAL
package, loaded before the run:

```julia
using Pkg; Pkg.add("ArchGDAL")   # once
using ArchGDAL, Circuitscape
Circuitscape.compute("job_with_tifs.ini")
```

ESRI ASCII grids (`.asc`) need no extra package.

## From R

Use the [JuliaCall](https://cran.r-project.org/package=JuliaCall) package. Set
`JULIA_NUM_THREADS` *before* `julia_setup()` if you want threads; the Julia
process is started once per R session.

```r
install.packages("JuliaCall")          # once
library(JuliaCall)

Sys.setenv(JULIA_NUM_THREADS = "8")    # optional, before julia_setup()
julia_setup()                          # julia_setup(installJulia = TRUE) installs Julia if needed

julia_install_package_if_needed("Circuitscape")   # once
julia_library("Circuitscape")

result <- julia_call("Circuitscape.compute", "job.ini")   # returns the resistance matrix
```

`julia_call` and `julia_eval` resolve names in Julia's `Main` module, and
Circuitscape exports nothing, so the function is always written
`Circuitscape.compute` (alternatively run `julia_command("using Circuitscape: compute")`
once and then call `julia_call("compute", ...)`).

A configuration can be built in R as a named list of strings and passed as a
Julia `Dict`:

```r
cfg <- julia_call("Circuitscape.init_config")
julia_assign("cfg", cfg)
julia_command('cfg["habitat_file"] = "resistance_map.asc"')
julia_command('cfg["point_file"] = "focal_nodes.asc"')
julia_command('cfg["scenario"] = "pairwise"')
julia_command('cfg["output_file"] = "output/results.out"')
julia_command('cfg["parallelize"] = "True"')
result <- julia_eval("Circuitscape.compute(cfg)")
```

For GeoTIFF rasters add `julia_install_package_if_needed("ArchGDAL")` and
`julia_library("ArchGDAL")` before `julia_library("Circuitscape")`.

## From Python

Use the [juliacall](https://pypi.org/project/juliacall/) package
(`pip install juliacall`). It installs Julia on first use if none is found.
Threads are set with the `PYTHON_JULIACALL_THREADS` environment variable
before `juliacall` is imported.

```python
import os
os.environ["PYTHON_JULIACALL_THREADS"] = "8"   # optional, before the import

from juliacall import Main as jl

jl.seval('import Pkg; Pkg.add("Circuitscape")')  # once
jl.seval("using Circuitscape")

result = jl.Circuitscape.compute("job.ini")      # Julia Matrix; np.asarray(result) converts
```

As in R, nothing is exported: functions are reached through the module object,
`jl.Circuitscape.compute`, `jl.Circuitscape.init_config`.

Building the configuration in Python:

```python
cfg = jl.Circuitscape.init_config()
cfg["habitat_file"] = "resistance_map.asc"
cfg["point_file"] = "focal_nodes.asc"
cfg["scenario"] = "pairwise"
cfg["output_file"] = "output/results.out"
cfg["parallelize"] = "True"
result = jl.Circuitscape.compute(cfg)
```

For GeoTIFF rasters run `jl.seval('import Pkg; Pkg.add("ArchGDAL")')` once and
`jl.seval("using ArchGDAL")` before `using Circuitscape`.

## Downstream packages

[Omniscape.jl](https://github.com/Circuitscape/Omniscape.jl) is the main
package built on Circuitscape. It runs omnidirectional connectivity analyses
by calling `Circuitscape.compute_omniscape_current` on in-memory arrays for
every moving window, and is the reference for embedding Circuitscape in
another Julia package. See the [API Reference](api.md).
