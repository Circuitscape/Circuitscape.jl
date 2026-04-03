# Calling Circuitscape from Other Programs

## From Julia

Circuitscape can be called directly from Julia:

```julia
using Circuitscape
result = compute("myconfig.ini")
```

You can also pass a configuration dictionary directly:

```julia
using Circuitscape
cfg = Circuitscape.init_config()
cfg["habitat_file"] = "resistance_map.asc"
cfg["point_file"] = "focal_nodes.asc"
cfg["scenario"] = "pairwise"
cfg["output_file"] = "output/results.out"
result = compute(cfg)
```

## From R

You can call Circuitscape from R using the [JuliaCall](https://cran.r-project.org/package=JuliaCall) R package.

## From Python

You can call Circuitscape from Python using the [juliacall](https://pypi.org/project/juliacall/) Python package.
