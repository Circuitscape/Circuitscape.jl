using Circuitscape, Documenter, tectonic_jll

makedocs(
    sitename = "Circuitscape.jl",
    format = Documenter.LaTeX(platform = "tectonic", tectonic = tectonic_jll.tectonic()),
    warnonly = [:cross_references],
    pages = ["Home" => "index.md",
             "User Guide" => "usage.md",
             "Inputs, Outputs and Options" => "options.md",
             "On Solvers and Computation Time" => "compute.md",
             "Calling Circuitscape from other Programs" => "otherprog.md",
             "Logging Options" => "logging.md"
    ]
)
