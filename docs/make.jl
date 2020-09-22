using Circuitscape, Documenter

makedocs(
    sitename = "Circuitscape.jl Documentation",
    pages = ["Home" => "index.md",
             "User Guide" => "usage.md",
             "Inputs, Outputs and Options" => "options.md", 
             "On Solvers and Computation Time" => "compute.md", 
             "Calling Circuitscape from other Programs" => "otherprog.md", 
             "Logging Options" => "logging.md"
    ]
)

deploydocs(
    repo = "github.com/Circuitscape/Circuitscape.jl.git",
    devbranch = "master",
    devurl = "latest"
)