using Circuitscape, Documenter

include("pages.jl")

# Doctests run here, once; the PDF build skips them.
makedocs(
    sitename = "Circuitscape.jl Documentation",
    pages = PAGES
)

deploydocs(
    repo = "github.com/Circuitscape/Circuitscape.jl.git",
    devbranch = "master",
    devurl = "latest"
)
