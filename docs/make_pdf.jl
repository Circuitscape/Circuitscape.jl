using Circuitscape, Documenter, tectonic_jll

include("pages.jl")

makedocs(
    sitename = "Circuitscape.jl",
    format = Documenter.LaTeX(platform = "tectonic", tectonic = tectonic_jll.tectonic()),
    warnonly = [:cross_references],
    doctest = false,   # make.jl already runs the doctests
    pages = PAGES
)
