module CircuitscapeREPLExt

# The interactive INI builder. Lives in an extension on the REPL stdlib so
# that scripts, Omniscape and R callers never load REPL; an interactive Julia
# session has it loaded already, so Circuitscape.start() works there with no
# extra step.

using REPL
using REPL.TerminalMenus
using Circuitscape

include("filepicker.jl")
include("run.jl")

Circuitscape.start() = start()

end # module CircuitscapeREPLExt
