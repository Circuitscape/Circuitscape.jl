using REPL
using REPL.TerminalMenus

mutable struct History
    path::String
end

h = History(homedir())
update!(h, path) = (h.path = path)

function filepicker()
    println()
    r = readdir()
    # length(r) == 1 && return readdir()[1]
    path = h.path
    printstyled("""
          Select a path:
          """, bold=true)
    pick(path)
end

function pick(path) 
    isfile(path) && return normpath(path)
    r = alldir(path)
    l = request(RadioMenu(r, pagesize=10)) 
    str = r[l]
    l == 1 && return manualfilepicker()
    l == 2 && (str = "..")
    println("---------------------------------")
    printstyled("Selected path: ", bold=true) 
    print("$(normpath(joinpath(path,str)))\n")
    println("---------------------------------")
    update!(h, path)
    pick(joinpath(path, str))
end

alldir(path) = vcat("↩", readdir(path))

function manualfilepicker() 
    println()
    println("Enter path: ")
    path = readline(stdin)
    ispath(path) && return path
    println("Invalid path, please reenter!")
    manualfilepicker()
end
    
function folderpicker()
    path = h.path
    printstyled("""
            Select a folder: 
            """, bold=true)
    pickfolder(path)
end

function pickfolder(path)
    r = vcat("↩", "DONE", readfolders(path))
    l = request(RadioMenu(r, pagesize=10))
    str = r[l]
    l == 1 && (str == "..")
    str == "DONE" && return path
    println("---------------------------------")
    printstyled("Selected path: ", bold=true) 
    print("$(normpath(joinpath(path,str)))\n")
    println("---------------------------------")
    update!(h, path)
    pickfolder(joinpath(path, str))
end
readfolders(path) = readdir(path) |> x -> filter(!isfile, x)




