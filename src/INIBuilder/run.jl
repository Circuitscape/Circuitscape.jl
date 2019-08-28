cfg = Circuitscape.init_config()

state = Dict{Int, Symbol}()
for i = 1:10
    state[i] = Symbol("state" * string(i))
end

function step1()
    println()
    println("Step 1: Choose Data Type")
    dt = ["raster", "network"]
    n = request(RadioMenu(vcat(dt)))
    cfg["data_type"] = dt[n]
    step2()
end

function step2()
    println()
    println("Step 2: Choose Modelling Mode")
    is_raster = cfg["data_type"] == "raster"
    if is_raster 
        mm = ["PREVIOUS STEP", "pairwise", "advanced", "one-to-all", 
                "all-to-one"]
    else
        mm = ["PREVIOUS STEP", "pairwise", "advanced"]
    end
    n = request(RadioMenu(vcat(mm)))
    n == 1 && step1()
    cfg["scenario"] = mm[n]
    step3()
end

function step3()
    println()
    println("Step 3a: Enter path to habitat file")
    opt = ["PREVIOUS STEP", "Enter path manually", "Use filepicker"]
    n = request(RadioMenu(opt))
    n == 1 && step2()
    path = n == 2 ? manualfilepicker() : filepicker()
    cfg["habitat_file"] = path
    println()
    println("Is this a resistance or conductance file?")
    is_res = ["PREVIOUS STEP", "resistances", "conductance"]
    n = request(RadioMenu(is_res))
    n == 1 && step2()
    cfg["habitat_map_is_resistances"] = is_res[n]
    step4() 
end

function step4()
    println()
    is_advanced = cfg["scenario"] == "advanced"
    if !is_advanced
        println()
        println("Step 4: Enter path to focal nodes:")
        opt = ["PREVIOUS STEP", "Enter path manually", "Use filepicker"]
        n = request(RadioMenu(opt))
        n == 1 && step3()
        path = n == 2 ? manualfilepicker() : filepicker()
        cfg["point_file"] = path
    else
        println()
        println("Step 4a: Enter path to source file") 
        opt = ["PREVIOUS STEP", "Enter path manually", "Use filepicker"]
        n = request(RadioMenu(opt))
        n == 1 && step3()
        path = n == 2 ? manualfilepicker() : filepicker()
        cfg["source_file"] = path
        println()
        println("Step 4b: Enter path to ground file") 
        n = ["PREVIOUS STEP", "Enter path manually", "Use filepicker"]
        n == 1 && step3()
        path = n == 2 ? manualfilepicker() : filepicker()
        cfg["ground_file"] = path
    end
    step5()
end

function step5()
    println()
    println("Step 5: Choose solver")
    opt = ["PREVIOUS STEP", "cg+amg", "cholmod"]
    n = request(RadioMenu(opt))
    n == 1 && step4()
    cfg["solver"] = opt[n]
    println()
    step6()
end

function step6()
    println()
    println("Step 6: Choose number of parallel processes")
    opt = collect(1:Sys.CPU_THREADS) |> x -> string.(x)
    n = request(RadioMenu(opt))
    n > 1 && (cfg["parallelize"] = "true")
    cfg["max_parallel"] = string(n)
    step7()
end

function step7()
    println()
    println("Step 7: Choose outputs")
    opt = ["PREVIOUS STEP", "Pick outputs"]
    n = request(RadioMenu(opt))
    n == 1 && step6()
    opt = ["Current maps", "Voltage maps"]
    n = request(MultiSelectMenu(opt))
    1 in n && (cfg["write_cur_maps"] = "true")
    2 in n && (cfg["write_volt_maps"] = "true")
    step8()
end

function step8()
    println()
    println("Step 8: Choose output file name")
    opt = ["PREVIOUS STEP", "Enter output file name"]
    n = request(RadioMenu(opt))
    n == 1 && step7() 
    name = readline(stdin) 
    step9(name)
end

function step9(name)
    println()
    println("Step 9: Choose output folder")
    opt = ["PREVIOUS STEP", "Enter path manually", "Use folderpicker"]
    n = request(RadioMenu(opt))
    n == 1 && step3()
    path = n == 2 ? manualfolderpicker() : folderpicker()
    cfg["output_file"] = normpath(joinpath(path, name))
    step10(name, path)
end

function step10(name, path)
    println()
    println("Step 10: Would you like to run Circuitscape?")
    opt = ["Yes", "Later"]
    n = request(RadioMenu(opt))
    n == 1 && Circuitscape.compute(cfg)
    if n == 2 
        println()
        Circuitscape.write_config(cfg)
        println("$name.ini written to $(normpath(path))")
        println()
        println("Would you to build another problem?")
        l = RadioMenu(["Yes", "No"])
        n = request(l)
        n == 1 && step1()
        n == 2 && return
    end
end
    
function start()
    logo = raw"""

    `-````````````````````````````````````````````````
    -dd```````````````````````````````````````````````
    `.dh:``````````````````````````````````````.``````
    ```:ddddddddddhhdd+:-..``.-:ddddddhhhhddddhdd.````
    ``````-dddddddd.dddddddddhhhhhddddddddddddddddd```
    ````````````````ddddddddddddddddddddddddddddhd.```
    ```````````.ddddddddddddddddddddddddddddhhhd.`````
    ````````.ddddddddddddd```````````/dddddd/dddd:````
    ````.ddddddddd-----.``````````````````````````````
    ````.dddddd.````.ooo```````/o/```````:oo-`````````
    ``.````````````/oooo/o```oo/oo/o```oooo.:/o```````
    `ooo///o.````:o/`````oooo/`````-oooo.`````ooooooo-
    ``.`````/oo/o/```````oo/````````oo.````````````..`
    `````````oo/````````oo````````.o.`````````````````


    """
    blue    = "\033[34m"
    red     = "\033[31m"
    green   = "\033[32m"
    magenta = "\033[35m"
    normal  = "\033[0m\033[0m"
    
    welcome_message = raw"""
    Welcome to Circuitscape. 
    """
    logo = replace(logo, "d" => "$(red)d$(normal)")
    logo = replace(logo, "h" => "$(red)h$(normal)")
    logo = replace(logo, "o" => "$(magenta)o$(normal)")
    print(logo)
    print(welcome_message)

    step1()
end
