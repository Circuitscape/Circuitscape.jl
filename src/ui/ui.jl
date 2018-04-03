using Blink
using WebIO
using Tachyons
using CSSUtil
using JSExpr

include("pairwise_ui.jl")
include("advanced_ui.jl")
include("output_ui.jl")

w = Window()

function showsome(uis, which)
    s = Scope()
    s["visible"] = which
    s.dom = Node(:div, id="cont", uis...)
    onjs(s["visible"], JSExpr.@js function (visbl)
                 @var cont = this.dom.querySelector("cont")
                 @var cs = cont.children
                 for i = 1:cs.length
                     if visbl.indexOf(i) >= 0
                         cs[i].style.display = "block"
                     else          
                         cs[i].style.display = "none"
                     end
                 end
             end)
    s
end

function generate_ui(w)

    heading = Node(:div, tachyons_css, "Circuitscape 5.0") |> 
                    class"f-subheadline lh-title tc red"

    section1 = Node(:div, tachyons_css, "Data Type and Modelling Mode") |> 
                    class"f4 lh-title"

    dt = get_data_type()
    # dt["value"][] = "Raster"

    mode = Observable{Any}("")
    points_input = Observable{Any}(Node(:div))

    # Next drop down
    mod_mode = map(dt["value"]) do v
        @show v
        points_input[] = pairwise_input_ui()
        if v == "Network"
            get_mod_mode_network()
        else
            get_mod_mode_raster()
        end
    end

    # Get the input raster/graph 
    input_section = Node(:div, tachyons_css, "Input Resistance Data") |> 
                    class"f4 lh-title"
    input1, input2 = input_ui()
    input = vbox(input1, 
                 input2)
    
    input_graph = Observable("")
    is_res = Observable(false)
    on(input1["filepath"]) do x
        input_graph[] = x
    end
    on(input2["check"]) do x
           is_res[] = x
    end
    
    # Focal points or advanced mode
    pair = pairwise_input_ui()
    adv = advanced_input_ui()
    
    on(mod_mode) do x
        v = x["value"]
        on(v) do s
            if s == "Advanced"
                points_input[] = adv
            else
                points_input[] = pair
            end
        end
    end
    dt["value"][] = "Raster"

    #=focal = Observable("")
    source = Observable("")
    ground = Observable("")
    on(points_input) do x
        on(x["focal"]) do y
            focal[] = y
            @show focal[]
        end
        on(x["source"]) do z
            source[] = z
            @show source[]
        end
        on(x["ground"]) do v
            ground[] = v
            @show ground[]
        end
    end=#
    
    # Output options
    write_cur_maps = Observable(false)
    write_volt_maps = Observable(false)
    out = Observable("")
    output = output_ui()

    on(output["cur"]) do x
        write_cur_maps[] = x
    end
    on(output["volt"]) do x
        write_volt_maps[] = x
    end
    on(output["out"]) do x
        out[] = x
    end

    page = vbox(heading, 
                hline(style = :solid, w=5px)(style = Dict(:margin => 20px)), 
                section1,
                dt, 
                mod_mode, 
                hline(style = :solid, w=3px)(style = Dict(:margin => 10px)),
                input_section,
                input, 
                hline(style = :solid, w=3px)(style = Dict(:margin => 10px)),
                points_input,
                hline(style = :solid, w=3px)(style = Dict(:margin => 10px)),
                output)|> class"pa3 system-sans-serif"

    body!(w, page)

end

function get_data_type()

    data_type_prompt = "Step 1: Choose your input data type: "

    # Select between raster and network
    data_type = vbox(Node(:div, data_type_prompt, 
                          attributes = Dict(:style => "margin-top: 12px")) |> class"b",
                 Node(:select, "Select Data Type", 
                 Node(:option, "Raster"), 
                 Node(:option, "Network"), id = "dt", 
                 attributes = Dict(:style => "margin-top: 12px; margin-bottom: 12px")))

    s = Scope()
    s.dom = data_type
    s["value"] = Observable("Raster")
    onimport(s, JSExpr.@js function ()
                 @var el = this.dom.querySelector("#dt")
                 el.onchange = (function ()
                        $(s["value"])[] = el.value
                    end)
             end)

    s
end

function get_mod_mode_network()
    
    mod_mode_prompt = "Step 2: Choose a Modelling Mode: "
    mod_mode = vbox(Node(:div, mod_mode_prompt, 
                          attributes = Dict(:style => "margin-top: 12px")) |> class"b",
                 Node(:select, "Select Modelling Mode", 
                 Node(:option, "Pairwise"), 
                 Node(:option, "Advanced"), id = "modelling", 
                 attributes = Dict(:style => "margin-top: 12px; margin-bottom: 12px")))

    s = Scope()
    s.dom = mod_mode
    s["value"] = Observable("Pairwise")
    onimport(s, JSExpr.@js function ()
                 @var el = this.dom.querySelector("#modelling")
                 el.onchange = (function ()
                        $(s["value"])[] = el.value
                    end)
             end)

    s
end

function get_mod_mode_raster()
    
    mod_mode_prompt = "Step 2: Choose a Modelling Mode: "
    mod_mode = vbox(Node(:div, mod_mode_prompt, 
                          attributes = Dict(:style => "margin-top: 12px")) |> class"b",
                 Node(:select, "Select Modelling Mode", 
                 Node(:option, "Pairwise"), 
                 Node(:option, "Advanced"), 
                 Node(:option, "One To All "), 
                 Node(:option, "All To One"), id = "modelling", 
                 attributes = Dict(:style => "margin-top: 12px; margin-bottom: 12px")))

    s = Scope()
    s.dom = mod_mode
    s["value"] = Observable("Pairwise")
    onimport(s, JSExpr.@js function ()
                 @var el = this.dom.querySelector("#modelling")
                 el.onchange = (function ()
                        $(s["value"])[] = el.value
                    end)
             end)
    s
end
generate_ui(w)
