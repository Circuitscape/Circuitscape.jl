using Blink
using WebIO
using Tachyons
using CSSUtil
using JSExpr

include("pairwise_ui.jl")
include("advanced_ui.jl")

w = Window()

function generate_ui(w)

    heading = Node(:div, tachyons_css, "Circuitscape 5.0") |> 
                    class"f-subheadline lh-title tc red"

    section1 = Node(:div, tachyons_css, "Data Type and Modelling Mode") |> 
                    class"f3 lh-title"

    dt = get_data_type()

    mod_mode_network = get_mod_mode_network()
    mod_mode_raster = get_mod_mode_raster()

    # Next drop down
    mod_mode = map(dt["value"]) do v
        if v == "Network"
            mod_mode_network
        else
            mod_mode_raster
        end
    end

    input_section = Node(:div, tachyons_css, "Input Resistance Data") |> 
                    class"f3 lh-title"
    
    #=input_data = map(mod_mVode["value"]) do v
        if v == "Pairwise"
            pairwise_input_ui()
        elseif v == "Advanced"
            advanced_input_ui()
        elseif v == "One To All"
            onetoall_input_ui()
        elseif v == "All To One"
            alltoone_inputui()
        end
    end=#

    input = input_ui()

    page = vbox(heading, 
                hline(style = :solid, w=5px)(style = Dict(:margin => 20px)), 
                section1,
                dt, 
                mod_mode, 
                hline(style = :solid, w=3px)(style = Dict(:margin => 10px)),
                input_section,
                input)|> class"pa3 system-sans-serif"

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
                 Node(:option, "Advanced"), id = "mod", 
                 attributes = Dict(:style => "margin-top: 12px; margin-bottom: 12px")))

    s = Scope()
    s.dom = mod_mode
    onimport(s, JSExpr.@js function ()
                 @var el = this.dom.querySelector("#mod")
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
                 Node(:option, "All To One"), id = "mod", 
                 attributes = Dict(:style => "margin-top: 12px; margin-bottom: 12px")))

    s = Scope()
    s.dom = mod_mode
    onimport(s, JSExpr.@js function ()
                 @var el = this.dom.querySelector("#mod")
                 el.onchange = (function ()
                        $(s["value"])[] = el.value
                    end)
             end)
    s
end
generate_ui(w)
