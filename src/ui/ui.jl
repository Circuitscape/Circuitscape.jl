using Blink
using WebIO
using Tachyons
using CSSUtil
using JSExpr

w = Window()

function generate_ui(w)

    heading = Node(:div, tachyons_css, "Circuitscape 5.0") |> 
                    class"f-subheadline tc red"

    line = Node(:hr)

    # Select between raster and network
    data_type = vbox(Node(:div, "Select Data Type: "),
                 Node(:select, "Select Data Type", 
                 Node(:option, "Raster"), 
                 Node(:option, "Network"), id = "dt", attributes = Dict(:style => "margin: 12px")))

    s = get_data_type!(data_type)
    ui1 = network_ui()
    ui2 = raster_ui()
    nextui = map(s["value"]) do v
        if v == "Network"
            ui1
        else
            ui2
        end
    end

    # boxes = [Node(:input, "Enter file path", attributes = Dict(:style => "margin: 12px")) for i = 1:10] 

    page = vbox(heading, 
                hline(style = :solid)(style = Dict(:margin => 20px)), 
                s, nextui) |> class"pa3 system-sans-serif"
    

    body!(w, page)

    dt
end

function get_data_type!(data_type)

    s = Scope()
    s.dom = data_type
    onimport(s, JSExpr.@js function ()
                 @var el = this.dom.querySelector("#dt")
                 el.onchange = (function ()
                        $(s["value"])[] = el.value
                    end)
             end)

    on(s["value"]) do x
        return s, x
    end
    s
end

generate_ui(w)
