function advanced_input_ui()

    # Title of section
    title = Node(:div, tachyons_css, "Advanced Mode Options") |> class"f4 lh-title"

    # Sources
    sources = vbox(Node(:div, "Current source file: ", 
                      attributes = Dict(:style => "margin-top: 12px")) |> class"b",
                 Node(:input, id = "source", attributes = Dict(:type => :file, 
                                                :style => "margin-top: 12px")))

    # Grounds
    grounds = vbox(Node(:div, "Ground point file: ", 
                      attributes = Dict(:style => "margin-top: 12px")) |> class"b",
                 Node(:input, id = "ground", attributes = Dict(:type => :file, 
                                                :style => "margin-top: 12px")))
    adv = vbox(title,
             sources,
             grounds)

    # Define scope
    #=s = Scope()
    s.dom = adv

    onimport(s, JSExpr.@js function ()
                 @var el1 = this.dom.querySelector("#source")
                 @var el2 = this.dom.querySelector("#ground")
                 el1.onchange = (function ()
                         $(s["source"])[] = el1.files[0].path
                         $(s["ground"])[] = el2.files[0].path
                 end)
                 el2.onchange = (function ()
                         $(s["source"])[] = el1.files[0].path
                         $(s["ground"])[] = el2.files[0].path
                 end)
                end)=#
    adv
end
