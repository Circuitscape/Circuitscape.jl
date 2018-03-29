function advanced_input_ui()

    # Title of section
    title = Node(:div, tachyons_css, "Advanced Mode Options") |> class"f4 lh-title"

    # Sources
    sources = vbox(Node(:div, "Current source file: ", 
                      attributes = Dict(:style => "margin-top: 12px")) |> class"b",
                 Node(:input, attributes = Dict(:type => :file, 
                                                :style => "margin-top: 12px")))

    # Grounds
    grounds = vbox(Node(:div, "Ground point file: ", 
                      attributes = Dict(:style => "margin-top: 12px")) |> class"b",
                 Node(:input, attributes = Dict(:type => :file, 
                                                :style => "margin-top: 12px")))
    vbox(title,
         sources,
         grounds)
end
