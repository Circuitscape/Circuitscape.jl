function output_ui()
    
    # Title of section
    title = Node(:div, tachyons_css, "Output Options") |> class"f4 lh-title"

    base_name = vbox(Node(:div, "Base output file name: ", 
                      attributes = Dict(:style => "margin-top: 12px")) |> class"b",
                 Node(:input, attributes = Dict(:type => :text, 
                                                :style => "margin-top: 12px")))

    options = vbox(hbox(Node(:input, attributes = Dict(:type => :checkbox, 
                           :style => "margin-top: 12px; margin-right: 5px")), 
                        Node(:div, "Current maps", 
                             attributes = Dict(:style => "margin-top: 12px"))),
                  hbox(Node(:input, attributes = Dict(:type => :checkbox, 
                           :style => "margin-top: 12px; margin-right: 5px")), 
                        Node(:div, "Voltage maps", 
                             attributes = Dict(:style => "margin-top: 12px"))))

    vbox(title,
         base_name, 
         options)
end
