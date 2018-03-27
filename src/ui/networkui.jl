
function network_ui()
    
    # First Node
    graph = vbox(Node(:div, "Select pair list from file: "),
                 Node(:input, attributes = Dict(:type => :file, 
                                                :style => "margin-top: 12px")))

    graph_is_res = hbox(Node(:input, "graph_is_res",
                            attributes = Dict(:type => :checkbox, 
                           :style => "margin-top: 12px; margin-right: 5px")), 
                        Node(:div, "Data represents resistances instead of conductances", 
                             attributes = Dict(:style => "margin-top: 12px")))

    focal = vbox(Node(:div, "Select focal node locations from file: ", 
                      attributes = Dict(:style => "margin-top: 12px")),
                 Node(:input,  attributes = Dict(:type => :file, 
                                                 :style => "margin-top: 12px")))

    output = vbox(Node(:div, "Enter base name for output file: ",
                       attributes = Dict(:style => "margin-top: 12px")),
                  Node(:input, attributes = Dict(:type => :text, 
                                                 :style => "margin-top: 12px")))

    el = vbox(graph, 
              graph_is_res, 
              focal, 
              output)

end
