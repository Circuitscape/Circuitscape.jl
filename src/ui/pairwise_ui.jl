
function input_ui()
    

    graph = vbox(Node(:div, "Raster resistance map or network/graph: ", 
                      attributes = Dict(:style => "margin-top: 12px")) |> class"b",
                 Node(:input, attributes = Dict(:type => :file, 
                                                :style => "margin-top: 12px")))

    graph_is_res = hbox(Node(:input, "graph_is_res",
                            attributes = Dict(:type => :checkbox, 
                           :style => "margin-top: 12px; margin-right: 5px")), 
                        Node(:div, "Data represents resistances instead of conductances", 
                             attributes = Dict(:style => "margin-top: 12px")))

    vbox(graph,
         graph_is_res)
end

function pairwise_input_ui()
    additional_input = Node(:div, tachyons_css, "Pairwise Mode Options") |> class"f4 lh-title"
    focal = vbox(Node(:div, "Select focal node locations from file: ", 
                      attributes = Dict(:style => "margin-top: 12px")),
                 Node(:input,  attributes = Dict(:type => :file, 
                                                 :style => "margin-top: 12px")))
    vbox(additional_input, 
         focal)
end

function output_nodes()
    output = vbox(Node(:div, "Enter base name for output file: ",
                       attributes = Dict(:style => "margin-top: 12px")),
                  Node(:input, attributes = Dict(:type => :text, 
                                                 :style => "margin-top: 12px")))

    el = vbox(graph,
              graph_is_res, 
              focal, 
              output)

end
