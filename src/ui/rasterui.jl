function raster_ui()

    # First Node
    rmap = vbox(Node(:div, "Select raster map from file: "), 
                 Node(:input, style = Dict(:style => "margin: 12px")))

    rmap
end
