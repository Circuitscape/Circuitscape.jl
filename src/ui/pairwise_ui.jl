
function input_ui()
    

    graph = vbox(Node(:div, "Raster resistance map or network/graph: ", 
                      attributes = Dict(:style => "margin-top: 12px")) |> class"b",
                 Node(:input, id = "input", attributes = Dict(:type => :file, 
                                                :style => "margin-top: 12px")))

    graph_is_res = hbox(Node(:input, "graph_is_res", id = "check",
                            attributes = Dict(:type => :checkbox, 
                           :style => "margin-top: 12px; margin-right: 5px")), 
                        Node(:div, "Data represents resistances instead of conductances", 
                             attributes = Dict(:style => "margin-top: 12px")))

	# Get file path 
	s1 = Scope()
	s1.dom = graph
	onimport(s1, JSExpr.@js function ()
				@var el = this.dom.querySelector("#input")
				el.onchange = (function ()
				   $(s1["filepath"])[] = el.files[0].path
				end)
			end)

	s2 = Scope()
	s2.dom = graph_is_res
	onimport(s2, JSExpr.@js function ()
				    @var el = this.dom.querySelector("#check")	
                    el.onchange = (function ()
                        $(s2["check"])[] = el.checked
                    end)
                end)

    s1, s2
end

function pairwise_input_ui()
    additional_input = Node(:div, tachyons_css, "Pairwise Mode Options") |> class"f4 lh-title"
    focal = vbox(Node(:div, "Select focal node locations from file: ", 
                      attributes = Dict(:style => "margin-top: 12px")),
                 Node(:input, id = "focal", attributes = Dict(:type => :file, 
                                                 :style => "margin-top: 12px")))
    pair = vbox(additional_input, 
                 focal)

    #=s = Scope()
    s.dom = pair
    onimport(s, JSExpr.@js function ()
                 @var el = this.dom.querySelector("#focal")
                 el.onchange = (function ()
                     $(s["focal"])[] = el.files[0].path
                    end)
             end)=#

    pair
end

