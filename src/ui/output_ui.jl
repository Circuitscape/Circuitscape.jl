function output_ui()
    
    # Title of section
    title = Node(:div, tachyons_css, "Output Options") |> class"f4 lh-title"

    base_name = vbox(Node(:div, "Base output file name: ",
                      attributes = Dict(:style => "margin-top: 12px")) |> class"b",
                 Node(:input, id = "out", attributes = Dict(:type => :text, 
                                                :style => "margin-top: 12px")))

    options = vbox(hbox(Node(:input, id = "cur", attributes = Dict(:type => :checkbox, 
                           :style => "margin-top: 12px; margin-right: 5px")), 
                        Node(:div, "Current maps", 
                             attributes = Dict(:style => "margin-top: 12px"))),
                  hbox(Node(:input, id = "volt", attributes = Dict(:type => :checkbox,
                           :style => "margin-top: 12px; margin-right: 5px")), 
                        Node(:div, "Voltage maps", 
                             attributes = Dict(:style => "margin-top: 12px"))))

    output_ui = vbox(title,
                     base_name, 
                     options)

    s = Scope()

	s.dom = output_ui
	onimport(s, JSExpr.@js function ()
				    @var el1 = this.dom.querySelector("#cur")	
				    @var el2 = this.dom.querySelector("#volt")	
                    @var el3 = this.dom.querySelector("#out")
                    el1.onchange = (function ()
                        $(s["cur"])[] = el1.checked
                        $(s["volt"])[] = el2.checked 
                    end)
                    el2.onchange = (function ()
                        $(s["cur"])[] = el1.checked
                        $(s["volt"])[] = el2.checked 
                    end)
                    el3.onchange = (function ()
                        $(s["out"])[] = el3.value
                    end)
                end)
    s
end
