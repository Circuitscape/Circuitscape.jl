function compare_all_output(str)

    gen_list, list_to_comp = generate_lists(str)
    @show gen_list

    for f in gen_list
        info("Testing $f")

        # Raster output files
        if endswith(f, "asc")
            r = read_aagrid("output/$f")
            x = get_comp(list_to_comp, f)
            @test compare_aagrid(r, x)
            info("Test $f passed")

        # Network output files
        elseif contains(f, "Network")

            # Branch currents
            if contains(f, "branch")
                r = read_branch_currents("output/$f")
                x = !startswith(f, "mg") ? get_network_comp(list_to_comp, f) : readdlm("output_verify/$f")
                @test compare_branch(r, x)
                info("Test $f passed")

            # Node currents
            else
                r = read_node_currents("output/$f")
                x = !startswith(f, "mg") ? get_network_comp(list_to_comp, f) : readdlm("output_verify/$f")
                @test compare_node(r, x)
                info("Test $f passed")
            end
        end
    end
            
end

list_of_files(str, pref) = readdir(pref) |> y -> filter(x -> startswith(x, str), y)
generate_lists(str) = list_of_files(str, "output/"), list_of_files(str, "output_verify/")
read_branch_currents(str) = readdlm(str)
read_node_currents(str) = readdlm(str)

read_aagrid(file) = readdlm(file, skipstart = 0) # Will change to 6 

compare_aagrid{T}(r::Matrix{T}, x::Matrix{T}) = sumabs2(x - r) < 1e-6

function get_comp(list_to_comp, f)
    outfile = ""
    if f in list_to_comp
        outfile = "output_verify/$f"
    end
    readdlm(outfile; skipstart = 6)
end

function get_network_comp(list_to_comp, f)
    s = split(f, ['_', '.'])
    for i = 1:size(s, 1)
        if isnumber(s[i])
            f = replace(f, "_$(s[i])", "_$(parse(Int, s[i]) - 1 |> string)", 1)
        end
    end
    @assert isfile("output_verify/$f")
    readdlm("output_verify/$f")
end

function compare_branch(r, x)
    x[:,1] = x[:,1] + 1
    x[:,2] = x[:,2] + 1
    sumabs2(sortrows(r) - sortrows(x)) < 1e-6
end

function compare_node(r, x)
    x[:,1] = x[:,1] + 1
    sumabs2(sortrows(r) - sortrows(x)) < 1e-6
end
