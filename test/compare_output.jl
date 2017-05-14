function compare_all_output(str)

    gen_list, list_to_comp = generate_lists(str)

    for f in gen_list
        info("Testing $f")
        if endswith(f, "asc")
            r = read_aagrid("output/$f")
            x = get_comp(list_to_comp, f)
            @test compare_aagrid(r, x)
            info("Test $f passed")
        end
    end
            
end

list_of_files(str, pref) = readdir(pref) |> y -> filter(x -> startswith(x, str), y)
generate_lists(str) = list_of_files(str, "output/"), list_of_files(str, "output_verify/")

read_aagrid(file) = readdlm(file, skipstart = 0) # Will change to 6 

compare_aagrid{T}(r::Matrix{T}, x::Matrix{T}) = sumabs2(x - r) < 1e-6

function get_comp(list_to_comp, f)
    outfile = ""
    if f in list_to_comp
        outfile = "output_verify/$f"
    end
    readdlm(outfile; skipstart = 6)
end

