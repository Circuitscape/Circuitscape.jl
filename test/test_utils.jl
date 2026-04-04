using Circuitscape
using Test
using DelimitedFiles

import Circuitscape: parse_config, _compute, SINGLE, CHOLMOD, PARDISO, TRUELIST, cswarn

function compute_cholmod(str, batch_size = 5)
    cfg = parse_config(str)
    T = cfg["precision"] in SINGLE ? Float32 : Float64
    V = cfg["use_64bit_indexing"] in TRUELIST ? Int64 : Int32
    if T == Float32
        cswarn("Cholmod supports only double precision. Changing precision to double.")
        T = Float64
    end
    cfg["solver"] = "cholmod"
    cfg["cholmod_batch_size"] = string(batch_size)
    _compute(T, V, cfg)
end

function compute_pardiso(str, batch_size = 5)
    cfg = parse_config(str)
    T = cfg["precision"] in SINGLE ? Float32 : Float64
    V = cfg["use_64bit_indexing"] in TRUELIST ? Int64 : Int32
    if T == Float32
        cswarn("Pardiso supports only double precision. Changing precision to double.")
        T = Float64
    end
    cfg["solver"] = "pardiso"
    cfg["cholmod_batch_size"] = string(batch_size)
    _compute(T, V, cfg)
end

function compute_single(str)
    cfg = parse_config(str)
    cfg["precision"] = "single"
    V = cfg["use_64bit_indexing"] in TRUELIST ? Int64 : Int32
    T = Float32
    if cfg["solver"] in CHOLMOD || cfg["solver"] in PARDISO
        cswarn("Cholmod and Pardiso support only double precision. Changing precision to double.")
        T = Float64
    end
    _compute(T, V, cfg)
end

function compute_parallel(str, n_processes = 2)
    cfg = parse_config(str)
    cfg["parallelize"] = "true"
    cfg["max_parallel"] = "$(n_processes)"
    compute(cfg)
end

function model_problem(T, s)
    cellmap = ones(T, s, s)
    nodemap = reshape(1:lastindex(cellmap), s, s) |> collect
    G = Circuitscape.construct_graph(cellmap, nodemap, true, true)
    Circuitscape.laplacian(G)
end
model_problem(s::Integer) = model_problem(Float64, s)

function test_problem(str)
    base_path = joinpath(dirname(pathof(Circuitscape)), "..", "test", "input")
    str2 = replace(str, ".ini" => "")
    if occursin("sgVerify", str)
        config_path = joinpath(base_path, "raster", "pairwise", replace(str2, "sgVerify" => ""), str)
    elseif occursin("mgVerify", str)
        config_path = joinpath(base_path, "raster", "advanced", replace(str2, "mgVerify" => ""), str)
    elseif occursin("oneToAll", str)
        config_path = joinpath(base_path, "raster", "one_to_all", replace(str2, "oneToAllVerify" => ""), str)
    elseif occursin("allToOne", str)
        config_path = joinpath(base_path, "raster", "all_to_one", replace(str2, "allToOneVerify" => ""), str)
    elseif occursin("sgNetworkVerify", str)
        config_path = joinpath(base_path, "network", str)
    elseif occursin("mgNetworkVerify", str)
        config_path = joinpath(base_path, "network", str)
    end
end

function runtests(f = compute)

    str = if f == compute_single
            "Single"
          else
            "Double"
          end

    is_single = false
    tol = 1e-6
    str == "Single" && (is_single = true; tol = 1e-4)

    @testset "$str Precision Tests" begin
        @testset "Network Pairwise" begin
            for i = 1:3
                @info("Testing sgNetworkVerify$i")
                r = f("input/network/sgNetworkVerify$(i).ini")
                x = readdlm("output_verify/sgNetworkVerify$(i)_resistances.out")
                valx = x[2:end, 2:end]
                valr = r[2:end, 2:end]
                @test sum(abs2, valx - valr) < tol
                pts_x = x[2:end,1]
                pts_r = r[2:end,1]
                @test pts_x .+ 1 == pts_r
                compare_all_output("sgNetworkVerify$(i)", is_single)
                @info("Test sgNetworkVerify$i passed")
            end
        end

        @testset "Network Advanced" begin
            for i = 1:3
                @info("Testing mgNetworkVerify$i")
                r = f("input/network/mgNetworkVerify$(i).ini")
                x = readdlm("output_verify/mgNetworkVerify$(i)_voltages.txt")
                @. x[:,1] = x[:,1] + 1
                @test sum(abs2, x - r) < tol
                compare_all_output("mgNetworkVerify$(i)", is_single)
                @info("Test mgNetworkVerify$i passed")
            end
        end

        @testset "Raster Pairwise" begin
            for i = 1:17
                if i == 16 && Sys.WORD_SIZE == 32
                    continue
                end

                @info("Testing sgVerify$i")
                r = f("input/raster/pairwise/$i/sgVerify$(i).ini")
                x = readdlm("output_verify/sgVerify$(i)_resistances.out")
                _x = readdlm("output/sgVerify$(i)_resistances.out")
                @test sum(abs2, _x - r) < tol
                @test sum(abs2, x - r) < tol
                compare_all_output("sgVerify$(i)", is_single)
                @info("Test sgVerify$i passed")
            end
        end

        @testset "Raster Advanced" begin
            for i in 1:6
                @info("Testing mgVerify$i")
                r = f("input/raster/advanced/$i/mgVerify$(i).ini")
                compare_all_output("mgVerify$(i)")
                @info("Test mgVerify$i passed")
            end
        end

        @testset "Raster One to All" begin
            for i in 1:13
                @info("Testing oneToAllVerify$i")
                r = f("input/raster/one_to_all/$i/oneToAllVerify$(i).ini")
                x = readdlm("output_verify/oneToAllVerify$(i)_resistances.out")
                @test sum(abs2, x - r) < tol
                compare_all_output("oneToAllVerify$(i)", is_single)
                @info("Test oneToAllVerify$i passed")
            end
        end

        @testset "Raster All to One" begin
            for i in 1:12
                @info("Testing allToOneVerify$i")
                r = f("input/raster/all_to_one/$i/allToOneVerify$(i).ini")
                x = readdlm("output_verify/allToOneVerify$(i)_resistances.out")
                @test sum(abs2, x - r) < tol
                compare_all_output("allToOneVerify$(i)", is_single)
                @info("Test allToOneVerify$i passed")
            end
        end
    end
end

function compare_all_output(str, is_single = false)

    gen_list, list_to_comp = generate_lists(str)
    tol = is_single ?  1e-4 : 1e-6

    for f in gen_list
        !occursin("_", f) && continue
        occursin("resistances", f) && continue

        @info("Testing $f")

        if endswith(f, "asc")
            r = read_aagrid("output/$f")
            x = get_comp(list_to_comp, f)
            @test compare_aagrid(r, x, tol)
            @info("Test $f passed")

        elseif occursin("Network", f)
            if occursin("branch", f)
                r = read_branch_currents("output/$f")
                x = !startswith(f, "mg") ? get_network_comp(list_to_comp, f) : readdlm("output_verify/$f")
                @test compare_branch(r, x, tol)
                @info("Test $f passed")
            else
                r = read_node_currents("output/$f")
                x = !startswith(f, "mg") ? get_network_comp(list_to_comp, f) : readdlm("output_verify/$f")
                @test compare_node(r, x, tol)
                @info("Test $f passed")
            end
        end
    end

end

list_of_files(str, pref) = readdir(pref) |> y -> filter(x -> startswith(x, "$(str)_"), y)
generate_lists(str) = list_of_files(str, "output/"), list_of_files(str, "output_verify/")
read_branch_currents(str) = readdlm(str)
read_node_currents(str) = readdlm(str)

read_aagrid(file) = readdlm(file, skipstart = 6)

compare_aagrid(r::Matrix{T}, x::Matrix{T}, tol = 1e-6) where T = sum(abs2, x - r) < tol

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
        if all(isnumeric, s[i])
            f = replace(f, "_$(s[i])"=>"_$(parse(Int, s[i]) - 1 |> string)", count=1)
        end
    end
    @assert isfile("output_verify/$f")
    readdlm("output_verify/$f")
end

function compare_branch(r, x, tol = 1e-6)
    @. x[:,1] = x[:,1] + 1
    @. x[:,2] = x[:,2] + 1
    sum(abs2, sortslices(r, dims=1) - sortslices(x, dims=1)) < tol
end

function compare_node(r, x, tol = 1e-6)
    @. x[:,1] = x[:,1] + 1
    sum(abs2, sortslices(r, dims=1) - sortslices(x, dims=1)) < tol
end

function change_to_test_path()
    origpath = pwd()
    pkgpath = Base.pathof(Circuitscape)
    cd(joinpath(dirname(pkgpath), "..", "test"))
    origpath
end

function runtest(str)
    origpath = change_to_test_path()
    prob = test_problem(str)
    compute(prob)
    cd(origpath)
end
