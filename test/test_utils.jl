using Circuitscape
using Test
using DelimitedFiles

import Circuitscape: parse_config, CSConfig

"""Clean the test output directory."""
function clean_output()
    outdir = joinpath(@__DIR__, "output")
    rm(outdir, force=true, recursive=true)
    mkpath(outdir)
end

"""
    compute_with(str; solver, precision, parallel)

Run a Circuitscape job with overridden settings.
"""
function compute_with(str::String;
                      solver::String = "",
                      precision::String = "",
                      parallel::Bool = false)
    cfg = parse_config(str)
    d = Dict{String,String}(cfg)
    solver != "" && (d["solver"] = solver)
    precision != "" && (d["precision"] = precision)
    d["parallelize"] = parallel ? "true" : "false"
    compute(d)
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

"""
    runtests(; solver, precision, parallel, label)

Run the full Circuitscape test suite with the given settings.
"""
function runtests(; solver::String = "", precision::String = "",
                    parallel::Bool = false, label::String = "")
    if label == ""
        parts = String[]
        push!(parts, solver == "" ? "CG+AMG" : uppercase(solver))
        push!(parts, parallel ? "parallel" : "sequential")
        push!(parts, precision == "single" ? "Float32" : "Float64")
        label = join(parts, ", ")
    end

    is_single = precision == "single"
    tol = is_single ? 1e-4 : 1e-6
    f(str) = compute_with(str; solver, precision, parallel)

    @testset "$label" begin
        @testset "Network Pairwise" begin
            for i = 1:3
                r = f("input/network/sgNetworkVerify$(i).ini")
                x = readdlm("output_verify/sgNetworkVerify$(i)_resistances.out")
                valx = x[2:end, 2:end]
                valr = r[2:end, 2:end]
                @test check_resistances(valx, valr, tol, label="sgNetworkVerify$(i)")
                pts_x = x[2:end,1]
                pts_r = r[2:end,1]
                @test pts_x .+ 1 == pts_r
                compare_all_output("sgNetworkVerify$(i)", is_single)
            end
        end

        @testset "Network Advanced" begin
            for i = 1:3
                r = f("input/network/mgNetworkVerify$(i).ini")
                x = readdlm("output_verify/mgNetworkVerify$(i)_voltages.txt")
                @. x[:,1] = x[:,1] + 1
                @test check_resistances(x, r, tol, label="mgNetworkVerify$(i)")
                compare_all_output("mgNetworkVerify$(i)", is_single)
            end
        end

        @testset "Raster Pairwise" begin
            for i = 1:17
                if i == 16 && Sys.WORD_SIZE == 32
                    continue
                end

                r = f("input/raster/pairwise/$i/sgVerify$(i).ini")
                x = readdlm("output_verify/sgVerify$(i)_resistances.out")
                _x = readdlm("output/sgVerify$(i)_resistances.out")
                @test check_resistances(_x, r, tol, label="sgVerify$(i) (written)")
                @test check_resistances(x, r, tol, label="sgVerify$(i) (verify)")
                compare_all_output("sgVerify$(i)", is_single)
            end
        end

        @testset "Raster Advanced" begin
            for i in 1:6
                r = f("input/raster/advanced/$i/mgVerify$(i).ini")
                compare_all_output("mgVerify$(i)")
            end
        end

        @testset "Raster One to All" begin
            for i in 1:13
                r = f("input/raster/one_to_all/$i/oneToAllVerify$(i).ini")
                x = readdlm("output_verify/oneToAllVerify$(i)_resistances.out")
                @test check_resistances(x, r, tol, label="oneToAllVerify$(i)")
                compare_all_output("oneToAllVerify$(i)", is_single)
            end
        end

        @testset "Raster All to One" begin
            for i in 1:12
                r = f("input/raster/all_to_one/$i/allToOneVerify$(i).ini")
                x = readdlm("output_verify/allToOneVerify$(i)_resistances.out")
                @test check_resistances(x, r, tol, label="allToOneVerify$(i)")
                compare_all_output("allToOneVerify$(i)", is_single)
            end
        end
    end
end

"""Check resistance matrices element-wise and report which entries differ."""
function check_resistances(x, r, tol; label="")
    nfail = 0
    for j in axes(x, 2), i in axes(x, 1)
        if abs(x[i,j] - r[i,j]) > sqrt(tol)
            nfail += 1
            if nfail <= 10
                # Print row/col headers (focal point IDs) if available
                ri = size(r,1) >= i ? r[i,1] : i
                rj = size(r,2) >= j ? r[1,j] : j
                @warn "$label MISMATCH [$i,$j] (points $ri,$rj): expected=$(x[i,j]) got=$(r[i,j]) diff=$(abs(x[i,j]-r[i,j]))"
            end
        end
    end
    if nfail > 10
        @warn "$label $nfail total entries differ (showing first 10)"
    elseif nfail > 0
        @warn "$label $nfail entries differ"
    end
    nfail == 0
end

function compare_all_output(str, is_single = false)
    gen_list, list_to_comp = generate_lists(str)
    tol = is_single ?  1e-4 : 1e-6

    for f in gen_list
        !occursin("_", f) && continue
        occursin("resistances", f) && continue

        if endswith(f, "asc")
            r = read_aagrid("output/$f")
            x = get_comp(list_to_comp, f)
            @test compare_aagrid(r, x, tol)
        elseif occursin("Network", f)
            if occursin("branch", f)
                r = read_branch_currents("output/$f")
                x = !startswith(f, "mg") ? get_network_comp(list_to_comp, f) : readdlm("output_verify/$f")
                @test compare_branch(r, x, tol)
            else
                r = read_node_currents("output/$f")
                x = !startswith(f, "mg") ? get_network_comp(list_to_comp, f) : readdlm("output_verify/$f")
                @test compare_node(r, x, tol)
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
