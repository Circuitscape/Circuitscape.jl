function parse_config(path::String)
    cfg = Dict{String, String}()
    f = open(path, "r")
    for i in EachLine(f)
        if first(i) == '['
            continue
        end
        idx = search(i, '=')
        var = rstrip(i[1:idx-1])
        val = strip(i[idx+1:end])
        cfg[var] = val
    end
    close(f)
    cfg
end
