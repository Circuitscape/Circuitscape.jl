# const fmt = "[{date} | {level} | {name}]: {msg}"
const fmt = x -> Dates.format(x, "yyyy-mm-dd HH:MM:SS") 
const ui_interface = Ref{Function}((x,y) -> nothing)
const logging = Dict{String,Any}()
logging["log_to_file"] = false
logging["file_logger"] = SimpleLogger()

function csinfo(msg) 
    log_to_file = logging["log_to_file"]
    file_logger = logging["file_logger"]
    msg = string(fmt(Dates.now())) * " : " * msg 
    @info(msg)
    ui_interface[](msg, :info)
    if log_to_file
        with_logger(file_logger) do 
            @info(msg)
        end
    end
end

function cswarn(msg)
    log_to_file = logging["log_to_file"]
    file_logger = logging["file_logger"]
    msg = string(fmt(Dates.now())) * " : " * msg 
    @warn(msg)
    ui_interface[](msg, :warn)
    if log_to_file
        with_logger(file_logger) do 
            @warn(msg)
        end
    end
end

function update_logging!(cfg)

    log_level = cfg["log_level"]
    log_file = cfg["log_file"]

    if log_level in DEBUG
        Logging.LogLevel(Logging.Debug)
    elseif log_level in WARNING
        Logging.LogLevel(Logging.Warn)
    end

    if log_file != "None"
        logging["file_logger"] = SimpleLogger(open(log_file, "w+"))
        logging["log_to_file"] = true
        csinfo("Logs will recorded to file: log_file")
    else
        logging["file_logger"] = SimpleLogger()
        logging["log_to_file"] = false
    end
end
