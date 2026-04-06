# const fmt = "[{date} | {level} | {name}]: {msg}"
const fmt = x -> Dates.format(x, "yyyy-mm-dd HH:MM:SS")
const ui_interface = Ref{Function}((x,y) -> nothing)

mutable struct LogState
    log_to_file::Bool
    file_logger::SimpleLogger
end
const logging = LogState(false, SimpleLogger())

function csinfo(msg, suppress_messages::Bool = false)
    log_to_file = logging.log_to_file
    file_logger = logging.file_logger
    msg = string(fmt(Dates.now())) * " : " * msg
    !suppress_messages && @info(msg)
    ui_interface[](msg, :info)
    if log_to_file
        with_logger(file_logger) do
            @info(msg)
        end
    end
end

function cswarn(msg)
    log_to_file = logging.log_to_file
    file_logger = logging.file_logger
    msg = string(fmt(Dates.now())) * " : " * msg
    @warn(msg)
    ui_interface[](msg, :warn)
    if log_to_file
        with_logger(file_logger) do
            @warn(msg)
        end
    end
end

function update_logging!(cfg::CSConfig)

    log_level = cfg.log_level
    log_file = cfg.log_file

    if log_level == ll_debug
        Logging.LogLevel(Logging.Debug)
    elseif log_level == ll_warning
        Logging.LogLevel(Logging.Warn)
    end

    if log_file != ""
        logging.file_logger = SimpleLogger(open(log_file, "w+"))
        logging.log_to_file = true
        csinfo("Logs will recorded to file: $log_file", cfg.suppress_messages)
    else
        logging.file_logger = SimpleLogger()
        logging.log_to_file = false
    end
end
