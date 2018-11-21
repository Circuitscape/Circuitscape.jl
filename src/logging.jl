const fmt = "[{date} | {level} | {name}]: {msg}"
const ui_interface = Ref{Function}((x,y) -> nothing)

csinfo(msg) = (@info(string(Dates.now()) * " : " * msg); ui_interface[](msg, :info))
cswarn(msg) = (@warn(string(Dates.now()) * " : " * msg); ui_interace[](msg, :warn))

function update_logging!(cfg)

    log_level = cfg["log_level"]
    log_file = cfg["log_file"]

    if log_level in DEBUG
        Logging.LogLevel(Logging.Debug)
    elseif log_level in WARNING
        Logging.LogLevel(Logging.Warn)
    end
end
