# const fmt = "[{date} | {level} | {name}]: {msg}"
const fmt = x -> Dates.format(x, "yyyy-mm-dd HH:MM:SS") 
const ui_interface = Ref{Function}((x,y) -> nothing)

csinfo(msg) = (@info(string(fmt(Dates.now())) * " : " * msg); ui_interface[](msg, :info))
cswarn(msg) = (@warn(string(fmt(Dates.now())) * " : " * msg); ui_interace[](msg, :warn))

function update_logging!(cfg)

    log_level = cfg["log_level"]
    log_file = cfg["log_file"]

    if log_level in DEBUG
        Logging.LogLevel(Logging.Debug)
    elseif log_level in WARNING
        Logging.LogLevel(Logging.Warn)
    end
end
