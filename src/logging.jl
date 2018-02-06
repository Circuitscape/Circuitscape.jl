const fmt = "[{date} | {level} | {name}]: {msg}"
const logger = Memento.config("info", 
                        fmt = fmt)
csinfo(msg) = info(logger, msg)
cswarn(msg) = warn(logger, msg)

function update_logging!(cfg)

    log_level = cfg["log_level"]
    log_file = cfg["log_file"]

    if log_level in DEBUG
        setlevel!(logger, "debug")
    elseif log_level in WARNING
        setlevel!(logger, "warn")
    elseif log_level in CRITICAL
        setlevel!(logger, "critical")
    end
    if !(log_file in NONE)
        push!(logger, 
            DefaultHandler(log_file, DefaultFormatter(fmt)))
    end
end