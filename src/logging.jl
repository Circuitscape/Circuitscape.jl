const ui_interface = Ref{Function}((x,y) -> nothing)

struct CSLogger <: AbstractLogger
    console_logger::ConsoleLogger
    file_logger::Union{SimpleLogger, Nothing}
    suppress_messages::Bool
end

function CSLogger(; level::Logging.LogLevel = Logging.Info, suppress_messages::Bool = false, file_logger::Union{SimpleLogger, Nothing} = nothing)
    CSLogger(ConsoleLogger(stderr, level), file_logger, suppress_messages)
end

Logging.min_enabled_level(logger::CSLogger) = Logging.min_enabled_level(logger.console_logger)
Logging.shouldlog(logger::CSLogger, level, _module, group, id) = true
Logging.catch_exceptions(logger::CSLogger) = false

function Logging.handle_message(logger::CSLogger, level, message, _module, group, id, filepath, line; kwargs...)
    ts = Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS")
    msg = "$ts : $message"

    # Forward to UI callback
    ui_level = level >= Logging.Warn ? :warn : :info
    ui_interface[](msg, ui_level)

    # Log to console unless suppressed
    if !logger.suppress_messages || level >= Logging.Warn
        Logging.handle_message(logger.console_logger, level, msg, _module, group, id, filepath, line; kwargs...)
    end

    # Log to file if configured
    if logger.file_logger !== nothing
        with_logger(logger.file_logger) do
            if level >= Logging.Warn
                @warn(msg)
            else
                @info(msg)
            end
        end
    end
    nothing
end

function update_logging!(cfg::CSConfig)
    file_logger = if cfg.log_file != ""
        SimpleLogger(open(cfg.log_file, "w+"))
    else
        nothing
    end

    logger = CSLogger(
        level = cfg.log_level,
        suppress_messages = cfg.suppress_messages,
        file_logger = file_logger
    )
    global_logger(logger)

    if cfg.log_file != ""
        @info("Logs will recorded to file: $(cfg.log_file)")
    end
end
