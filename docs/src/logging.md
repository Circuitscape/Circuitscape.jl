# Logging Options

Circuitscape uses a custom `CSLogger` built on Julia's `AbstractLogger` interface.
The logger is configured automatically when you call `compute()`, based on settings
in your INI file.

## Suppressing Messages

Set `suppress_messages = True` in your INI file to suppress informational messages
during computation. Warnings will still be displayed.

## Logging to File

Set `log_file` in your INI file to write log messages to a file:

```
log_file = /path/to/logfile.log
```

When a log file is set, messages are written to both the console and the file.

## Disabling Log Output

To disable all informational log messages from Julia's side, use the built-in
logging system:

```julia
using Logging
Logging.disable_logging(Logging.Info)
```

To re-enable:

```julia
Logging.disable_logging(Logging.Debug)
```

## UI Callback

Circuitscape exposes a `ui_interface` callback for GUI integration. This is a
`Ref{Function}` that receives `(message, level)` for every log event, where
`level` is `:info` or `:warn`. This is used by downstream packages like Omniscape
to integrate Circuitscape's logging into their own UI.
