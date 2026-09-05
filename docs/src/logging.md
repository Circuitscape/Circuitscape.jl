# Logging Options

Circuitscape uses a custom `CSLogger` built on Julia's `AbstractLogger` interface.
The logger is installed as the global logger when you call `Circuitscape.compute()`, based on
the `log_level`, `log_file` and `suppress_messages` settings in your INI file.
Every message is prefixed with a timestamp.

## Log Level

`log_level` accepts `DEBUG`, `INFO` (the default), `WARNING`, `ERROR` and
`CRITICAL`, in any case. `WARN` is also accepted, and `CRITICAL` is treated as
`ERROR`, as in Circuitscape 4. Any other value is an error.

```
log_level = DEBUG
```

At `DEBUG`, every pair solve is logged (`Solving pair k of n`) and, when the run
finishes, a timing table produced with
[TimerOutputs](https://github.com/KristofferC/TimerOutputs.jl) is printed. It
breaks the run down into loading the raster, constructing the graph and the
preconditioner or factorization, `solve linear system` and `postprocess` per
task, and writing cumulative maps. This is the place to look when deciding
between solvers or checking that threads are being used.

## Suppressing Messages

Set `suppress_messages = True` in your INI file to suppress informational messages
on the console. Warnings and errors are still displayed, and a log file, if set,
still receives every message.

## Logging to File

Set `log_file` in your INI file to write log messages to a file:

```
log_file = /path/to/logfile.log
```

When a log file is set, messages are written to both the console and the file.
The file is created (or truncated) at the start of each run.

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
