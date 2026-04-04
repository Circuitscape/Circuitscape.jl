# Logging Options

## Disabling Log Output

To disable all informational log messages, use Julia's built-in logging system:

```julia
using Logging
Logging.disable_logging(Logging.Info)
```

This disables logging at a global level. To re-enable logging:

```julia
Logging.disable_logging(Logging.Debug)
```

## Logging to File

Set `log_file` in your INI file to a file path to write log messages to a file:

```
log_file = /path/to/logfile.log
```

## Suppressing Messages

Set `suppress_messages = True` in your INI file to suppress informational messages during computation. Warnings will still be displayed.

## Log Level

Set `log_level` in your INI file to control verbosity:

- **`DEBUG`** — Most verbose. Includes solver residual norms for verifying solution accuracy, in addition to all lower levels.
- **`INFO`** (default) — Reports solver selection, timing, progress through pair solves, and job completion.
- **`WARNING`** — Only warnings (e.g., precision overrides when using direct solvers with single precision).
- **`CRITICAL`** — Only critical errors.
