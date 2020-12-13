### The Log Window

To disable all logging, write the following code: 
```julia
using Logging
Logging.disable_logging(Logging.Info)
```
This disables logging on a global level. If you 
want to enable logging again, run

```julia
Base.CoreLogging._min_enabled_level = Logging.LogLevel(0)
```