@enum DataType dt_raster dt_network
@enum Scenario sc_pairwise sc_advanced sc_onetoall sc_alltoone
@enum SolverType st_cg_amg st_cholmod st_pardiso st_accelerate st_amgx
@enum Precision pr_single pr_double
@enum RemovePolicy rp_keepall rp_rmvsrc rp_rmvgnd rp_rmvall

"""
    CSConfig

Typed, immutable configuration for one Circuitscape run: one field per INI key,
with the INI's string values parsed into enums, booleans and numbers. Built by
`parse_config(path)` or `CSConfig(dict::Dict{String,String})`; every parser
rejects values it does not recognise rather than falling back to a default.
Keys and option values are case-insensitive (see [`normalize_config`](@ref));
file paths keep their case.
`Dict{String,String}(cfg)` converts back to INI strings, and
`CSConfig(cfg; field = value, ...)` returns a copy with fields replaced.

Field names match the INI keys documented in the options reference, which
also gives their defaults. The field defaults are the one statement of the
default configuration; [`init_config`](@ref) is them written as INI strings.
`residual_tolerance = 0.0` means automatic (see [`residual_tolerance`](@ref)),
and `scenario` has no usable default and must be set.
"""
Base.@kwdef struct CSConfig
    version::String = "unknown"
    data_type::DataType = dt_raster
    scenario::Scenario = sc_pairwise
    habitat_file::String = ""
    habitat_map_is_resistances::Bool = true
    connect_four_neighbors_only::Bool = false
    connect_using_avg_resistances::Bool = false
    use_polygons::Bool = false
    polygon_file::String = ""
    source_file::String = ""
    ground_file::String = ""
    ground_file_is_resistances::Bool = true
    use_unit_currents::Bool = false
    use_direct_grounds::Bool = false
    remove_src_or_gnd::RemovePolicy = rp_keepall
    use_mask::Bool = false
    mask_file::String = ""
    solver::SolverType = st_cg_amg
    parallelize::Bool = false
    precision::Precision = pr_double
    use_64bit_indexing::Bool = false   # false = 32-bit indices when the problem fits
    cholmod_batch_size::Int = 32
    residual_tolerance::Float64 = 0.0   # 0 = automatic: 1e-4 double, 1e-3 single
    low_memory_mode::Bool = false
    preemptive_memory_release::Bool = false
    use_variable_source_strengths::Bool = false
    variable_source_file::String = ""
    use_included_pairs::Bool = false
    included_pairs_file::String = ""
    point_file::String = ""
    use_reclass_table::Bool = false
    reclass_file::String = ""
    output_file::String = ""
    write_cur_maps::Bool = false
    write_volt_maps::Bool = false
    write_cum_cur_map_only::Bool = false
    write_max_cur_maps::Bool = false
    set_null_currents_to_nodata::Bool = false
    set_null_voltages_to_nodata::Bool = false
    set_focal_node_currents_to_zero::Bool = false
    compress_grids::Bool = false
    log_transform_maps::Bool = false
    write_as_tif::Bool = false
    log_file::String = ""
    log_level::Logging.LogLevel = Logging.Info
    suppress_messages::Bool = false
end

# The INI spellings each boolean and enumerated value is read from, in any
# case; the first is the one written. Every parser rejects values it does not
# recognise.
# Falling back to a default silently turned typos into different runs:
# `solver = cholmodd` ran cg+amg, `scenario = pairwsie` ran pairwise,
# `write_cur_maps = Ture` wrote nothing (issues #246, #341, #436).
const INI_NAMES = [
    Bool => [true => ["True", "1"], false => ["False", "0"]],
    DataType => [dt_raster => ["raster"], dt_network => ["network"]],
    Scenario => [sc_pairwise => ["pairwise"],
                 sc_advanced => ["advanced"],
                 sc_onetoall => ["one-to-all", "one_to_all"],
                 sc_alltoone => ["all-to-one", "all_to_one"]],
    SolverType => [st_cg_amg => ["cg+amg", "amg+cg"],
                   st_cholmod => ["cholmod", "cholesky", "cholfact"],
                   st_pardiso => ["pardiso", "mklpardiso"],
                   st_accelerate => ["accelerate", "apple_accelerate"],
                   st_amgx => ["amgx"]],
    Precision => [pr_single => ["single"], pr_double => ["double"]],
    # Circuitscape 4 wrote "not entered" for the default policy
    RemovePolicy => [rp_keepall => ["keepall", "not entered"], rp_rmvsrc => ["rmvsrc"],
                     rp_rmvgnd => ["rmvgnd"], rp_rmvall => ["rmvall"]],
    # Circuitscape 4 took Python's level names; keep reading them.
    Logging.LogLevel => [Logging.Debug => ["DEBUG"], Logging.Info => ["INFO"],
                         Logging.Warn => ["WARNING", "WARN"],
                         Logging.Error => ["ERROR", "CRITICAL"]],
]
const Named = Union{Bool, Enum, Logging.LogLevel}
# A method per type, so that parsing a field does no table lookup at run time.
for (T, names) in INI_NAMES
    @eval _names(::Type{$T}) = $names
end

# Equality ignoring ASCII case, without allocating; every spelling is ASCII.
_lower(c::UInt8) = UInt8('A') <= c <= UInt8('Z') ? c + 0x20 : c
_same(a, b) = ncodeunits(a) == ncodeunits(b) &&
    all(i -> _lower(codeunit(a, i)) == _lower(codeunit(b, i)), 1:ncodeunits(a))

# `_from_ini(T, key, s)` reads the INI string `s` of `key` as a `T`;
# `_to_ini(v)` writes `v` back.
function _from_ini(::Type{T}, key, s) where {T<:Named}
    for (v, names) in _names(T)
        any(n -> _same(s, n), names) && return v
    end
    valid = join((first(names) for (_, names) in _names(T)), ", ")
    # init_config leaves `scenario` as "not entered"
    _same(s, "not entered") && throw(ArgumentError("$key is not set; expected one of: $valid"))
    throw(ArgumentError("$key = \"$s\" is not recognised; expected one of: $valid"))
end
_to_ini(v::T) where {T<:Named} = first(only(names for (x, names) in _names(T) if x == v))

# Paths and free text keep their case. "None" is how Circuitscape 4 (Python)
# wrote an unset path.
_from_ini(::Type{String}, key, s) = _same(s, "None") ? "" : String(s)
_to_ini(v::String) = v == "" ? "None" : v

_from_ini(::Type{Int}, key, s) =
    @something(tryparse(Int, s), throw(ArgumentError("$key = \"$s\" is not an integer")))
_to_ini(v::Int) = string(v)

# residual_tolerance, the only real-valued option: 0 means automatic.
_from_ini(::Type{Float64}, key, s) =
    any(x -> _same(s, x), ("auto", "none", "")) ? 0.0 :
    @something(tryparse(Float64, s), throw(ArgumentError("$key = \"$s\" is not a number or auto")))
_to_ini(v::Float64) = v > 0 ? string(v) : "auto"

# The sections of the INI file `write_config` writes, and the keys in each, in
# the order Circuitscape 4 wrote them. Sections mean nothing when reading.
# Every CSConfig field appears exactly once.
const INI_SECTIONS = [
    "Circuitscape Mode" => [:data_type, :scenario],
    "Version" => [:version],
    "Habitat raster or graph" => [:habitat_file, :habitat_map_is_resistances],
    "Options for reclassification of habitat data" => [:use_reclass_table, :reclass_file],
    "Connection Scheme for raster habitat data" =>
        [:connect_four_neighbors_only, :connect_using_avg_resistances],
    "Short circuit regions (aka polygons)" => [:use_polygons, :polygon_file],
    "Options for advanced mode" =>
        [:ground_file_is_resistances, :source_file, :remove_src_or_gnd, :ground_file,
         :use_unit_currents, :use_direct_grounds],
    "Mask file" => [:use_mask, :mask_file],
    "Options for one-to-all and all-to-one modes" =>
        [:use_variable_source_strengths, :variable_source_file],
    "Options for pairwise and one-to-all and all-to-one modes" =>
        [:included_pairs_file, :use_included_pairs, :point_file],
    "Calculation options" =>
        [:solver, :residual_tolerance, :precision, :parallelize, :use_64bit_indexing,
         :cholmod_batch_size, :low_memory_mode, :preemptive_memory_release],
    "Output options" =>
        [:write_cum_cur_map_only, :log_transform_maps, :output_file, :write_max_cur_maps,
         :write_volt_maps, :set_null_currents_to_nodata, :set_null_voltages_to_nodata,
         :compress_grids, :write_cur_maps, :set_focal_node_currents_to_zero, :write_as_tif],
    "Logging Options" => [:log_level, :log_file, :suppress_messages],
]

# Where `init_config` differs from the CSConfig defaults: the INI Builder's
# placeholders, which `validate` treats as unset, and Circuitscape 4 keys that
# are accepted, so that old INI files load without warnings, and do nothing.
const INI_PLACEHOLDERS = Dict(
    "scenario" => "not entered",
    "habitat_file" => "(Browse for a resistance file)",
    "polygon_file" => "(Browse for a short-circuit region file)",
    "source_file" => "(Browse for a current source file)",
    "ground_file" => "(Browse for a ground point file)",
    "included_pairs_file" => "(Browse for a file with pairs to include or exclude)",
    "point_file" => "(Browse for file with locations of focal points or regions)",
    "reclass_file" => "(Browse for file with reclassification data)",
    "output_file" => "(Choose a base name for output files)",
    "print_timings" => "False",
    "print_rusages" => "False",
    "screenprint_log" => "False",
    "profiler_log_file" => "None",
)

# An unset path: empty (as `None` is read), `None` in a CSConfig built
# directly, or an INI Builder placeholder; every placeholder above starts
# with "(" (checked in the tests).
_is_unset(path) = isempty(path) || _same(path, "None") || startswith(path, "(")

Base.Dict{String,String}(cfg::CSConfig) =
    Dict{String,String}(String(f) => _to_ini(getfield(cfg, f)) for f in fieldnames(CSConfig))

const INI_DEFAULTS = merge(Dict{String,String}(CSConfig()), INI_PLACEHOLDERS)

# INI keys are matched in any case, ignoring surrounding whitespace.
_ini_key(k) = lowercase(strip(k))

"""
    normalize_config(dict) -> Dict{String,String}

`dict` with every key lowercased and every key and value stripped of
surrounding whitespace, so that `Solver = cholmod` sets `solver`. A key given
twice in different case with different values is an `ArgumentError`. Values
keep their case here: option values are matched in any case when they are
parsed, and paths and free text (the `String` fields of [`CSConfig`](@ref))
keep theirs.
"""
function normalize_config(dict)
    d = Dict{String,String}()
    for (k, v) in dict
        key = _ini_key(k)
        val = String(strip(v))
        haskey(d, key) && d[key] != val &&
            throw(ArgumentError("configuration key $key is given twice, in different case, " *
                                "with different values: \"$(d[key])\" and \"$val\""))
        d[key] = val
    end
    d
end

function CSConfig(dict::Dict{String,String})
    # Known keys are lowercase and stripped, so only a dict with other keys
    # (never one built from init_config with update!) needs normalising.
    if !all(k -> haskey(INI_DEFAULTS, k), keys(dict))
        dict = normalize_config(dict)
        unknown = sort!([k for k in keys(dict) if !haskey(INI_DEFAULTS, k)])
        isempty(unknown) ||
            @warn("Ignoring unrecognised configuration keys: $(join(unknown, ", "))")
    end
    _parse_fields(dict)
end

# One positional constructor call with every field's parser known statically,
# as in `_copy_with` below. A key absent from `dict` takes its init_config
# value; values are stripped of surrounding whitespace.
@generated function _parse_fields(dict::Dict{String,String})
    args = [:(_from_ini($(fieldtype(CSConfig, f)), $(String(f)),
                        strip(get(dict, $(String(f)), INI_DEFAULTS[$(String(f))]))))
            for f in fieldnames(CSConfig)]
    :(CSConfig($(args...)))
end

# Inputs read only when their flag is set.
const OPTIONAL_FILES = (:use_polygons => :polygon_file, :use_mask => :mask_file,
                        :use_included_pairs => :included_pairs_file,
                        :use_variable_source_strengths => :variable_source_file,
                        :use_reclass_table => :reclass_file)

# Inputs that may be rasters, and so may need GDAL.
const RASTER_FILES = (:habitat_file, :point_file, :polygon_file, :mask_file,
                      :source_file, :ground_file)

# Options parsed for INI compatibility with Circuitscape 4 but not implemented.
const UNIMPLEMENTED = (:low_memory_mode, :preemptive_memory_release)

"""
    validate(cfg::CSConfig)

Fail fast, before any data is read, on configurations that cannot run:
input files that do not exist, an output directory that does not exist, and
options that Circuitscape.jl parses for compatibility with the Python version
but has never implemented. Each message names the INI key so the user knows
what to change.
"""
function validate(cfg::CSConfig)
    problems = String[]

    function need(key)
        path = getfield(cfg, key)
        if _is_unset(path)
            push!(problems, "$key is not set")
        elseif !isfile(path)
            push!(problems, "$key = \"$path\" does not exist")
        end
    end

    need(:habitat_file)
    if is_advanced(cfg)
        need(:source_file)
        need(:ground_file)
    else
        need(:point_file)
    end
    for (flag, key) in OPTIONAL_FILES
        getfield(cfg, flag) && need(key)
    end

    if _is_unset(cfg.output_file)
        push!(problems, "output_file is not set")
    elseif basename(cfg.output_file) == ""
        # A trailing slash names a directory, not a file: the outputs would
        # be written as `dir/_cum_curmap.asc`, `dir/.ini` and so on.
        push!(problems, "output_file \"$(cfg.output_file)\" names a directory; give a file name in it, such as \"$(joinpath(cfg.output_file, "results.out"))\"")
    else
        outdir = dirname(cfg.output_file)
        outdir == "" || isdir(outdir) ||
            push!(problems, "output directory \"$outdir\" (from output_file) does not exist")
    end

    # GDAL is optional: only GeoTIFF (or other non-ASCII) rasters need it.
    if !gdal_loaded()
        for key in RASTER_FILES
            path = getfield(cfg, key)
            _is_unset(path) || !needs_gdal(path) ||
                push!(problems, gdal_missing_message("$key = \"$path\""))
        end
        cfg.write_as_tif &&
            push!(problems, gdal_missing_message("write_as_tif = True"))
    end

    for key in UNIMPLEMENTED
        getfield(cfg, key) &&
            push!(problems, "$key is not implemented in Circuitscape.jl; remove it or set it to False")
    end

    isempty(problems) ||
        throw(ArgumentError("Invalid configuration:\n  - " * join(problems, "\n  - ")))
    cfg
end

# Predicates on a CSConfig. These replace the RasterFlags / NetworkFlags /
# OutputFlags structs that used to be rebuilt from cfg at every entry point:
# with a single representation there is no second copy for an option to be
# parsed into and then never read from.
is_raster(cfg::CSConfig) = cfg.data_type == dt_raster
is_network(cfg::CSConfig) = cfg.data_type == dt_network
is_pairwise(cfg::CSConfig) = cfg.scenario == sc_pairwise
is_advanced(cfg::CSConfig) = cfg.scenario == sc_advanced
is_onetoall(cfg::CSConfig) = cfg.scenario == sc_onetoall
is_alltoone(cfg::CSConfig) = cfg.scenario == sc_alltoone
policy(cfg::CSConfig) = Symbol(_to_ini(cfg.remove_src_or_gnd))

"""
    CSConfig(cfg::CSConfig; field = value, ...)

Copy of `cfg` with the given fields replaced.
"""
CSConfig(cfg::CSConfig; kwargs...) = _copy_with(cfg, values(kwargs))

# The field names to override are part of the NamedTuple's type, so the copy
# is emitted as one positional call with every argument's type known to
# inference: 0 bytes, instead of building a 46-entry NamedTuple at run time
# and dispatching the keyword constructor on it (~10 us and 17 KB per copy,
# once per moving window in `compute_omniscape_current`). The positional
# constructor converts each value to its field type, as the keyword one does.
@generated function _copy_with(cfg::CSConfig, kw::NamedTuple{names}) where {names}
    for n in names
        n in fieldnames(CSConfig) ||
            return :(throw(ArgumentError($("CSConfig has no field $n"))))
    end
    args = [n in names ? :(kw.$n) : :(getfield(cfg, $(QuoteNode(n))))
            for n in fieldnames(CSConfig)]
    :(CSConfig($(args...)))
end

"""
    parse_config(path::AbstractString) -> CSConfig

Read a Circuitscape INI file, with IniFile.jl, as Circuitscape 4 (Python's
ConfigParser) did: each line is `key = value` or `key: value`, keys and option
values are case-insensitive (file paths keep their case), blank lines and
comment lines (starting with `#` or `;`) are skipped, and a `;` after
whitespace starts a comment. Any other line is an error naming the file and
line. Sections do not matter: all keys are read, from every section. When a key
is given more than once the later one wins; a section header that appears again
continues its first section, as in Circuitscape 4. Every key is parsed on top
of the defaults from [`init_config`](@ref), so keys that are absent take their
default value. A value that an enumerated or boolean option does not accept
raises an `ArgumentError` naming the key and the accepted values; unknown keys
produce a warning and are ignored. File existence is not checked here but in
[`validate`](@ref).
"""
function parse_config(path::AbstractString)
    ini = read(Inifile, path; lowercase_keys = true, inline_comment_prefixes = (";",),
               strict = true)
    CSConfig(merge!(init_config(), merged(ini)))
end

"""
    init_config() -> Dict{String,String}

The default configuration as INI strings, one entry per key `compute` accepts.
Modify the entries you need and pass the dictionary to [`compute`](@ref), or
convert it with `CSConfig(dict)`. File paths default to INI Builder
placeholders such as `"(Browse for a resistance file)"`, which
[`validate`](@ref) treats as unset. The Circuitscape 4 keys that
Circuitscape.jl ignores are present too, so that old INI files load without
warnings; the options reference lists them.
"""
init_config() = copy(INI_DEFAULTS)

"""
    update!(cfg::Dict{String,String}, new) -> cfg

Merge the INI keys and values in `new` (any iterable of `key => value`
pairs, typically a `Dict{String,String}`) into the configuration dictionary
`cfg`, overwriting existing entries; keys are lowercased and stripped, as in
an INI file. Mutates and returns `cfg`. This is how a
caller layers its own settings on top of [`init_config`](@ref) before passing
the result to [`compute`](@ref) or [`compute_omniscape_current`](@ref); no
values are parsed or checked here (that happens in `CSConfig(cfg)`).
"""
function update!(cfg, new)
    for (key,val) in new
        cfg[_ini_key(key)] = val
    end
    cfg
end

"""
    output_prefix(cfg::CSConfig) -> String
    output_prefix(output_file::AbstractString) -> String

The prefix every output file name starts with: `output_file` without its
extension (`splitext`), so `output/run.out` gives `output/run` and
`out/run` stays `out/run`. Splitting on the substring `".out"` (as was done
before) truncated any path whose *directory* contained `.out`.
"""
output_prefix(cfg::CSConfig) = output_prefix(cfg.output_file)
output_prefix(output_file::AbstractString) = first(splitext(String(output_file)))

"""
    write_config(cfg::CSConfig, path = cfg.output_file)
    write_config(cfg::Dict{String,String}, path = ...)

Write the configuration, in Circuitscape INI form, to `path`; by default the
file named by `cfg.output_file`, as in every earlier release. The file holds
every option, so `parse_config(path)` gives back `cfg` (the `[Version]`
section records the version of Circuitscape.jl that wrote it, not
`cfg.version`).

`compute` records the options of each run in `<prefix>.ini` (see
[`output_prefix`](@ref)), as Circuitscape 4 did, except when that file is
the INI being run, which is never overwritten.
"""
function write_config(cfg::CSConfig, path::AbstractString = cfg.output_file)
    ini = Inifile()
    for (section, keys) in INI_SECTIONS, key in keys
        value = key === :version ? string(pkgversion(@__MODULE__)) :
                                   _to_ini(getfield(cfg, key))
        IniFile.set(ini, section, String(key), value)
    end
    write(path, ini)
    nothing
end

# Keep backward compat for Dict-based write_config
write_config(cfg::Dict{String,String}) = write_config(CSConfig(cfg))
write_config(cfg::Dict{String,String}, path::AbstractString) = write_config(CSConfig(cfg), path)
