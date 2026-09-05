# Options and Flags

## INI Argument Reference

All Circuitscape configuration is done through an `.ini` file. Below is a complete reference of all arguments, organized by section, with their types, default values, and descriptions.

### Circuitscape Mode

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `data_type` | String | `raster` | Input data type. Values: `raster`, `network`. |
| `scenario` | String | `not entered` | Modeling mode. Values: `pairwise`, `advanced`, `one-to-all`, `all-to-one`. |

### Habitat Raster or Graph

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `habitat_file` | Path | — | Path to the resistance/conductance map file (raster or network). ESRI ASCII grids (`.asc`, optionally gzipped) are read natively; GeoTIFF and other GDAL formats require `using ArchGDAL` (see *GeoTIFF Support*). |
| `habitat_map_is_resistances` | Boolean | `True` | If `True`, habitat map values are resistances. If `False`, values are conductances. |

### Connection Scheme for Raster Habitat Data

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `connect_four_neighbors_only` | Boolean | `False` | If `True`, connect each cell to 4 cardinal neighbors only. If `False`, connect to all 8 neighbors. |
| `connect_using_avg_resistances` | Boolean | `False` | If `True`, use average resistance for cell connections. If `False`, use average conductance. |

### Short-Circuit Regions (Polygons)

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `use_polygons` | Boolean | `False` | If `True`, read a short-circuit region file. Cells within each region are collapsed into a single node with zero resistance. |
| `polygon_file` | Path | — | Path to the short-circuit region map. Must have the same cell size and extent as the resistance grid. |

### Options for Advanced Mode

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `source_file` | Path | — | Path to the current source file (raster or text list). Values specify source strength in amps. |
| `ground_file` | Path | — | Path to the ground point file (raster or text list). Values specify resistance or conductance to ground. |
| `ground_file_is_resistances` | Boolean | `True` | If `True`, ground file values are resistances. If `False`, conductances. Set to `False` with value 0 to connect directly to ground. |
| `use_unit_currents` | Boolean | `False` | If `True`, all current sources are set to 1 Amp regardless of values in the source file. |
| `use_direct_grounds` | Boolean | `False` | If `True`, all ground nodes are tied directly to ground (R=0) regardless of values in the ground file. |
| `remove_src_or_gnd` | String | `keepall` | When a source and ground are at the same node: `keepall`, `rmvsrc`, `rmvgnd`, or `rmvall`. |

### Options for Pairwise, One-to-All, and All-to-One Modes

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `point_file` | Path | — | Path to focal node file (raster or text list). Each focal node must have a unique positive integer ID. |
| `use_included_pairs` | Boolean | `False` | If `True`, only run calculations on a subset of focal node pairs specified in the included pairs file. |
| `included_pairs_file` | Path | — | Path to file specifying pairs to include or exclude from calculations. |

### Options for One-to-All and All-to-One Modes

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `use_variable_source_strengths` | Boolean | `False` | If `True`, read per-node source strengths from a file instead of using 1 Amp for all. |
| `variable_source_file` | Path | `None` | Path to text file with focal node IDs and corresponding source strengths. |

### Mask File

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `use_mask` | Boolean | `False` | If `True`, apply a mask to the resistance map. Cells with negative, zero, or NODATA values in the mask are dropped. |
| `mask_file` | Path | `None` | Path to the raster mask file. |

### Output Options

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `output_file` | Path | — | Base path and filename for all output files. |
| `write_cur_maps` | Boolean | `False` | If `True`, write current maps for each iteration and a cumulative map. |
| `write_volt_maps` | Boolean | `False` | If `True`, write voltage maps for each iteration. |
| `write_cum_cur_map_only` | Boolean | `False` | If `True`, calculate current maps for each iteration but only write the cumulative (summed) map to disk. |
| `write_max_cur_maps` | Boolean | `False` | If `True`, write a map showing the maximum current value at each cell across all iterations. |
| `set_null_currents_to_nodata` | Boolean | `False` | If `True`, cells that are NODATA in the resistance map (and so have no graph node) are written as NODATA in current maps instead of 0. |
| `set_null_voltages_to_nodata` | Boolean | `False` | If `True`, cells that are NODATA in the resistance map are written as NODATA in voltage maps instead of 0. |
| `set_focal_node_currents_to_zero` | Boolean | `False` | If `True`, the cells of the *active* focal nodes carry zero current in each raster current map: the two focal nodes of the pair in pairwise mode, every focal node in one-to-all and all-to-one. Cumulative and maximum maps then show only current that flows *through* a focal region while other pairs are active. Raster current maps only; network node currents and resistances are unaffected. |
| `compress_grids` | Boolean | `False` | Accepted for compatibility with Circuitscape 4. Output grids are currently written uncompressed; gzipped *input* grids (`.asc.gz`) are read. |
| `log_transform_maps` | Boolean | `False` | If `True`, log10-transform values in output current maps. Cells with zero current are set to NODATA. |
| `write_as_tif` | Boolean | `False` | If `True`, write current and voltage maps as GeoTIFF (LZW-compressed) instead of ESRI ASCII grids. Requires the optional GDAL support: `Pkg.add("ArchGDAL")`, then `using ArchGDAL` before running. |

### Calculation Options

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `solver` | String | `cg+amg` | Linear solver to use. Values: `cg+amg` (iterative, recommended for large grids), `cholmod` (direct, uses more memory), `accelerate` (direct, macOS only, requires [AppleAccelerate.jl](https://github.com/JuliaLinearAlgebra/AppleAccelerate.jl)), `pardiso` (direct, requires [Pardiso.jl](https://github.com/JuliaSparse/Pardiso.jl)). Also accepted: `amg+cg`, `cholesky`, `cholfact`, `mklpardiso`, `apple_accelerate`. Any other value is an error. |
| `precision` | String | `double` | Floating-point precision. Values: `double`, `single` (either case). Single precision is supported and tested for `cg+amg`, `cholmod` and `accelerate`; `pardiso` is double-only and switches to double with a warning. See the accuracy note under *Convergence Checks* in [On Solvers and Computation Time](compute.md). |
| `use_64bit_indexing` | Boolean | `True` | If `True`, use 64-bit integer indexing. Required for very large grids. |
| `cholmod_batch_size` | Integer | `1000` | Number of right-hand sides the direct solvers (`cholmod`, `pardiso`, `accelerate`) solve at once in pairwise mode. Postprocessing of each batch runs across threads. Ignored by `cg+amg`. |
| `residual_tolerance` | Float | `auto` | Every linear solve is accepted only if its true relative residual `‖Gv − b‖ / ‖b‖` is below this. `auto` means `1e-4` in double precision and `1e-3` in single. With `cg+amg`, a solve that stops short is retried with tighter Krylov tolerances before the run is aborted. |
| `parallelize` | Boolean | `False` | If `True`, run pair solves (or per-focal-node solves) in parallel using Julia threads. Start Julia with `julia -t N` for N threads; see *Parallelism* below. |
| `low_memory_mode` | Boolean | `False` | Circuitscape 4 option, not implemented in Circuitscape.jl. Setting it to `True` is an error so that it cannot be relied on silently. |
| `preemptive_memory_release` | Boolean | `False` | Circuitscape 4 option, not implemented in Circuitscape.jl. Setting it to `True` is an error. |

### Reclassification

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `use_reclass_table` | Boolean | `False` | If `True`, replace values in the habitat raster using the lookup table in `reclass_file` before the run. Raster mode only. |
| `reclass_file` | Path | — | Text file with two whitespace-separated columns, `old new`, one pair per line. Every cell whose value equals `old` is set to `new`; other cells, including NODATA, are left alone. Values are replaced in a single pass, so `1 2` followed by `2 3` sends 1 to 2, not to 3. Lets one large raster be reused with many resistance assignments without exporting a new file each time. |

### Logging

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `log_file` | Path | `None` | Path to log file. If `None`, no file logging. |
| `log_level` | String | `INFO` | Logging level. Values: `DEBUG`, `INFO`, `WARNING`, `ERROR`, `CRITICAL` (case-insensitive; `WARN` is also accepted; `CRITICAL` is treated as `ERROR`). Any other value is an error. `DEBUG` logs every pair solve and prints a timing table at the end of the run. |
| `suppress_messages` | Boolean | `False` | If `True`, suppress informational messages on the console. Warnings and errors are still shown, and a `log_file` still receives everything. |

The Circuitscape 4 keys `profiler_log_file`, `print_timings`, `print_rusages` and `screenprint_log` are accepted without a warning and have no effect.

---

## Configuration Validation

Since 6.0 a configuration is checked in full before any data is read, and a problem is reported as an error rather than silently turning into a different run.

**Values.** `solver`, `scenario`, `data_type`, `precision`, `log_level`, `remove_src_or_gnd` and every boolean option accept only the values listed in the tables above. Anything else is an error naming the key and the valid values. In earlier versions `solver = cholmodd` ran `cg+amg`, `scenario = pairwsie` ran pairwise and `write_cur_maps = Ture` wrote nothing; each of these is now refused:

```
ArgumentError: solver = "cholmodd" is not recognised; expected one of: cg+amg, cholmod, pardiso, accelerate
ArgumentError: write_cur_maps = "Ture" is not recognised; expected one of: True, False
```

A missing `scenario` is an error (`scenario is not set; expected one of: pairwise, advanced, one-to-all, all-to-one`) rather than an implicit pairwise run.

**Files.** Every input file the run will actually use must exist: `habitat_file` always; `point_file` in pairwise, one-to-all and all-to-one; `source_file` and `ground_file` in advanced mode; `polygon_file`, `mask_file`, `included_pairs_file`, `variable_source_file` and `reclass_file` only when the matching `use_*` option is `True`. The directory of `output_file` must exist. An INI Builder placeholder such as `(Browse for a resistance file)` counts as unset. All problems are collected and reported in one message, so a configuration is fixed in one round trip:

```
ArgumentError: Invalid configuration:
  - habitat_file = "resistance.asc" does not exist
  - point_file = "nodes.asc" does not exist
  - output directory "output" (from output_file) does not exist
  - low_memory_mode is not implemented in Circuitscape.jl; remove it or set it to False
```

**Never-implemented options.** `low_memory_mode` and `preemptive_memory_release` are Circuitscape 4 options that Circuitscape.jl parses but has never implemented. Setting either to `True` is an error, so that a large run cannot rely on them by mistake.

**GeoTIFF without GDAL.** If ArchGDAL is not loaded, a configuration that names a GeoTIFF input (detected by file content, not extension) or sets `write_as_tif = True` is refused up front with the installation hint (`Pkg.add("ArchGDAL")`, then `using ArchGDAL`). See *GeoTIFF Support* in [On Solvers and Computation Time](compute.md).

**Unknown keys** produce a warning (`Ignoring unrecognised configuration keys: ...`) and are otherwise ignored, so a misspelt key is visible without breaking old INI files.

**Circuitscape 4 spellings** still work: `remove_src_or_gnd = not entered` means `keepall`, Python log level names (`WARNING`, `CRITICAL`, any case) are read, and `residual_tolerance = None` means automatic.

## Data Type

Set `data_type` to `raster` or `network` to choose whether you will be analyzing raster grid or network data.

## Modeling Mode

Circuitscape is run in one of four modes (set via `scenario`). Pairwise and advanced modes are available for both raster and network data types. The one-to-all and all-to-one modes are available for raster data only.

## Resistance Map or Network/Graph

The resistance file (`habitat_file`) specifies the ability of each cell in a landscape or link in a network to carry current. File formats are described in the _Input file formats_ section below.

Most users code their data in terms of resistances (with higher values denoting greater resistance to movement). Set `habitat_map_is_resistances = False` to specify conductances instead (conductance is the reciprocal of resistance; higher values indicate greater ease of movement).

Zero and infinite values for conductances and resistances represent special cases. Infinite resistances are coded as NODATA values in input resistance grids, or as zero or NODATA in input conductance grids; these are treated as complete barriers and are disconnected from all other cells. For raster analyses, cells with zero resistance (infinite conductance) can be specified using a separate short-circuit region file as described below.

## Focal Nodes (Pairwise, One-to-All, and All-to-One Modes)

The focal node file (`point_file`) specifies locations of nodes between which effective resistance and current flow are to be calculated. **Each focal node should have a unique positive integer ID.** Files may be text lists specifying coordinates or raster grids. When a grid is used, it must have the same cell size and extent as the resistance grid. Cells that do not contain focal nodes should be coded with NODATA values.

For raster analyses, focal nodes may occur at points (single cells on the resistance grid) or across regions. For the latter, a single ID would occupy more than one cell in a grid or more than one pair of coordinates in a text list (and falling within more than one cell in the underlying resistance grid). Cells within a single region are collapsed into a single node. The difference from short-circuit regions is that a focal region will be "burned in" to the resistance grid only for pairwise calculations that include that focal node. Focal regions need not be contiguous. For large grids, focal regions may require more computation time. When calculating resistances on large raster grids and not creating voltage or current maps, focal points will run much more quickly.

### Parallelism

Circuitscape uses Julia's native threads. Two things are needed: start Julia with threads (`julia -t N`, or set the `JULIA_NUM_THREADS` environment variable before starting Julia) **and** set `parallelize = True` in the INI file. With `parallelize = False`, the default, the run is serial however many threads Julia has. The log line `Starting up Circuitscape to use N threads in parallel` confirms both are in effect.

What runs in parallel depends on the mode and solver:

- **Pairwise with `cg+amg`** (raster and network): the pair solves. All pairs of a connected component are collected and split into equal-sized chunks, about four per thread, each solved on its own task with its own preconditioner workspace. Every pair costs about the same, so equal counts are equal work. This includes the default configuration with no maps requested, where the resistance shortcut solves one system per focal node rather than one per pair: those solves used to run on a single task and are now spread across the threads too.
- **Pairwise with a direct solver** (`cholmod`, `pardiso`, `accelerate`): the factorization is computed once per connected component and `cholmod_batch_size` right-hand sides are solved at a time inside the solver library; the postprocessing of each batch (node currents, map accumulation and writing) runs across threads.
- **One-to-all and all-to-one**: one task per focal node. Each is an independent solve; the per-node current maps are accumulated on the main thread afterwards.
- **Advanced mode** is a single solve per connected component and is not parallelized.

Julia 1.12 starts an interactive thread in addition to the default worker threads. Circuitscape schedules its tasks on the default thread pool only, so `Threads.nthreads()`, the number given by `-t`, is what matters and the interactive thread changes nothing. BLAS is set to a single thread at startup because the workload is sparse.

Set `log_level = DEBUG` to get a timing table at the end of the run with the `solve linear system` and `postprocess` time per task.

## Advanced Mode Options

### Current Source File

The source file (`source_file`) specifies locations and strengths, in amps, of current sources. Either a raster or a text list may be used. Rasters must have the same cell size, projection, and extent as the resistance grid, and cells that do not contain current sources should be coded with NODATA values. Current sources may be positive or negative (i.e., they may inject current into the grid or pull current out). Similarly, grounds may either serve as a sink for current or may contribute current if there are negative current sources in the grid.

### Ground Point File

The ground file (`ground_file`) specifies locations of ground nodes and resistances or conductances of resistors tying them to ground. Either a raster or a text list may be used. Rasters must have the same cell size, projection, and extent as the resistance grid, and cells that do not contain grounds should be coded with NODATA values. If a direct (R = 0) ground connection conflicts with a current source, the ground will be removed unless `remove_src_or_gnd` is set to `rmvgnd` or `keepall`.

Set `ground_file_is_resistances = False` if your ground point file specifies connections to ground in terms of conductance instead. To tie cells directly to ground, keep `ground_file_is_resistances = True` and set values in the ground point file to zero.

### Unit Currents and Direct Grounds

Set `use_unit_currents = True` to force all current sources to 1 Amp, regardless of the value specified in the source file. Set `use_direct_grounds = True` to tie all ground nodes directly to ground (R=0) regardless of the ground file values.

### Source/Ground Conflicts

When a cell is connected both to a current source and to ground, `remove_src_or_gnd` determines the behavior: `keepall` (default), `rmvsrc`, `rmvgnd`, or `rmvall`. When using `keepall`, if a source is tied directly to ground (zero resistance), the ground connection will be removed.

## Output Options

### Current Maps

When `write_cur_maps = True`, current maps will be generated for every pair of focal nodes in pairwise mode, or for the source/ground configuration in advanced mode. Current maps have the same dimension as the input files, with values at each cell representing the amount of current flowing through the node. In pairwise mode, a current map file will be created for each focal node pair, and a cumulative (additive) map will also be written. (Note that for a given pair of focal nodes, current maps are identical regardless of which node is the source and which is the ground due to symmetry.) For advanced mode, a single map will be written. Such maps can be used to identify areas which contribute most to connectivity between focal points (McRae et al. 2008).

![](assets/image24.png)

**Fig. 1.** Current map used to predict important connectivity areas between core habitat patches (green polygons, entered as focal regions) for mountain lions. Warmer colors indicate areas with higher current density. "Pinch points," or areas where connectivity is most tenuous, are shown in yellow. Research Collaborators: Brett Dickson and Rick Hopkins, Live Oak Associates.

### Voltage Maps

When `write_volt_maps = True`, voltage maps are written. In pairwise mode, these give node voltages observed for each focal node pair if one node were connected to a 1 amp current source and the other to ground. In advanced mode, voltage maps show voltages resulting from the source and ground configurations.

## Calculation Options

### Cell Connectivity

For raster operations, Circuitscape creates a graph by connecting cells to their neighbors. Set `connect_four_neighbors_only = True` for 4 cardinal neighbors only (default is 8, including diagonals).

### Average Resistance vs Conductance

Set `connect_using_avg_resistances = True` to connect cells by their average resistance instead of average conductance (the default).

The distinction is particularly important when connecting cells with zero or infinite values. When average resistances are used, first-order neighbors are connected by resistors with resistance: _Rab_ = (_Ra_ + _Rb_) / 2, and second-order (diagonal) neighbors by: _Rab_ = sqrt(2) * (_Ra_ + _Rb_) / 2. When average conductances are used, first-order neighbors are connected by: _Gab_ = (_Ga_ + _Gb_) / 2, and second-order (diagonal) neighbors by: _Gab_ = (_Ga_ + _Gb_) / (2 * sqrt(2)).

## Mapping Options

### Maximum Current Maps

In pairwise, one-to-all, and all-to-one modes, current maps are created for every iteration. By default, Circuitscape writes a cumulative map showing the sum of values at each node or grid cell across all iterations. Set `write_max_cur_maps = True` to also write a map showing the maximum current value at each cell across iterations.

### Cumulative Maps Only

Set `write_cum_cur_map_only = True` to calculate current maps for each iteration but only write the cumulative (and optionally maximum) map to disk. This saves disk space when many iterations are performed.

### Compress Output Grids

Set `compress_grids = True` to compress output ASCII grids using gzip. This can be useful when many large maps will be written.

### Log-Transform Current Maps

Set `log_transform_maps = True` to apply a log10 transform to current densities in output maps. Cells with zero current will be set to NODATA values.

### Set Focal Node Currents to Zero

When `set_focal_node_currents_to_zero = True`, focal nodes have zero current in raster current maps when they are activated. For pairwise mode, cumulative maps will still show currents flowing through focal regions from other pairs being activated. This helps show the importance of each focal region for connecting others (see Dickson et al. 2013). In one-to-all and all-to-one modes every focal node is active on each solve, so all focal cells are zeroed. Applies to raster current maps; network node currents are unaffected.

## Optional Input Files

### Mask File

Set `use_mask = True` and provide `mask_file` to apply a raster mask. Cells with negative, zero, or NODATA values in the mask will be dropped from the resistance map (treated as complete barriers). Positive integer cells will be retained. File should only contain integers and be in raster format.

### Short-Circuit Region Map

Set `use_polygons = True` and provide `polygon_file` to load short-circuit regions. These act as areas of zero resistance, providing patches through which current gets a "free ride." Each region should have a unique positive integer identifier; cells within each region are merged into a single node, including non-adjacent cells (regions need not be contiguous). Non-region areas should be stored as NODATA values. The file must have the same cell size and extent as the resistance grid.

### Variable Source Strengths

In one-to-all and all-to-one modes, set `use_variable_source_strengths = True` and provide `variable_source_file` with focal node IDs and corresponding source strengths. The file should be a text list with two columns (ID followed by source strength). All nodes not in the list will default to 1 Amp.

### Include/Exclude Focal Node Pairs

Set `use_included_pairs = True` and provide `included_pairs_file` to restrict calculations to a subset of focal node pairs. Users can either identify pairs to include or pairs to exclude, as specified in the first line of the file. This affects all modes except advanced mode. Files should be in tab-delimited text with a .txt extension.

## Input Raster Format

Raster input maps should be stored in Arc/Info ASCII grid or GeoTIFF format, as exported by standard GIS packages. ASCII grids (optionally gzipped, with `xllcorner`/`yllcorner` or `xllcenter`/`yllcenter`, `cellsize` or `dx`/`dy`) are read natively; GeoTIFF needs the optional ArchGDAL package (see *GeoTIFF Support* in [On Solvers and Computation Time](compute.md)). Set `write_as_tif = True` to produce GeoTIFF output instead of ASCII grids. For focal nodes, the value stored in each grid location refers to the focal node ID, and a single ID can occupy more than one cell (IDs must be positive integers). For current sources, the grid value specifies the source strength in amps. For grounds, the grid value specifies either the resistance or conductance of the resistor tying each ground node to ground.

The ASCII raster format is as follows:

**Header:**

    ncols        <Number of columns>
    nrows        <Number of rows>
    xllcorner    <X coordinate of lower left corner>
    yllcorner    <Y coordinate of lower left corner>
    cellsize     <size of each cell>
    NODATA_value <Code for cells with no habitat, focal nodes, sources or grounds>

**Body (grid data):**

Numeric data only. Columns are delimited with tabs and rows are delimited with new line characters.

**Examples**

Below is a 10 x 10 resistance map. Cells with infinite resistance are assigned NODATA values (-9999):

    ncols         10
    nrows         10
    xllcorner     1
    yllcorner     1
    cellsize      1
    NODATA_value  -9999
    130    168    153    -9999  14     12    13     107    140    171
    104    3      2      -9999  13     158   12     14     13     114
    124    2      2      12     -9999  -9999 13     161    4      5
    184    5      4      14     13     14    -9999  13     4      4
    105    143    103    169    -9999  115   10     -9999  166    14
    187    1      163    188    121    142   14     175    -9999  10
    198    11     110    115    149    2     2      164    3      -9999
    100    11     193    14     12     4     2      1      11     13
    -9999  11     12     11     10     12    167    157    181    157
    -9999  -9999  122    134    12     157   192    184    190    172

Below is a 10 x 10 focal region map. Groups of cells have been coded as focal regions that will be treated as "core area polygons" to be connected in circuit analyses. All cells within each focal region will be collapsed into a single node (even the non-contiguous cell in region #1) when that region is activated in pairwise, one-to-all, or all-to-one analyses. This format is identical to the short-circuit region file format.

    ncols                10
    nrows                10
    xllcorner            1
    yllcorner            1
    cellsize             1
    NODATA_value -9999
    -9999  -9999  -9999  -9999  -9999  -9999  -9999  -9999  -9999  -9999
    -9999  1      1      -9999  -9999  -9999  -9999  -9999  -9999  -9999
    -9999  1      1      -9999  -9999  -9999  -9999  -9999  3      3
    -9999  1      1      -9999  -9999  -9999  -9999  -9999  3      3
    -9999  -9999  -9999  -9999  -9999  -9999  -9999  -9999  -9999  -9999
    -9999  1      -9999  -9999  -9999  -9999  -9999  -9999  -9999  -9999
    -9999  -9999  -9999  -9999   2      2     -9999  -9999  -9999  -9999
    -9999  -9999  -9999  -9999   2      2      2     -9999  -9999  -9999
    -9999  -9999  -9999  -9999  -9999  -9999  -9999  -9999  -9999  -9999
    -9999  -9999  -9999  -9999  -9999  -9999  -9999  -9999  -9999   -9999

Note that regions 1 and 2 are well-connected by a low-resistance corridor in the resistance map above. Region 3 is connected to the other two regions only if cells are connected to their eight neighbors. In the four-neighbor case, region 3 would be completely isolated.

## Text List File Format

For network/graph operations, resistor networks, focal nodes, current sources, and grounds should be stored as text lists (saved with a ".txt" extension). To specify a network of resistors, three columns are used. The first and second columns give the node IDs being connected by a resistor, and the third column gives the resistance value. For example, the simple circuit:

![](assets/SimpleNetworkWithNumbers2.png)

can be defined by the following text list:

        0    1    1
        1    2    1
        1    3    1
        2    4    1
        3    4    1

**Please note:** typically, there should just be one entry for each pair of connected nodes. If there are two entries for a single pair in the form of (node1, node2, value1) and (node2, node1, value2), these will be considered parallel resistors and their conductances will be summed. For example, if the above text list had an extra entry for node pair (4, 3) like this:

        0    1    1
        1    2    1
        1    3    1
        2    4    1
        3    4    1
        4    3    1

then the resistance between nodes 3 and 4 in the resulting graph would be 1/2 ohm.

For advanced mode, current sources and grounds are also stored as text lists. The above circuit can be expanded to include a current source and grounds with two extra input files. For example, we can add a 1 Amp current source at node 0 with a file that looks like this:

        0    1

To tie node 4 directly to ground (i.e. to connect it to ground with a wire that has a resistance of 0 Ohms) and connect the remaining nodes to ground with resistors, we can use a file that looks like this:

        0    99
        1    33
        2    49.5
        3    49.5
        4    0

The resulting circuit would look like this (from McRae et al. 2008):

![](assets/AdvancedNetwork.png)

For **raster** operations, you can also store focal nodes, current sources, and grounds as text lists (saved with a ".txt" extension). For each node referenced in a text list, a value and X and Y coordinates are specified as shown below.

        Value1 X1 Y1
        Value2 X2 Y2
        …

Note: X and Y are geographical coordinates, not row and column numbers.

Example text list (a partial list of the cell locations in the focal region map above; coordinates are for cell centroids):

        1    2.5    9.5
        1    3.5    9.5
        1    2.5    8.5
        1    3.5    8.5
        1    2.5    7.5
        1    3.5    7.5
        1    2.5    5.5
        2    6.5    4.5
        ...

For focal nodes, the value field references the focal node ID; values must be positive integers, and a single ID can occupy more than one pair of coordinates (and more than one cell in the underlying resistance grid). For current sources, the value field references the source strength in amps. For grounds, the value field references either the resistance or conductance of the resistor tying each ground node to ground.

## Include/Exclude File Format

This file is used when `use_included_pairs = True`, and affects all modes except advanced mode. There are two file formats that can be used. The first is the simplest, and gives a list of pairs to include in calculations, or pairs to exclude, as specified in the first line of the file. For example, if there are five focal nodes, numbered 1-5, and the following list is entered, only pairs (1,2), (1,3), and (1,5) will be analyzed:

    mode    include
    1        2
    1        3
    1        5

Similarly, if the first line in the above file read:

    mode     exclude

all pairs except (1,2), (1,3), and (1,5) would be analyzed.

The second method uses a matrix identifying which pairs of focal nodes to connect. The file specifies minimum and maximum values in the matrix to consider a pair connected. This method can be useful when used with a distance matrix to only run analyses between points separated by a minimum distance, or by a distance equal to or less than a maximum distance. Note: any focal node not in the matrix will be dropped from analyses. Entries on the diagonal are ignored. For example, in the following matrix, only pairs with entries between 2 and 50 are connected. Pairs (1,2), (2,4), and (3,4) will not be analyzed. Focal node 5 will be dropped entirely:

    min    2
    max    50
    0     1     2     3     4     5
    1     0     100   6.67  7     1
    2     100   0     11    1     60
    3     6.67  11    0     -1    100
    4     7     1     -1    0     0
    5     1     60     100  0     0

Make sure to include a zero in the upper-left corner of the matrix.

Files should be in tab-delimited text with a .txt extension.

## Output Files

Every output name starts with `output_file` with its `.out` extension removed; below this prefix is written `<out>`. The configuration as run is written, in INI form, to `output_file` itself, so the options are recorded next to the results. Focal node IDs in file names are the IDs from the focal node file.

### Which files are written

Raster modes write ESRI ASCII grids (`.asc`), or GeoTIFF (`.tif`) when `write_as_tif = True`. When the input raster carries a projection it is written alongside each `.asc` as a `.prj` sidecar.

| File | Mode | Written when |
|------|------|--------------|
| `<out>_resistances.out`, `<out>_resistances_3columns.out` | pairwise | always |
| `<out>_curmap_<i>_<j>.asc` | pairwise | `write_cur_maps = True` and `write_cum_cur_map_only = False` |
| `<out>_cum_curmap.asc` | pairwise | `write_cur_maps = True` or `write_cum_cur_map_only = True` |
| `<out>_max_curmap.asc` | pairwise | `write_max_cur_maps = True`, together with `write_cur_maps` or `write_cum_cur_map_only` |
| `<out>_voltmap_<i>_<j>.asc` | pairwise | `write_volt_maps = True` |
| `<out>_curmap_<id>.asc` (one per focal node) | one-to-all, all-to-one | `write_cur_maps = True` or `write_cum_cur_map_only = True` (`write_cum_cur_map_only` does not suppress the per-node maps in these modes) |
| `<out>_cum_curmap.asc`, `<out>_max_curmap.asc` | one-to-all, all-to-one | as in pairwise |
| `<out>_voltmap_<id>.asc` | one-to-all, all-to-one | `write_volt_maps = True` |
| `<out>_curmap.asc` | advanced | `write_cur_maps = True` or `write_cum_cur_map_only = True` |
| `<out>_voltmap.asc` | advanced | `write_volt_maps = True` |

One-to-all and all-to-one do not write a resistances file; their result (one value per focal node) is the return value of `compute`. `set_focal_node_currents_to_zero`, `log_transform_maps` and `set_null_currents_to_nodata` are applied to each per-pair or per-node map before it is accumulated, so they affect the cumulative and maximum maps as well.

Network modes write tab-separated text. Branch current files omit branches carrying less than `1e-6` A.

| File | Mode | Written when |
|------|------|--------------|
| `<out>_resistances.out`, `<out>_resistances_3columns.out` | pairwise | always |
| `<out>_node_currents_<i>_<j>.txt`, `<out>_branch_currents_<i>_<j>.txt` | pairwise | `write_cur_maps = True`. Since 6.0 these honour the flag as the raster path always did; before, they were written on every run. (They are also produced when any other map option is on, since every map option enables the per-pair postprocessing; `write_cum_cur_map_only` does not suppress them.) |
| `<out>_node_currents_cum.txt`, `<out>_branch_currents_cum.txt` | pairwise | `write_cur_maps = True` |
| `<out>_voltages_<i>_<j>.txt` | pairwise | `write_volt_maps = True` |
| `<out>_node_currents.txt`, `<out>_branch_currents.txt` | advanced | `write_cur_maps = True` or `write_cum_cur_map_only = True` |
| `<out>_voltages.txt` | advanced | `write_volt_maps = True` |

`write_max_cur_maps` has no effect in network mode. Node current files have two columns (node ID, current); branch current files three (node, node, current); voltage files two (node ID, voltage).

### Current and Voltage Data

Current and voltage data for networks are written in text list formats. Raster voltage and current maps are written in ESRI ASCII raster format (`.asc`), or GeoTIFF if `write_as_tif = True`. ASCII grids are read and written natively; when the input raster carries a projection, it is written alongside each `.asc` output as a `.prj` sidecar, as GIS tools expect.

### Resistance Files

Resistance data are written in both matrix and 3-column formats.

Here are pairwise resistances written to the output directory for the eight neighbor case (using per-cell resistances and average resistances for cell connection calculations). The first row and column contain the focal node IDs:

      0          1            2            3
      1          0            11.93688471  15.03634473
      2          11.93688471  0            11.57640568
      3          15.03634473  11.57640568  0

Here are pairwise resistances written to the output directory for the four neighbor case, in which focal node 3 was completely isolated (-1 indicates infinite resistance):

      0          1            2            3
      1          0            33.55792693  -1
      2          33.55792693  0            -1
      3          -1           -1           0

For convenience, resistances are also written to a separate file in a 3-column format, e.g.:

      1      2       33.55792693
      1      3       -1
      2      3       -1
