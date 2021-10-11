# Options and Flags

To start the user interface using Windows, run Circuitscape as you would any other installed program. In Mac OS X, double-click on Circuitscape.app. The user interface shown in the Introduction section above will appear. You can also call Circuitscape from the Circuitscape for ArcGIS Toolbox, available from the Circuitscape website, or from the command line.

## Step 1: Choose your input data type

The first step is to choose whether you will be analyzing network or raster data.

## Step 2: Choose a modeling mode

As described above, Circuitscape is run in one of four modes. Pairwise and advanced modes are available for both raster and network data types. The one-to-all and all-to-one modes are available for raster data only.

## Raster resistance map or network/graph

The resistance file specifies the ability of each cell in a landscape or link in a network to carry current (See Figs. 5 and 8). File formats are described in the _Input file formats_ section below.

### Data represent conductances instead of resistances

Most users code their network or raster data in terms of resistances (with higher values denoting greater resistance to movement), which is common in connectivity modeling. Check this box if you want to specify conductances instead (conductance is the reciprocal of resistance; higher values indicate greater ease of movement).

Note that zero and infinite values for conductances and resistances represent special cases. Infinite resistances are coded as NODATA values (see file format description in the _Input file formats_ section below) in input resistance grids, or as zero or NODATA in input conductance grids; these are treated as complete barriers, and are disconnected from all other cells. For raster analyses, cells with zero resistance (infinite conductance) can be specified using a separate short-circuit region file as described below.

## Pairwise, one-to-all, and all-to-one mode options

### Focal node location and data type

This file specifies locations of nodes between which effective resistance and current flow are to be calculated (See Figs. 6 and 9). **Each focal node should have a unique positive integer ID.** Files may be text lists specifying coordinates or appropriate raster grid formats. When a grid is used, it must have the same cell size and extent as the resistance grid. The value stored in each grid cell location refers to the focal node ID. Cells that do not contain focal nodes should be coded with NODATA values. When a text list is used, the value field references the focal node ID.

For raster analyses, focal nodes may occur at points (single cells on the resistance grid) or across regions (Fig. 8). For the latter, a single ID would occupymore than one cell in a grid or more than one pair of coordinates in a text list (and falling within more than one cell in the underlying resistance grid). Cells within a single region would then be collapsed into a single node, as they are when short-circuit region files are used (see below). The difference is that a focal region will be "burned in" to the resistance grid only for pairwise calculations that include that focal node. As with short-circuit regions, focal regions need not be made up of contiguous cells. For large grids or large numbers of focal nodes, focal regions may require more computation time. When calculating resistances on large raster grids and not creating voltage or current maps, focal points will run much more quickly.

### Number of parallel processors to use

On Mac OS X and Linux systems, Circuitscape can run iterations in parallel for pairwise mode when focal points, not focal regions, are used. Choose how many processors you would like to devote to Circuitscape runs.

## Advanced mode options

### Current source file

This file specifies locations and strengths, in amps, of current sources (Figs. 7 and 11). Either a raster or a text list may be used. Rasters must have the same cell size, projection, and extent as the resistance grid, and cells that do not contain current sources should be coded with NODATA values. Note: current sources may be positive or negative (i.e., they may inject current into the grid or pull current out. Similarly, grounds may either serve as a sink for current or may contribute current if there are negative current sources in the grid).

### Ground point file

This file specifies locations of ground nodes and resistances or conductances of resistors tying them to ground (Figs. 7 and 11). Either a raster or a text list may be used. Rasters must have the same cell size, projection, and extent as the resistance grid, and cells that do not contain grounds should be coded with NODATA values. Note that if a direct (R = 0) ground connection conflicts with a current source, the ground will be removed unless the 'remove source' option in the Options Window is chosen.

### Data represent conductances instead of resistances to ground

The default (unchecked) setting is to specify resistances to ground. Checking this box means that your ground point file specifies connections to ground in terms of conductance instead. To tie cells directly to ground, use resistances as the data type and set values in the corresponding ground point file to zero.

## Output options

### Base output filename

Choose a directory path and base file name for output files. Resistances, current maps, voltage maps, and configuration files (which save user interface settings and have a .ini extension) will all use this base name, along with appropriate suffixes and extensions.

### Create current maps

When checked, current maps will be generated for every pair of focal nodes in the pairwise mode, or for the current source and ground configuration specified in the advanced mode. Current maps have the same dimension as the original input files, with values at each node (cell) representing the amount of current flowing through the node. In the pairwise mode, a current map file will be created for each focal node pair, and a cumulative (additive) file will be also written. (Note that for a given pair of focal nodes, current maps are identical regardless of which node is the source and which is the ground due to symmetry). For the advanced modeling mode, a single map will be written showing current densities at each cell resulting from the current source and ground configurations in the input files. These files can be displayed in a GIS as in Fig. 13\. Such maps can be used to identify areas which contribute most to connectivity between focal points (McRae et al. 2008).

![](https://raw.github.com/Circuitscape/Circuitscape/master/docs/4.0/images/image24.png)

**Fig. 1.** Current map used to predict important connectivity areas between core habitat patches (green polygons, entered as focal regions) for mountain lions. Warmer colors indicate areas with higher current density. "Pinch points," or areas where connectivity is most tenuous, are shown in yellow. Quantile classification schemes or "histogram equalize" stretches tend to work well for current maps when using ArcGIS. Research Collaborators: Brett Dickson and Rick Hopkins, Live Oak Associates.

### Create voltage maps

For the pairwise modeling mode, voltage maps give node voltages that would be observed for each focal node pair if one node were connected to a 1 amp current source and the other to ground. For the advanced modeling mode, voltage maps show voltages at each cell resulting from the current source and ground configurations in the input files.

### Calculation Options

#### Connect raster cells to FOUR neighbors instead of EIGHT

For raster operations, Circuitscape creates a graph (network) by connecting cells to their four (Fig. 4) or eight immediate neighbors. The default is eight (four cardinal and four diagonal neighbors), but check this box if you want to connect cells to their four cardinal neighbors only.

#### Use average conductance instead of resistance for connections between cells

For raster operations, this choice determines whether cells are connected by their average resistance or by their average conductance. Most users will want the default (unchecked).

The distinction is particularly important when connecting cells with zero or infinite values. When average resistances are used, first-order neighbors connected by resistors with resistance given by: _Rab_ = (_Ra_ + _Rb_) / 2, and second-order (diagonal) neighbors are connected by resistors with resistance given by: _Rab_ = sqrt(2) * (_Ra_ + _Rb_) / 2, where _Ra_ and _Rb_ are the resistances of the neighboring cells. When average conductances are used, first-order neighbors connected by resistors with conductance (the reciprocal of resistance) given by: _Gab_ = (_Ga_ + _Gb_) / 2, and second-order (diagonal) neighbors are connected by resistors with resistance given by: _Gab_ = (_Ga_ + _Gb_) / (2 * sqrt(2)), where _Ga_ and _Gb_ are the conductances of the neighboring cells. (As noted above, resistance and conductance are reciprocals of each other, i.e., _Gab_ = 1 / _Rab_.)

#### Advanced mode: use unit currents (i=1) for all current sources

All current sources will be set to 1 Amp, regardless of the value specified in the current source input file.

#### Advanced mode: use direct connections to ground (R=0) for all ground points

All ground cells will be tied directly to ground, regardless of the value specified in the input ground file.

#### Advanced mode: when a source and ground are at the same node:

Whenever a cell is connected both to a current source and to ground, this choice will determine whether the source is removed, the ground is removed, both are removed, or both are retained. For the latter, if a source is tied directly to ground (i.e., with zero resistance), the ground connection will be removed.

### Mapping Options

#### Write maximum of current maps

In pairwise, one-to-all, and all-to-one modes, current maps are created for every iteration. By default, Circuitscape will also write a cumulative map showing the sum of values at each node or grid cell across all iterations. If this option is checked, an an extra map that shows the maximum current value at each node or cell across iterations.

#### Write cumulative & max current maps only

Maps of current flow between each pair of focal nodes (or for each focal node in one-to-all and all-to-one modes) will be calculated, but only one summed map of current from all calculations (and a map of maximum values if that option is checked) will be written to disk.

#### Compress output grids

Output ASCII grids are automatically compressed using the gzip file format. This can be useful when many large maps will be written.

#### Log-transform current maps

Values in output current maps will reflect a log10 transform of current densities, which can be useful for visualizing them in some GIS packages (e.g., ArcView 3.X). Cells with zero current will be re-coded with NODATA values.

#### Set focal node currents to zero

When running raster data in pairwise, all-to-one, and one-to all modes, focal nodes will have zero current in output maps when they are activated. For pairwise mode, cumulative maps will still show currents flowing through focal regions that results from other pairs being activated. This helps to show current flowing through a focal region as it moves between other focal regions in cumulative current maps. This current passing through a focal region can give an idea of the importance of the focal region for connecting other focal region pairs (for an example, see Fig. 5 in Dickson et al. 2013).

### Optional Input Files

#### Read raster mask file

When checked, a dialog will open to select a raster mask file. Cells with negative, zero, or NODATA values in the mask will be dropped from the corresponding resistance map (i.e., treated as complete barriers). Positive integer cells will be retained. File should only contain integers and be in raster format.

#### Load a raster short-circuit region map

Short-circuit regions act as areas of zero resistance, essentially providing patches through which current is given a "free ride" as it flows across the landscape. Each short-circuit region should have a unique positive integer identifier; cells within each region are merged into a single node with all other cells in the region, including non-adjacent cells (i.e., regions need not be contiguous). Non-short-circuit-region areas should be stored as NODATA values. The file must have the same cell size and extent as the resistance grid.

#### One-to-all and All-to-One modes: Read source strength file

When checked, a dialog will open to select a text list of focal node IDs and corresponding source strengths. For any focal node in this list, the amount of current injected into that node when it is a source node will be set to the strength specified in the list. All nodes not in the list will default to 1 Amp. This should be in the same file format as the Text List File Format given below, but with two columns (ID followed by source strength). File should have a .txt extension.

#### Read file with focal node pairs to include/exclude

This option allows users to only perform calculations on a subset of focal node pairs. Users can either identify pairs to include in calculations, or pairs to exclude, as specified in the first line of the file.

This affects all modes except the Advanced Mode. Files should be in tab-delimited text with a .txt extension. See formatting information in the _Input file formats_ section below.


## Input raster format

Raster input maps should be stored in Arc/Info ASCII grid of GeoTIFF format, as exported by standard GIS packages. For focal nodes, the value stored in each grid location refers to the focal node ID, and a single ID can occupy more than one cell (IDs must be positive integers). For current sources, the grid value specifies the source strength in amps. For grounds, the grid value specifies either the resistance or conductance of the resistor tying each ground node to ground, as specified in the Options window.

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

Below is a 10 x10 resistance map. Cells with infinite resistance are assigned NODATA values (-9999):

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

Below is a 10 x 10 focal region map. Here, groups of cells have been coded as focal regions- so these will be treated as "core area polygons" to be connected in circuit analyses. All cells within each focal region will be collapsed into a single node (even the non-contiguous cell in region #1) when that region is activated in pairwise, one-to-all, or all-to-one analyses. This format is identical to the short-circuit region file format.

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

## Text list file format

For network/graph operations, resistor networks, focal nodes, current sources, and grounds should be stored as text lists (saved with a ".txt" extension). To specify a network of resistors, three columns are used. The first and second columns give the node IDs being connected by a resistor, and the third column gives the resistance value. For example, the simple circuit:

![](https://raw.github.com/Circuitscape/Circuitscape/master/docs/4.0/images/SimpleNetworkWithNumbers2.png)

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

![](https://raw.github.com/Circuitscape/Circuitscape/master/docs/4.0/images/AdvancedNetwork.png)

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

For focal nodes, the value field references the focal node ID; values must be positive integers, and a single ID can occupy more than one pair of coordinates (and more than one cell in the underlying resistance grid). For current sources, the value field references the source strength in amps. For grounds, the value field references either the resistance or conductance of the resistor tying each ground node to ground, as set in the Options window.

## Include/exclude file format

This file will be loaded when the 'Read file with focal node pairs to include/exclude' option is checked, and affects all modes except the advanced mode. There are two file formats that can be used. The first is the simplest, and gives a list of pairs to include in calculations, or pairs to exclude, as specified in the first line of the file. For example, if there are five focal nodes, numbered 1-5, and the following list is entered, only pairs (1,2), (1,3), and (1,5) will be analyzed:

    mode    include
    1        2
    1        3
    1        5

Similarly, if the first line in the above file read:

    mode     exclude

all pairs except (1,2), (1,3), and (1,5) would be analyzed.

The second method uses a matrix identifying which pairs of focal nodes to connect. The file specifies minimum and maximum values in the matrix to consider a pair connected. This method can be useful when used with a distance matrix to only run analyses between points separated by a minimum distance, or by a distance equal to or less than a maximum distance. Note: any focal node not in the matrix will be dropped from analyses. Entries on the diagonal are ignored. For example, in the following matrix, only pairs with entries between 2 and 50 are connected. Pairs (1,2), (2,4), and (3,4) will not be analyzed.  
Focal node 5 will be dropped entirely:

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

# 9\. Output files

## Current and voltage data

Current and voltage data for networks will be written in text list formats.

When using the ArcGIS toolbox, current and voltage maps will be written in the raster format chosen by the user.

When not using the ArcGIS toolbox, raster voltage and current maps are written using the ASCII raster format described above.

## Resistance files

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