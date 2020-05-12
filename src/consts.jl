# Various constants used to determine problem type

const RASTER = ["raster", "Raster"]
const PAIRWISE = ["pairwise", "Pairwise"]
const ADVANCED = ["advanced", "Advanced"]
const ONETOALL = ["one-to-all", "one_to_all"]
const ALLTOONE = ["all-to-one", "all_to_one"]
const SINGLE = ["single", "Single"]
const DOUBLE = ["double", "Double"]

# Solver constants
const AMG = ["cg+amg", "amg+cg"]
const CHOLMOD = ["cholmod", "cholesky", "cholfact"]

# Constants used in IO
const AAGRID = 2
const TXTLIST = 3
const PAIRS_AAGRID = 4
const PAIRS_LIST = 5
const TRUELIST = ["True", "true", "1"]

const FILE_TYPE_NPY = 1
const FILE_TYPE_AAGRID = 2
const FILE_TYPE_TXTLIST = 3
const FILE_TYPE_INCL_PAIRS_AAGRID = 4
const FILE_TYPE_INCL_PAIRS = 5
const FILE_TYPE_GEOTIFF = 6

const FILE_HDR_GZIP = "\x1f\x8b\x08"
const FILE_HDR_NPY = "\x93NUMPY"
const FILE_HDR_AAGRID = "ncols"
const FILE_HDR_INCL_PAIRS_AAGRID = "min"
const FILE_HDR_INCL_PAIRS = "mode"

# Constants for logging
const NONE = ["NONE", "None", "none"]
const INFO = ["INFO", "info", "Info"]
const DEBUG = ["DEBUG", "debug", "Debug"]
const WARNING = ["WARNING", "WARN", "warning", "warn", "Warning", "Warn"]
const CRITICAL = ["CRITICAL", "critical", "Critical"]
