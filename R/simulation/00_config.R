# 00_config.R
#
# Settings shared by the simulation modules: the folders that results and
# figures are written to, the rule that classifies the outcome of a replicate,
# the simulation horizons, the default number of parallel workers and the
# choice of engine. All times are in hours.
#
# Every setting written as `if (!exists("NAME")) NAME <- value` is a default: a
# script that assigns NAME before sourcing this file keeps its own value. The
# results folder is chosen this way. Assign RESULTS_DIR before sourcing, either
# the folder of shipped results or a separate folder for a re-simulation. A
# re-simulation needs a folder of its own because run_batch() resumes from any
# batch it finds under RESULTS_DIR and would otherwise run nothing. MANIFEST is
# derived from RESULTS_DIR when this file is sourced. run_batch() reads
# RESULTS_DIR at each call and MANIFEST whenever it appends a manifest row, so
# a script that changes RESULTS_DIR afterwards must reassign MANIFEST as well.
# Paths may be relative to the working directory. Sourcing this file creates
# RESULTS_DIR and FIGURE_DIR if they do not exist.
#
# Part of the code for "Unstable Heteroresistance as a Stepping Stone to Stable
# Resistance via Compensatory Mutations" (Berryhill and Andersson).
#
# Defines: PIPE_ROOT, RESULTS_DIR, FIGURE_DIR, MANIFEST, CLASSIFY_HORIZON,
# CLIMB_THRESHOLD, EXTINCT_BELOW, SHORT_HORIZON, LONG_HORIZON, N_WORKERS,
# USE_ENGINE, ENGINE_R, CODE_VERSION and new_batch_id().

# Paths ----

# PIPE_ROOT is the folder that contains R/simulation/. ENGINE_R below and the
# default results and figures folders are resolved against it.
if (!exists("PIPE_ROOT")) PIPE_ROOT <- getwd()

# Batches are written to RESULTS_DIR/<scenario>/<tag>/. MANIFEST is a CSV file
# to which run_batch() appends one row per call, unless it is called with
# manifest = FALSE.
if (!exists("RESULTS_DIR")) RESULTS_DIR <- file.path(PIPE_ROOT, "results")
if (!exists("FIGURE_DIR"))  FIGURE_DIR  <- file.path(PIPE_ROOT, "figures")
MANIFEST <- file.path(RESULTS_DIR, "manifest.csv")
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGURE_DIR,  showWarnings = FALSE, recursive = TRUE)

# Classification rule ----

# The outcome frequencies reported in the paper depend on these three values.
# classify_outcome() in 04_classify.R applies the two thresholds.
#
# CLASSIFY_HORIZON: simulated time, in hours, of a replicate that is to be
# classified. A 300 h run has 12 serial transfers and is classified on the 12
# rows recorded at them, the last at 288 h, the classification horizon.
if (!exists("CLASSIFY_HORIZON")) CLASSIFY_HORIZON <- 300
# CLIMB_THRESHOLD: a surviving replicate whose peak mean copy number is at
# least this value is a climb; otherwise it is low-copy.
if (!exists("CLIMB_THRESHOLD"))  CLIMB_THRESHOLD  <- 5
# EXTINCT_BELOW: a replicate whose total density on its last row is below this
# value, in cells per mL, is an extinction.
if (!exists("EXTINCT_BELOW"))    EXTINCT_BELOW    <- 1

# Other horizons ----

# SHORT_HORIZON is one transfer interval, so a run of this length ends before
# its first transfer. LONG_HORIZON is 125 transfer intervals.
if (!exists("SHORT_HORIZON")) SHORT_HORIZON <- 24
if (!exists("LONG_HORIZON"))  LONG_HORIZON  <- 3000

# Parallel workers ----

# Default worker count for scripts that pass N_WORKERS to run_batch() or
# run_sweep(). The number of workers does not change any result (see
# 05_run.R).
if (!exists("N_WORKERS")) N_WORKERS <- 32

# Engine selection ----

# USE_ENGINE is "cpp" for the compiled engine of the hrsim package or "R" for
# the R reference engine. 07_engines.R reads it to bind run_model. ENGINE_R is
# the path of the R reference engine relative to PIPE_ROOT.
if (!exists("USE_ENGINE")) USE_ENGINE <- "cpp"
ENGINE_R <- "R/simulation/03_engine_R.R"

# Run records ----

# Written to the code_version column of every manifest row.
CODE_VERSION <- "pipeline-2.0.0"

#' Time-stamped identifier for an output folder
#'
#' @param tag Character scalar appended after the time stamp; "" for none.
#' @return Character scalar "YYYYMMDD-HHMMSS_<tag>", or the time stamp alone
#'   when `tag` is empty, from the system clock at the time of the call.
new_batch_id <- function(tag = "") {
  ts <- format(Sys.time(), "%Y%m%d-%H%M%S")
  if (nzchar(tag)) paste0(ts, "_", tag) else ts
}
