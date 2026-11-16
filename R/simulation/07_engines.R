# 07_engines.R
#
# Binds the simulation engines to the names the other modules call. After this
# file has been sourced:
#   run_model_R    is the R reference engine of 03_engine_R.R, always;
#   run_model      is the active engine, the one run_single() calls: the
#                  compiled engine hrsim::run_sim() when USE_ENGINE is "cpp",
#                  otherwise the R reference;
#   run_model_cpp  is the compiled engine, defined only when USE_ENGINE is
#                  "cpp".
# The R reference is always available so that the two engines can be compared
# and so that a batch can be run on it: assign run_model_R to the global
# run_model before the batch and put the previous value back afterwards.
#
# USE_ENGINE, ENGINE_R and PIPE_ROOT come from 00_config.R, so this file is
# sourced after it. With USE_ENGINE equal to "cpp" the hrsim package must be
# installed; assign USE_ENGINE <- "R" before sourcing this file to work without
# it.
#
# Part of the code for "Unstable Heteroresistance as a Stepping Stone to Stable
# Resistance via Compensatory Mutations" (Berryhill and Andersson).
#
# Defines: run_model_R, run_model and, with the compiled engine, run_model_cpp.

# R reference engine ----

# The file is sourced into an environment of its own, so that sourcing it does
# not assign the global run_model.
.R_env <- new.env()
source(file.path(PIPE_ROOT, ENGINE_R), local = .R_env)
run_model_R <- .R_env$run_model
rm(.R_env)

# Active engine ----

if (USE_ENGINE == "cpp") {
  if (!requireNamespace("hrsim", quietly = TRUE))
    stop("USE_ENGINE = 'cpp' but the hrsim package is not installed. See README.md.",
         call. = FALSE)
  suppressMessages(library(hrsim))
  run_model     <- hrsim::run_sim
  run_model_cpp <- hrsim::run_sim
} else {
  run_model <- run_model_R
}
