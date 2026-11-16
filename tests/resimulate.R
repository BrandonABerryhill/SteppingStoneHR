# resimulate.R
#
# A small end-to-end test of the re-simulation mode of SteppingStone.Rmd
# (`resimulate: true`). It runs a cut-down set of batches through the same
# functions the notebook calls in that mode, those of R/simulation/08_notebook.R
# and the modules they call, and compares the replicates it obtains with
# results shipped in results/. Part of the code for "Unstable Heteroresistance
# as a Stepping Stone to Stable Resistance via Compensatory Mutations"
# (Berryhill and Andersson).
#
# What it runs, with every setting at the notebook's value except the numbers
# of seeds:
#   - the comparison of the R reference and the compiled engine with the
#     random layer switched off, to 300 h; the test halts if the largest
#     relative error exceeds the notebook's tolerance;
#   - Scenario 1, Scenario 2 and control frequency batches on the compiled
#     engine, seeds 1 to 20 each, to 300 h, transfer rows kept;
#   - the R reference engine for Scenario 1, seeds 1 and 2, to 300 h. These
#     two replicates are run and classified but not compared with shipped
#     results; the test checks only that both were returned and that the
#     compiled engine is the active engine again afterwards;
#   - one Scenario 1 run to 3,000 h on the compiled engine, seed 1;
#   - parameter sets 1 and 2 of the full 300-set design of the sensitivity
#     sweep, 100 seeds each.
#
# What it compares. Each comparison applies identical() to every column the
# two data frames share, after checking that they have the same number of
# rows:
#   - Scenario 1, seeds 1 to 20: with the transfer rows of the same seed in
#     results/S1/S1_traj300/sims.rds;
#   - Scenario 2, seeds 1 to 20: with the transfer rows before 300 h (the 12
#     transfers from 24 h to 288 h) of the same seed in
#     results/S2/S2_long/sims.rds;
#   - control, seeds 1 to 20: with the same seed in
#     results/omega0/omega0_freq_cpp/sims.rds;
#   - the 3,000 h run: with seed 1 in results/S1/S1_long/sims.rds;
#   - sweep sets 1 and 2: with rows 1 and 2 of
#     results/sensitivity/sweep_S1/sweep.rds, on the six parameter columns
#     and the three frequency columns.
# It prints one PASS or FAIL line per comparison and a summary line, and exits
# with status 1 if any comparison fails or the engine comparison exceeds the
# tolerance, and with status 0 otherwise. An error in any step also ends the
# script with a non-zero status.
#
# How the settings are obtained. The settings of the notebook (conditions,
# horizons, output grids, sweep ranges, master seed and the tolerance of the
# engine comparison) are not repeated here. The script reads SteppingStone.Rmd,
# takes the code of its chunk labelled setup-settings and evaluates it, so a
# change to that chunk reaches this test as well. The script stops if it finds
# no such chunk, or more than one, or if the chunk does not define the
# settings used below. The notebook's parameters (the shipped results folder
# and the default number of workers, the number of cores minus 2) are read
# from its YAML header with knitr::knit_params(), which evaluates the `!r`
# expression of n_workers as knitting would.
#
# How to run it. From the repository root, with the hrsim package installed
# (from hrsim_0.1.0.tar.gz) and the packages future, furrr, lhs and knitr:
#   Rscript tests/resimulate.R
#   Rscript tests/resimulate.R --n-workers=6
# The optional argument sets the number of parallel workers; the number of
# workers does not change any result. Every batch is written to a new folder
# made with tempfile() under the R session's temporary directory, which R
# removes when the session ends; nothing is written into the repository.
#
# How long it takes. Timed twice with 12 workers, the whole test took about
# 12 minutes each time. Most of the time goes to the R reference engine (the
# engine comparison and the two R-reference replicates) and to the sweep sets.
#
# Seeds. Replicate seeds are passed to the engine and govern all of that
# replicate's random draws; each replicate's random stream depends only on
# its own seed, not on the number of workers. The seeds used are 1 to 20 for
# the frequency batches, 1 and 2 for the R reference, 1 for the 3,000 h run
# and 1 to 100 for each sweep set. The sweep design is drawn, as in the
# notebook, with set.seed(MASTER_SEED + 7919) immediately before the Latin
# hypercube sample; the settings chunk calls set.seed(MASTER_SEED) first,
# which keeps the session on R's default generator. The engine comparison
# passes seed 1, which has no effect without the random layer.
#
# Defines no functions for use elsewhere; the helpers below serve this script
# only.

started <- Sys.time()

# Repository and arguments ----

if (!file.exists("SteppingStone.Rmd") || !dir.exists(file.path("R", "simulation")))
  stop("Run this script from the repository root, the folder that holds ",
       "SteppingStone.Rmd.", call. = FALSE)

needed <- c("hrsim", "future", "furrr", "lhs", "knitr")
absent <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(absent) > 0)
  stop("Install these packages before running the test: ",
       paste(absent, collapse = ", "), call. = FALSE)

rmd <- readLines("SteppingStone.Rmd", warn = FALSE)
notebook_params <- lapply(knitr::knit_params(rmd), function(p) p$value)

args <- commandArgs(trailingOnly = TRUE)
workers_arg <- grep("^--n-workers=", args, value = TRUE)
N_WORKERS <- if (length(workers_arg) == 1) {
  as.integer(sub("^--n-workers=", "", workers_arg))
} else {
  as.integer(notebook_params$n_workers)
}
if (length(args) > length(workers_arg) || length(N_WORKERS) != 1 ||
    is.na(N_WORKERS) || N_WORKERS < 1)
  stop("The only argument accepted is --n-workers=N with N a positive ",
       "integer.", call. = FALSE)

# Settings of the notebook ----

# The chunk starts at the line "```{r setup-settings}" (possibly followed by
# chunk options) and ends at the next line that is a closing fence.
chunk_start <- grep("^```\\{r setup-settings[,} ]", rmd)
if (length(chunk_start) != 1)
  stop("Expected one chunk labelled setup-settings in SteppingStone.Rmd, ",
       "found ", length(chunk_start), ".", call. = FALSE)
chunk_end <- chunk_start + which(rmd[-seq_len(chunk_start)] == "```")[1]
if (is.na(chunk_end))
  stop("The setup-settings chunk of SteppingStone.Rmd has no closing fence.",
       call. = FALSE)
eval(parse(text = rmd[(chunk_start + 1):(chunk_end - 1)]), envir = globalenv())

settings <- c("CONDITIONS", "H_CLASS", "NOUT_CLASS", "H_LONG", "NOUT_LONG",
              "SWEEP_SETS", "SWEEP_SEEDS", "SWEEP_RANGES", "SWEEP_LOG",
              "MASTER_SEED", "DET_RELERR_TOL")
missing_settings <- settings[!vapply(settings, exists, logical(1))]
if (length(missing_settings) > 0)
  stop("The setup-settings chunk did not define: ",
       paste(missing_settings, collapse = ", "), call. = FALSE)

# Modules ----

# 00_config.R keeps RESULTS_DIR and FIGURE_DIR when they are assigned before
# it is sourced, so every batch, the manifest and the sweep tables go to the
# temporary folder.
shipped_dir <- notebook_params$results_dir
RESULTS_DIR <- tempfile("smoke_resimulate_")
FIGURE_DIR  <- file.path(RESULTS_DIR, "figures")
USE_ENGINE  <- "cpp"
for (module in c("00_config.R", "01_parameters.R", "02_transitions.R",
                 "03_engine_R.R", "04_classify.R", "05_run.R",
                 "06_sensitivity.R", "07_engines.R", "08_notebook.R"))
  source(file.path("R", "simulation", module))

cat(sprintf("Re-simulation smoke test: %d workers, batches written to %s\n",
            N_WORKERS, RESULTS_DIR))

# Reporting ----

results <- character(0)

report <- function(ok, label, detail = NULL) {
  status <- if (isTRUE(ok)) "PASS" else "FAIL"
  results[label] <<- status
  cat(sprintf("%s  %s\n", status, label))
  if (!isTRUE(ok) && length(detail) > 0) cat(paste0("      ", detail, "\n"), sep = "")
}

elapsed <- function(t0) as.numeric(difftime(Sys.time(), t0, units = "secs"))

#' Compare two replicates column by column
#'
#' @param new Data frame of a re-simulated replicate.
#' @param ref Data frame of the shipped replicate it is compared with.
#' @return Character vector describing each difference found: a missing
#'   replicate, a different number of rows, no shared column, or the names of
#'   the shared columns that are not identical(). Empty when the two agree.
replicate_differences <- function(new, ref) {
  if (is.null(new)) return("re-simulated replicate missing")
  if (is.null(ref)) return("shipped replicate missing")
  if (nrow(new) != nrow(ref))
    return(sprintf("%d rows against %d shipped", nrow(new), nrow(ref)))
  shared <- intersect(names(new), names(ref))
  if (length(shared) == 0) return("no shared column")
  differ <- shared[!vapply(shared, function(k) identical(new[[k]], ref[[k]]),
                           logical(1))]
  if (length(differ) == 0) return(character(0))
  sprintf("columns that differ: %s", paste(differ, collapse = ", "))
}

transfer_rows <- function(d) d[d$cycle_end == 1, , drop = FALSE]

read_shipped <- function(...) readRDS(file.path(shipped_dir, ...))

# Engine comparison ----

t0 <- Sys.time()
gate <- tryCatch(
  engine_gate(t_end = H_CLASS, n_out = NOUT_CLASS, tol = DET_RELERR_TOL,
              halt = TRUE,
              record = file.path(RESULTS_DIR, "engine_verification.rds")),
  error = function(e) e)
if (inherits(gate, "error")) {
  report(FALSE, sprintf("engine comparison, 0 to %g h, random layer off", H_CLASS),
         conditionMessage(gate))
  cat("Halted after the engine comparison: 1 of 1 checks failed.\n")
  quit(save = "no", status = 1)
}
report(TRUE, sprintf(paste("engine comparison, 0 to %g h, random layer off:",
                           "maximum relative error %.3e over %d columns and",
                           "%d rows (tolerance %.1e)"),
                     gate$t_end, gate$relerr, gate$n_columns, gate$n_rows,
                     gate$tol))
cat(sprintf("      %.0f s\n", elapsed(t0)))

# Frequency batches, compiled engine ----

seeds_freq <- 1:20
t0 <- Sys.time()
freq <- list()
for (sc in CONDITIONS)
  freq[[sc]] <- load_or_run_batch(
    sc, paste0(sc, "_freq_cpp"), t_end = H_CLASS, seeds = seeds_freq,
    keep = "cycle_end", engine = "cpp", resimulate = TRUE,
    n_workers = N_WORKERS)
cat(sprintf("Frequency batches run: %s, seeds 1 to %d, %.0f s\n",
            paste(CONDITIONS, collapse = ", "), max(seeds_freq), elapsed(t0)))

ship_s1 <- read_shipped("S1", "S1_traj300", "sims.rds")
for (s in as.character(seeds_freq)) {
  ref <- ship_s1[[s]]
  if (!is.null(ref)) ref <- transfer_rows(ref)
  d <- replicate_differences(freq$S1$sims[[s]], ref)
  report(length(d) == 0,
         sprintf("S1 seed %s against the transfer rows of shipped S1_traj300", s), d)
}

ship_s2 <- read_shipped("S2", "S2_long", "sims.rds")
for (s in as.character(seeds_freq)) {
  ref <- ship_s2[[s]]
  if (!is.null(ref)) {
    ref <- transfer_rows(ref)
    ref <- ref[ref$time < H_CLASS, , drop = FALSE]
  }
  d <- replicate_differences(freq$S2$sims[[s]], ref)
  report(length(d) == 0,
         sprintf("S2 seed %s against the transfer rows to 288 h of shipped S2_long", s), d)
}
rm(ship_s2)

ship_o0 <- read_shipped("omega0", "omega0_freq_cpp", "sims.rds")
for (s in as.character(seeds_freq)) {
  d <- replicate_differences(freq$omega0$sims[[s]], ship_o0[[s]])
  report(length(d) == 0,
         sprintf("omega0 seed %s against shipped omega0_freq_cpp", s), d)
}
rm(ship_o0)

# R reference engine ----

t0 <- Sys.time()
ref_R <- load_or_run_batch(
  "S1", "S1_freq_R", t_end = H_CLASS, seeds = 1:2, keep = "cycle_end",
  engine = "R", resimulate = TRUE, n_workers = N_WORKERS)
report(identical(names(ref_R$sims), c("1", "2")) &&
         identical(run_model, run_model_cpp),
       "S1 R reference, seeds 1 and 2, run (not compared with shipped results)",
       c(sprintf("replicates returned: %s", paste(names(ref_R$sims), collapse = ", ")),
         sprintf("compiled engine active again: %s",
                 identical(run_model, run_model_cpp))))
cat(sprintf("      outcomes: %s; %.0f s\n",
            paste(sprintf("seed %s %s", ref_R$cls$seed, ref_R$cls$outcome),
                  collapse = ", "), elapsed(t0)))

# One 3,000 h run ----

t0 <- Sys.time()
long1 <- load_or_run_batch(
  "S1", "S1_long", t_end = H_LONG, seeds = 1, keep = "all",
  n_out = NOUT_LONG, engine = "cpp", resimulate = TRUE, n_workers = N_WORKERS)
d <- replicate_differences(long1$sims[["1"]],
                           read_shipped("S1", "S1_long", "sims.rds")[["1"]])
report(length(d) == 0,
       sprintf("S1 seed 1 to %g h against seed 1 of shipped S1_long", H_LONG), d)
cat(sprintf("      %.0f s\n", elapsed(t0)))

# Sweep sets 1 and 2 ----

t0 <- Sys.time()
sweep_run <- sensitivity_sweep(
  n_sets = SWEEP_SETS, seeds_per_set = SWEEP_SEEDS,
  ranges = SWEEP_RANGES, log = SWEEP_LOG,
  design_seed = MASTER_SEED + 7919L, resimulate = TRUE,
  n_workers = N_WORKERS, first_sets = 2)
sweep      <- sweep_run$sweep
ship_sweep <- read_shipped("sensitivity", "sweep_S1", "sweep.rds")
sweep_cols <- c(names(SWEEP_RANGES), "f_climb", "f_lowcopy", "f_extinct")
show_value <- function(x) {
  if (length(x) == 0) "missing" else format(x, digits = 17)
}
for (i in 1:2) {
  differ <- sweep_cols[!vapply(sweep_cols, function(k)
    identical(sweep[[k]][i], ship_sweep[[k]][i]), logical(1))]
  detail <- vapply(differ, function(k)
    sprintf("%s: %s against %s shipped", k, show_value(sweep[[k]][i]),
            show_value(ship_sweep[[k]][i])), character(1))
  report(length(differ) == 0,
         sprintf(paste("sweep set %d, %d seeds, against row %d of the shipped",
                       "sweep table (6 parameters, 3 frequencies)"),
                 i, SWEEP_SEEDS, i), detail)
}
cat(sprintf("      %.0f s\n", elapsed(t0)))

# Summary ----

n_fail <- sum(results == "FAIL")
cat(sprintf("%s: %d of %d checks passed, %d failed; %.1f min\n",
            if (n_fail == 0) "PASS" else "FAIL",
            length(results) - n_fail, length(results), n_fail,
            elapsed(started) / 60))
quit(save = "no", status = if (n_fail == 0) 0 else 1)
