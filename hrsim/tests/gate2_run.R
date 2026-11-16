# gate2_run.R
#
# Part of the tests of the hrsim package, the simulation engine for Berryhill
# and Andersson (2026), "Unstable Heteroresistance as a Stepping Stone to
# Stable Resistance via Compensatory Mutations".
#
# What it does
#   Runs stochastic replicates of the Scenario 1 condition with the compiled
#   engine, classifies each one, and appends one line per seed to a CSV file.
#   It is not called by run_all_gates.R.
#
#   The file is written one seed at a time, as soon as that seed finishes,
#   and seeds already present in it are skipped, so the script can be stopped
#   and run again until the target number of seeds is reached. A timer is
#   set just before the loop over the seeds, and a new seed is started only
#   while the time elapsed on it does not exceed the cutoff option, so a
#   single call does a limited amount of work.
#
#   By default the output file is the stored fixture that gate2_check.R
#   reads, tests/fixtures/gate2_results.csv, which this script then extends.
#   Set hrsim.gate2_results to write elsewhere. The second stored file,
#   tests/fixtures/gate2_screen48.csv, is the output of this script for seeds
#   1 to 9 with hrsim.gate2_horizon = 48.
#
# Classification
#   Each run is classified on its transfer rows (cycle_end == 1): `extinction`
#   if the total density on the last of them is below 1; otherwise `basin1`
#   (climb) if the mean copy number reached 5 on any of them; otherwise
#   `basin2` (low-copy). A run with no transfer row is labelled `undecided`.
#
# Tolerance
#   None here; the comparison is made in gate2_check.R.
#
# Seeds
#   1 to the target, in increasing order, each passed to run_sim() as the
#   seed of that replicate. A replicate's random stream depends only on its
#   own seed.
#
# Options
#   hrsim.gate2_results         output CSV (default as above)
#   hrsim.gate2_target          highest seed to run (default 200)
#   hrsim.gate2_horizon         simulated hours per replicate (default 300)
#   hrsim.gate2_start_cutoff_s  no new seed is started once more than this
#                               many seconds have passed on the timer
#                               (default 60)
#
# Run from the package root with hrsim installed.

suppressMessages(library(hrsim))
ref_dir <- getOption("hrsim.ref", "tests/ref")
source(file.path(ref_dir, "load_canonical.R"))

RESULTS      <- getOption("hrsim.gate2_results", "tests/fixtures/gate2_results.csv")
HORIZON      <- getOption("hrsim.gate2_horizon", 300)
TARGET       <- as.integer(getOption("hrsim.gate2_target", 200L))
START_CUTOFF <- getOption("hrsim.gate2_start_cutoff_s", 60)

pop_cols <- c(paste0("R", 0:9), paste0("R", 0:9, "C"))
cn <- copy_number

# Classifies one run from its transfer rows. Returns the label, the peak mean
# copy number, the total density on the last transfer row, and the index of
# the first transfer row whose mean copy number is at least 5 (0 if none).
classify <- function(d) {
  tr <- d[d$cycle_end == 1, , drop = FALSE]
  if (nrow(tr) == 0) return(list(basin = "undecided", mcn_peak = NA, total_end = NA, first5 = 0L))
  P   <- as.matrix(tr[, pop_cols, drop = FALSE]); tot <- rowSums(P)
  mcn <- ifelse(tot > 0, as.numeric(P %*% cn) / tot, 0)
  total_end <- tot[length(tot)]; mcn_peak <- max(mcn)
  first5 <- which(mcn >= 5)[1]; if (is.na(first5)) first5 <- 0L
  basin <- if (total_end < 1) "extinction" else if (mcn_peak >= 5) "basin1" else "basin2"
  list(basin = basin, mcn_peak = mcn_peak, total_end = total_end, first5 = first5)
}

done <- if (file.exists(RESULTS)) read.csv(RESULTS)$seed else integer(0)
if (!file.exists(RESULTS))
  cat("seed,basin,mcn_peak,total_end,first5_transfer,n_transfer\n", file = RESULTS)
remaining <- setdiff(1:TARGET, done)

t0 <- Sys.time(); ran <- 0L
for (s in remaining) {
  if (as.numeric(difftime(Sys.time(), t0, units = "secs")) > START_CUTOFF) break
  d  <- run_sim(p, transitions, y0, t_end = HORIZON, seed = s, stochastic = TRUE, n_out = -1L)
  cl <- classify(d)
  cat(sprintf("%d,%s,%.4f,%.6g,%d,%d\n",
              s, cl$basin, cl$mcn_peak, cl$total_end, cl$first5, sum(d$cycle_end == 1)),
      file = RESULTS, append = TRUE)
  ran <- ran + 1L
  cat(sprintf("  seed %d -> %-10s mcn_peak=%.2f  first>=5 @transfer %d  total_end=%.4g\n",
              s, cl$basin, cl$mcn_peak, cl$first5, cl$total_end))
}
total_done <- length(read.csv(RESULTS)$seed)
cat(sprintf("Gate2 runner: +%d this call; %d/%d done.\n", ran, total_done, TARGET))
