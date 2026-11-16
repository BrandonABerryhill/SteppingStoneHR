# gate2_check.R
#
# Part of the tests of the hrsim package, the simulation engine for Berryhill
# and Andersson (2026), "Unstable Heteroresistance as a Stepping Stone to
# Stable Resistance via Compensatory Mutations".
#
# What it verifies
#   Reads the outcomes of stochastic Scenario 1 runs of the compiled engine
#   that gate2_run.R stored in tests/fixtures/gate2_results.csv, and compares
#   the frequency of each outcome with the reference percentages held in
#   `ref` below. This script runs no simulation. The outcome labels are those
#   written by gate2_run.R: `basin1` is a climb to high copy number, `basin2`
#   is low-copy survival and `extinction` is extinction.
#
# Tolerance
#   For each outcome the 95% confidence interval of the observed proportion
#   is taken from prop.test() (a Wilson score interval with continuity
#   correction). The outcome is flagged as overlapping when its reference
#   percentage lies inside that interval. The interval is wide when the file
#   holds few seeds, so a small sample gives a weak test.
#
# Seeds
#   Those listed in the `seed` column of the results file. The stored file
#   holds seeds 1 to 9; its `n_transfer` column gives the number of transfer
#   rows on which each seed was classified: 12 for seeds 1, 2 and 7, and 1
#   for the other six, whose lines equal those of
#   tests/fixtures/gate2_screen48.csv.
#
# Run from the package root.

RESULTS <- "tests/fixtures/gate2_results.csv"
if (!file.exists(RESULTS)) { cat("Gate 2: no results yet (run gate2_run.R)\n"); quit(save = "no") }
r <- read.csv(RESULTS); n <- nrow(r)
ref <- c(basin1 = 83, basin2 = 13, extinction = 4)
tab <- table(factor(r$basin, levels = names(ref)))

cat(sprintf("== Gate 2: stochastic basin frequencies (n=%d seeds, 300h, S1) ==\n", n))
all_overlap <- TRUE
for (b in names(ref)) {
  k  <- as.integer(tab[[b]])
  ci <- tryCatch(prop.test(k, n)$conf.int * 100, error = function(e) c(NA, NA))
  ov <- isTRUE(ref[[b]] >= ci[1] && ref[[b]] <= ci[2])
  all_overlap <- all_overlap && ov
  cat(sprintf("  %-11s obs=%5.1f%% (%2d/%d)  95%% CI [%5.1f, %5.1f]  ref=%2g%%  overlap=%s\n",
              b, 100 * k / n, k, n, ci[1], ci[2], ref[[b]], ov))
}
# The line printed here is for information and does not enter the result
# line. first5_transfer is the index of the first transfer row on which the
# mean copy number was at least 5 (0 if none), as written by gate2_run.R.
b1 <- r[r$basin == "basin1" & r$first5_transfer > 0, "first5_transfer"]
if (length(b1)) cat(sprintf("  basin1 resolved (mcn>=5) by transfer: median %g, max %g (of %d transfers)\n",
                            median(b1), max(b1), suppressWarnings(max(r$n_transfer))))
cat(sprintf("\nGATE 2 (%s sample): %s\n",
            if (n >= 200) "full" else sprintf("n=%d", n),
            if (all_overlap) "consistent with 83/13/4" else "CHECK — see overlap flags"))
