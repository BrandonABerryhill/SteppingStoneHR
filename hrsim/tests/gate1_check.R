# gate1_check.R
#
# Part of the tests of the hrsim package, the simulation engine for Berryhill
# and Andersson (2026), "Unstable Heteroresistance as a Stepping Stone to
# Stable Resistance via Compensatory Mutations".
#
# What it verifies
#   With the Monte Carlo layer disabled (stochastic = FALSE) the compiled
#   engine is a deterministic Euler integrator. Its output for the Scenario 1
#   condition is compared with trajectories of the R reference engine stored
#   in tests/fixtures/det_24h.csv (401 rows, no transfer) and
#   tests/fixtures/det_300h.csv (413 rows, of which 12 are transfer rows).
#   Both engines record on the n_out = 400 grid. The 24 columns other than
#   time (the 20 states, r, A, B and cycle_end) are compared cell by cell
#   under the tolerance below. The time column is not part of that test: the
#   two time columns only have to agree within 1e-6 h on every row, and the
#   comparison stops with an error if they do not.
#
# Tolerance
#   A cell passes when |C++ - R| <= 1e-12 + rtol * |R|, with rtol = 1e-9 over
#   24 h and rtol = 1e-6 over 300 h. The bounds are set several orders of
#   magnitude above the relative errors between the two engines reported in
#   the Methods of the paper, so that the test checks the arithmetic and is
#   not sensitive to the last bits; the 300 h bound is the looser one because
#   rounding differences can accumulate over 3e7 steps where the 24 h run
#   takes 2.4e6. The largest relative error found is printed.
#
# Seeds
#   None. seed = 1L is passed but no random number is drawn when
#   stochastic = FALSE.
#
# Run from the package root with hrsim installed.

suppressMessages(library(hrsim))
ref_dir <- getOption("hrsim.ref", "tests/ref")
source(file.path(ref_dir, "load_canonical.R"))
source("tests/helpers_compare.R")

read_fix <- function(f) as.matrix(read.csv(f, check.names = FALSE))

# Runs the compiled engine to t_end hours and compares it with one fixture.
# Returns TRUE when every cell outside the time column is within tolerance.
gate1_one <- function(t_end, fixture, rtol) {
  ref  <- read_fix(fixture)
  test <- as.matrix(run_sim(p, transitions, y0, t_end = t_end, seed = 1L,
                            stochastic = FALSE, n_out = 400L))
  cmp <- compare_trajectories(ref, test, rtol = rtol, atol = 1e-12)
  cat(sprintf("Gate 1 @ %dh: rows R=%d C++=%d | pass=%s n_fail=%d max_relerr=%.3e\n",
              t_end, nrow(ref), nrow(test), cmp$pass, cmp$n_fail, cmp$max_relerr))
  if (!cmp$pass)
    cat(sprintf("   worst: t=%g state=%s ref=%.10g test=%.10g\n",
                cmp$worst_time, cmp$worst_state, cmp$worst_pair["ref"], cmp$worst_pair["test"]))
  cmp$pass
}

g24  <- gate1_one(24,  "tests/fixtures/det_24h.csv",  rtol = 1e-9)
g300 <- if (file.exists("tests/fixtures/det_300h.csv"))
          gate1_one(300, "tests/fixtures/det_300h.csv", rtol = 1e-6) else {
          cat("Gate 1 @ 300h: SKIPPED (det_300h.csv not present yet)\n"); NA }

cat(sprintf("\nGATE 1: %s\n",
    if (isTRUE(g24) && isTRUE(g300)) "PASS" else
    if (isTRUE(g24) && is.na(g300)) "24h PASS (300h pending)" else "FAIL"))
