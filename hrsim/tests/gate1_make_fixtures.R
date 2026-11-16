# gate1_make_fixtures.R
#
# Part of the tests of the hrsim package, the simulation engine for Berryhill
# and Andersson (2026), "Unstable Heteroresistance as a Stepping Stone to
# Stable Resistance via Compensatory Mutations".
#
# Defines the two functions that write the deterministic reference
# trajectories read by gate1_check.R. Each runs the R reference engine,
# run_model(), for the Scenario 1 condition with stochastic = FALSE and
# n_out = 400, and writes the result as CSV under tests/fixtures/. Sourcing
# this file loads the R reference implementation through load_canonical.R,
# creates the folder tests/fixtures/ if it does not exist and defines the two
# functions; it runs no simulation and writes no fixture. Call one of the
# functions to rewrite its fixture.
#
#   make_24h()          writes det_24h.csv  (24 h, 2.4e6 Euler steps)
#   make_300h_direct()  writes det_300h.csv (300 h, 3e7 Euler steps in
#                       interpreted R, which takes correspondingly longer)
#
# No random number is drawn, so the output does not depend on a seed. Run
# from the package root.

ref_dir <- getOption("hrsim.ref", "tests/ref")
source(file.path(ref_dir, "load_canonical.R"))
dir.create("tests/fixtures", recursive = TRUE, showWarnings = FALSE)

make_24h <- function() {
  d <- run_model(p, transitions, y0, t_end = 24, n_out = 400, stochastic = FALSE)
  write.csv(d, "tests/fixtures/det_24h.csv", row.names = FALSE)
  cat(sprintf("det_24h.csv: %d rows, %d transfer rows\n", nrow(d), sum(d$cycle_end == 1)))
}

make_300h_direct <- function() {
  d <- run_model(p, transitions, y0, t_end = 300, n_out = 400, stochastic = FALSE)
  write.csv(d, "tests/fixtures/det_300h.csv", row.names = FALSE)
  cat(sprintf("det_300h.csv: %d rows, %d transfer rows\n", nrow(d), sum(d$cycle_end == 1)))
}
