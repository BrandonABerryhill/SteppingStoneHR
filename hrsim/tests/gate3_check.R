# gate3_check.R
#
# Part of the tests of the hrsim package, the simulation engine for Berryhill
# and Andersson (2026), "Unstable Heteroresistance as a Stepping Stone to
# Stable Resistance via Compensatory Mutations".
#
# What it verifies
#   Properties that must hold whatever the trajectory. I1 to I3 are tested on
#   every recorded row, I4 on the transitions fired in a first step:
#     I1  every state density is >= 0 (the engine floors them each step)
#     I2  the antibiotic is >= 0 (also floored each step)
#     I3  every recorded value is finite
#     I4  the transitions fired in a step add as many cells as they remove,
#         so the net change they cause sums to zero over the 20 states
#   The resource r and beta-lactamase B are not floored by the engine. Their
#   minima over the recorded rows are printed at the end of each trajectory
#   line for information and are not tested. The output line labelled I5 is a
#   fixed note that negative values of r or B are not errors; it reports no
#   value.
#
# Tolerance
#   None: I1 to I4 are exact conditions.
#
# Seeds
#   I1 to I3 are checked on one 24 h run with the Monte Carlo layer disabled
#   and on five stochastic 48 h runs with seeds 1 to 5, all for the Scenario 1
#   condition on the n_out = 400 grid. I4 is checked on the first step for
#   seeds 1 to 2000, through hrsim::debug_step1_fired(). Each seed's random
#   stream depends only on that seed.
#
# Run from the package root with hrsim installed.

suppressMessages(library(hrsim))
ref_dir <- getOption("hrsim.ref", "tests/ref")
source(file.path(ref_dir, "load_canonical.R"))
source("tests/helpers_args.R")

pop_cols <- c(paste0("R", 0:9), paste0("R", 0:9, "C"))

# Tests I1 to I3 on one trajectory and prints one line. Returns the three
# results as a logical vector.
check_traj <- function(d, tag) {
  P <- as.matrix(d[, pop_cols, drop = FALSE])
  i1 <- all(P >= 0)
  i2 <- all(d$A >= 0)
  i3 <- all(is.finite(as.matrix(d)))
  cat(sprintf("  [%s] rows=%d  I1 pop>=0:%s  I2 A>=0:%s  I3 finite:%s  | min(r)=%.4g min(B)=%.4g\n",
              tag, nrow(d), i1, i2, i3, min(d$r), min(d$B)))
  c(i1, i2, i3)
}

cat("== Gate 3: trajectory invariants ==\n")
res <- c()
res <- c(res, check_traj(run_sim(p, transitions, y0, 24, stochastic = FALSE, n_out = 400L), "det 24h"))
for (s in 1:5)
  res <- c(res, check_traj(run_sim(p, transitions, y0, 48, seed = s, stochastic = TRUE, n_out = 400L),
                           sprintf("stoch 48h s=%d", s)))

# I4: rebuild the net change per state from the set of fired transitions.
args <- debug_args(p, transitions, y0)
from <- transitions$from; to <- transitions$to
bad_net <- 0L
for (s in 1:2000) {
  fired <- as.logical(do.call(hrsim::debug_step1_fired, c(args, list(seed = s))))
  if (any(fired)) {
    net <- tabulate(to[fired], 20) - tabulate(from[fired], 20)
    if (sum(net) != 0L) bad_net <- bad_net + 1L
  }
}
i4 <- (bad_net == 0L)
cat(sprintf("  [I4] sum(net)==0 over 2000 seeded steps: %s (violations=%d)\n", i4, bad_net))
cat(sprintf("  [I5] r/B unclamped: by construction (never floored); negatives above are not errors\n"))

cat(sprintf("\nGATE 3: %s\n", if (all(res) && i4) "PASS" else "FAIL"))
