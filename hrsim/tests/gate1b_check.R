# gate1b_check.R
#
# Part of the tests of the hrsim package, the simulation engine for Berryhill
# and Andersson (2026), "Unstable Heteroresistance as a Stepping Stone to
# Stable Resistance via Compensatory Mutations".
#
# What it verifies
#   The part of a step that gate1_check.R cannot see, because that test runs
#   with the Monte Carlo layer disabled: the transition probabilities, and
#   the rule that turns them into firings. Both parts use the Scenario 1
#   condition at its initial state. The two parts are labelled C.3 and C.4 in
#   the printed output.
#
#   Part C.3 compares the quantities computed by hrsim::debug_step1() (psi,
#   growth, dr, dB, dA, probs) element by element with the values the R
#   reference expressions give, stored in tests/fixtures/step1_internals.rds
#   by gate1b_dump_step1.R. It also checks that the transition indices come
#   back 1-based and equal to R's, and that the number of probabilities that
#   are exactly zero equals the number of transitions whose source state is
#   empty.
#
#   Part C.4 calls hrsim::debug_step1_fired() once per seed and compares the
#   fraction of seeds in which each transition fired with its probability.
#
# Tolerance
#   C.3 passes when the largest relative error over psi, growth, dr, dB and
#   probs is at most 1e-9, the absolute error of dA is at most 1e-12 (its
#   reference value is zero, so it is not compared by relative error), and
#   the two checks on the indices and on the zero probabilities hold. A
#   single step has no accumulated error, so agreement is expected at
#   rounding level; the output also states whether the tighter 1e-12 was met.
#   C.4 passes when no transition with zero probability fired and, for every
#   other transition, the observed frequency lies within 5 binomial standard
#   errors, sqrt(p * (1 - p) / N), of its probability p.
#
# Seeds
#   C.3 draws no random number. C.4 uses seeds 1 to N, where N is the option
#   hrsim.firerate_n (default 200000). Each seed's random stream depends only
#   on that seed. The seeds are fixed, so the result is the same on every run.
#
# Run from the package root with hrsim installed.

suppressMessages(library(hrsim))
ref_dir <- getOption("hrsim.ref", "tests/ref")
source(file.path(ref_dir, "load_canonical.R"))
source("tests/helpers_args.R")

ref  <- readRDS("tests/fixtures/step1_internals.rds")
args <- debug_args(p, transitions, y0)
cpp  <- do.call(hrsim::debug_step1, args)

# C.3: element-by-element agreement ----

rel <- function(a, b) { d <- abs(a - b); d / pmax(abs(b), 1e-300) }
fields <- c("psi","growth","dr","dB","dA","probs")
worst <- sapply(fields, function(f) max(rel(cpp[[f]], ref[[f]])))
# The reference dA is zero at the initial state, so a relative error is not
# defined for it; the absolute error takes its place.
abs_dA <- abs(cpp$dA - ref$dA)
worst["dA"] <- abs_dA

cat("== Gate 1b C.3: elementwise kernel vs R dump ==\n")
for (f in fields) cat(sprintf("  %-7s max %s = %.3e\n", f,
    if (f=="dA") "abserr" else "relerr", worst[f]))

# The compiled code works with 0-based indices internally; what it returns
# must be the 1-based indices it was given.
idx_ok <- identical(as.integer(cpp$from_1based), as.integer(ref$from)) &&
          identical(as.integer(cpp$to_1based),   as.integer(ref$to))
cat(sprintf("  index integrity (from/to 1-based match): %s\n", idx_ok))

# In the condition tested here, Scenario 1 at its initial state, every rate
# constant and the resource term are positive, so a probability is exactly
# zero only where the source state is empty.
n_zero_probs <- sum(cpp$probs == 0)
n_empty_src  <- sum(ref$pop_from == 0)
cat(sprintf("  probs==0 count = %d ; empty-source (pop[from]==0) = %d ; match = %s\n",
            n_zero_probs, n_empty_src, n_zero_probs == n_empty_src))

tol_primary <- 1e-12; tol_fallback <- 1e-9
mx <- max(worst[c("psi","growth","dr","dB","probs")])
c3_pass <- mx <= tol_fallback && abs_dA <= 1e-12 && idx_ok && (n_zero_probs == n_empty_src)
c3_tier <- if (mx <= tol_primary) "<=1e-12" else if (mx <= tol_fallback) "<=1e-9 (fallback)" else ">1e-9 FAIL"
cat(sprintf("  -> C.3 max relerr %.3e (%s)\n", mx, c3_tier))

# C.4: firing frequency over many seeds ----

N <- as.integer(getOption("hrsim.firerate_n", 200000L))
NTR <- length(ref$probs)
counts <- integer(NTR)
for (s in 1:N) counts <- counts + as.integer(do.call(hrsim::debug_step1_fired,
                                                      c(args, list(seed = s))))
emp  <- counts / N
pk   <- ref$probs
se   <- sqrt(pmax(pk * (1 - pk), 0) / N)

# A transition whose probability is zero must never fire.
empty <- which(pk <= 0)
empty_fires <- sum(counts[empty])
cat(sprintf("\n== Gate 1b C.4: fire-rate over N=%d seeds ==\n", N))
cat(sprintf("  empty-source transitions: %d ; total fires among them: %d (must be 0)\n",
            length(empty), empty_fires))

# z is the distance between observed frequency and probability in standard
# errors. The five most probable transitions are printed; all are tested.
active <- which(pk > 0)
z <- ifelse(se[active] > 0, abs(emp[active] - pk[active]) / se[active], 0)
hi <- active[order(pk[active], decreasing = TRUE)][1:min(5, length(active))]
cat("  top transitions (k: p_theory  p_emp  z):\n")
for (k in hi) cat(sprintf("    %3d: %.6f  %.6f  z=%.2f\n", k, pk[k], emp[k], z[match(k, active)]))
firerate_pass <- (empty_fires == 0) && all(z < 5)
cat(sprintf("  max z over active transitions: %.2f (threshold 5)\n", max(z)))

cat(sprintf("\nGATE 1b: %s\n", if (c3_pass && firerate_pass) "PASS" else "FAIL"))
