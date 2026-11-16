# run_all_gates.R
#
# Driver for the verification tests of the hrsim package, the simulation
# engine for Berryhill and Andersson (2026), "Unstable Heteroresistance as a
# Stepping Stone to Stable Resistance via Compensatory Mutations". The tests
# are called gates. This script sources the four check scripts in turn, each
# in its own environment. Each prints what it measured and ends with a line
# that begins with "GATE".
#
# Run it with hrsim installed and the working directory at the package root:
#   Rscript -e 'source("tests/run_all_gates.R")'
#
# What is run
#   gate1b_check.R  the rates and transition probabilities of the first step
#                   against the R reference, and firing frequencies over
#                   seeds 1 to N
#   gate1_check.R   trajectories with the Monte Carlo layer disabled, over
#                   24 h and 300 h, against the R reference
#   gate3_check.R   invariants of recorded trajectories
#   gate2_check.R   outcome frequencies in a stored sample of stochastic
#                   runs; this script runs no simulation
# The tolerance of each test, the reason for it and the seeds used are given
# in the header of each script.
#
# Options
#   hrsim.ref         folder that holds the R reference implementation
#                     (default "tests/ref")
#   hrsim.firerate_n  number of seeds in the firing-frequency test of
#                     gate1b_check.R (default 200000)

if (is.null(getOption("hrsim.ref"))) options(hrsim.ref = "tests/ref")

banner <- function(s) cat(sprintf("\n==================== %s ====================\n", s))
run_check <- function(f) { banner(basename(f)); source(f, local = new.env(parent = globalenv())) }

run_check("tests/gate1b_check.R")
run_check("tests/gate1_check.R")
run_check("tests/gate3_check.R")
run_check("tests/gate2_check.R")

banner("SUMMARY")
cat("Gate 1b / Gate 1 / Gate 3: see PASS lines above (run live).\n")
cat("Gate 2: persisted sample in tests/fixtures/gate2_results.csv; extend toward the\n")
cat("        >=200-seed target with tests/gate2_run.R (multi-core host recommended).\n")
