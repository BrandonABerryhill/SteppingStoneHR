# load_canonical.R
#
# Loads the R reference implementation the compiled engine was verified
# against (the three other files in this folder) and builds the objects the
# tests use: the Scenario 1 condition with its default arguments, its initial
# state and the transition table. Part of the tests of the hrsim package, the
# simulation engine for Berryhill and Andersson (2026), "Unstable
# Heteroresistance as a Stepping Stone to Stable Resistance via Compensatory
# Mutations".
#
# Defines
#   p            parameter list of the Scenario 1 condition
#   y0           its initial state
#   transitions  the transition table
# and, through the sourced files, states, copy_number, base_y0(),
# make_condition(), build_transitions() and run_model().
#
# The option hrsim.ref names this folder; the default assumes the working
# directory is the package root.

.ref <- getOption("hrsim.ref", "tests/ref")
source(file.path(.ref, "01_parameters.R"))
source(file.path(.ref, "02_transitions.R"))
source(file.path(.ref, "03_engine_R.R"))

.cond       <- make_condition("S1")
p           <- .cond$p
y0          <- .cond$y0
transitions <- build_transitions()

# Guards against a reference that was changed or loaded incompletely: the
# model has 150 transitions and 20 states in the fixed order.
stopifnot(nrow(transitions) == 150L,
          length(y0$pop) == 20L,
          identical(names(y0$pop), states))
