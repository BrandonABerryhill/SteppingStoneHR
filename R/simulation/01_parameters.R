# 01_parameters.R
#
# The state space of the model, the founder population, and make_condition(),
# which returns the complete parameter set and initial state for one of the
# three conditions: Scenario 1, Scenario 2, or the control without
# compensation.
#
# Units. Time is in hours and the Euler time step is DT = 1e-5 h. Bacterial
# densities are in cells per mL, the resource r and the antibiotic A in µg per
# mL, and the β-lactamase B in units per mL (Supplementary Table 1). Each
# parameter is described where it is set in make_condition(). One that has a
# row in Supplementary Table 2 is given with its symbol, its dimensions as
# listed there, and the name of the row in quotation marks. The time step DT
# and the time between transfers have no row in that table; the Methods
# describe the fixed-step Euler integration and the serial transfer regimen.
#
# Part of the code for "Unstable Heteroresistance as a Stepping Stone to Stable
# Resistance via Compensatory Mutations" (Berryhill and Andersson).
#
# Defines: states, copy_number, base_y0() and make_condition().

# State space ----

# The twenty bacterial states: R0 to R9 carry 0 to 9 copies of the resistance
# gene, and R0C to R9C are the same with a compensatory mutation. Both engines,
# the transition table and the classifier index states by position in this
# vector, so the order must not change.
states      <- c(paste0("R", 0:9), paste0("R", 0:9, "C"))
# Copy number of each state, in the order of `states`.
copy_number <- c(0:9, 0:9)

# Initial state ----

#' Initial state of a run
#'
#' @param amax Numeric scalar, antibiotic concentration at time zero in µg per
#'   mL.
#' @return A list: `pop`, a numeric vector of length 20 named by `states`, the
#'   initial densities in cells per mL; `r`, the resource in µg per mL; `A`,
#'   the antibiotic in µg per mL; and `B`, the β-lactamase in units per mL.
base_y0 <- function(amax) {
  # Founder densities of R0 to R9. The compensated states start at zero.
  pop0 <- c(99, 34600, 2.01e7, 55300, 493, 211, 4.13, 2.44, 0.872, 0.434,
            rep(0, 10))
  names(pop0) <- states
  # The resource starts at 1000 µg per mL, the value that make_condition() also
  # sets as the maximum resource C. No β-lactamase is present at the start.
  list(pop = pop0, r = 1000, A = amax, B = 0)
}

# Condition builder ----

#' Parameters and initial state for one condition
#'
#' @param scenario Character scalar, one of "S1", "S2" or "omega0". "S1" is
#'   Scenario 1: every compensated state has an MIC of 25 µg per mL. "S2" is
#'   Scenario 2: the MIC of a compensated state is three times that of the
#'   non-compensated state with the same copy number. "omega0" is the control
#'   without compensation: the compensation rate is set to zero whatever
#'   `omega` is, so the compensated states are never populated.
#' @param cost_per_copy Numeric scalar, per cell per hour. Amount by which each
#'   copy of the resistance gene lowers the maximum growth rate. It has no row
#'   of its own in Supplementary Table 2; it enters the maximum growth rates
#'   `v`, whose rows are "Maximum Growth Rates" and "Maximum Growth Rates of
#'   Compensated States". The Methods (Sensitivity Analysis) call it the
#'   per-copy growth-rate cost.
#' @param omega Numeric scalar, per cell per hour. Compensation rate ω, used in
#'   both directions between a state and its compensated counterpart.
#' @param amax Numeric scalar, µg per mL. Antibiotic concentration at time zero
#'   and the dose AMAX added at each transfer.
#' @param mu Numeric scalar, per cell per hour. Base rate μ of stepwise
#'   recombination between copy numbers; the transition table multiplies it by
#'   a count that depends on the two copy numbers.
#' @param dil Numeric scalar between 0 and 1, dimensionless. Serial transfer
#'   rate: the fraction of every density and of the antibiotic carried over at
#'   a transfer.
#' @param comp_frac Numeric scalar between 0 and 1, dimensionless. Fraction of
#'   the per-copy cost that compensation recovers. It has no row of its own in
#'   Supplementary Table 2; it enters the maximum growth rates of the
#'   compensated states, whose row is "Maximum Growth Rates of Compensated
#'   States". The Supplementary Text (Heteroresistance Model) and the Methods
#'   (Sensitivity Analysis) describe it.
#' @return A list of two elements. `p` is the parameter list with the scalars
#'   `DT`, `transfer`, `C`, `e`, `kmon`, `AMAX`, `dil`, `d`, `KB`, `mu`,
#'   `gamma` and `omega`; the vectors `v`, `vm`, `kH`, `M` and `P`, each of
#'   length 20 in the order of `states`; the label `scenario`; and the number
#'   `cost_per_copy`. `y0` is the initial state returned by `base_y0(amax)`.
make_condition <- function(scenario = c("S1", "S2", "omega0"),
                           cost_per_copy = 0.0013,
                           omega = 1e-8,
                           amax  = 16,
                           mu        = 4e-4,
                           dil       = 0.01,
                           comp_frac = 0.5) {
  scenario <- match.arg(scenario)

  p <- list(
    # DT: Euler time step, in hours. transfer: time between serial transfers,
    # in hours.
    DT = 1e-5, transfer = 24,
    # C: maximum resource, µg per mL, the value the resource is reset to at
    # each transfer ("Maximum Resource").
    # e: conversion efficiency, µg per cell ("Conversion Efficiency").
    # kmon: Monod constant, µg per mL ("Monod Constant"), written k in
    # Supplementary Table 2 and in Supplementary Equation 2. In both engines
    # it is added to the resource r, which is also in µg per mL.
    C  = 1000, e = 5e-7, kmon = 1,
    # AMAX: antibiotic dose added at each transfer, µg per mL ("Antibiotic Dose
    # at Each Transfer").
    # dil: serial transfer rate, dimensionless ("Serial Transfer Rate").
    # d: antibiotic decay rate, per hour ("Antibiotic Decay Rate"). In both
    # engines d * A is a rate of change of the antibiotic A.
    # KB: β-lactamase association constant, mL per unit per hour
    # ("β-lactamase Association Constant"). In both engines KB * B * A is a
    # rate of change of A, where B is the β-lactamase concentration in units
    # per mL.
    AMAX = amax, dil = dil, d = 0, KB = 1,
    # mu and gamma: recombination rates μ and γ, per cell per hour
    # ("Recombination Transition Rates"). mu is the base rate of stepwise
    # changes in copy number; gamma is the rate of deletion to zero copies and
    # of doubling. 02_transitions.R assigns one of the two to each
    # recombination transition.
    mu = mu, gamma = 4e-7,
    # omega: compensation rate ω, per cell per hour ("Compensation Rate"). It
    # is the rate of the transition from each state to its compensated
    # counterpart and of the transition back. The control has no compensation.
    omega = if (scenario == "omega0") 0 else omega
  )

  # v: maximum growth rate of each state, per cell per hour ("Maximum Growth
  # Rates" and "Maximum Growth Rates of Compensated States"). Without
  # compensation it is 2 per hour at zero copies, and each copy lowers it by
  # cost_per_copy.
  v_R  <- 2 - cost_per_copy * (0:9)
  # Compensation recovers the fraction comp_frac of the per-copy cost.
  v_C  <- 2 - (1 - comp_frac) * cost_per_copy * (0:9)
  p$v  <- c(v_R, v_C)
  # vm: minimum growth rate v_min of Supplementary Equation 1, per cell per
  # hour ("Minimum Growth/Death Rate"). kH: Hill coefficient κ, dimensionless
  # ("Hill Coefficient"). Both are the same for every state.
  p$vm <- rep(-2, 20)
  p$kH <- rep(1, 20)

  # M: minimum inhibitory concentration MIC of each state, µg per mL ("MIC of
  # Each State" and "MIC of Each Compensated State", Scenario 1 and Scenario 2).
  # Without compensation it is 0.375 µg per mL at zero copies and rises by 1.5
  # µg per mL with each copy. A compensated state has 25 µg per mL in Scenario
  # 1 and three times the non-compensated value in Scenario 2. The control
  # uses the Scenario 1 values for the compensated states.
  M_R  <- 0.375 + 1.5 * (0:9)
  M_C  <- if (scenario == "S2") 3 * M_R else rep(25, 10)
  p$M  <- c(M_R, M_C)

  # P: β-lactamase production rate p of each state, units per µg per hour
  # ("β-lactamase Production Rate" and "β-lactamase Production Rate of
  # Compensated States"): 0.01 per copy, the same with and without
  # compensation.
  P_R  <- 0.01 * (0:9)
  p$P  <- c(P_R, P_R)

  # The scenario label names the folder a batch is written to. It and the value
  # of cost_per_copy are recorded in the manifest (05_run.R).
  p$scenario      <- scenario
  p$cost_per_copy <- cost_per_copy

  list(p = p, y0 = base_y0(amax))
}
