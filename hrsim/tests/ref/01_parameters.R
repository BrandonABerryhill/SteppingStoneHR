# 01_parameters.R
#
# State space, parameter values and initial state of the heteroresistance
# model, as used by the R reference implementation the compiled engine was
# verified against. Belongs to Berryhill and Andersson (2026), "Unstable
# Heteroresistance as a Stepping Stone to Stable Resistance via Compensatory
# Mutations". Supplementary Table 2 is the parameter table of that paper.
#
# Defines
#   states, copy_number   the 20 states and the gene copy number of each
#   base_y0()             initial state
#   make_condition()      parameter list and initial state of one condition
#
# Simulated time is in hours throughout, and the Euler step is DT = 1e-5 h.

# State space ----

# The order R0..R9, R0C..R9C is fixed: every per-state vector, the transition
# indices and the engine's output columns rely on it. R0..R9 carry 0 to 9
# copies of the beta-lactamase gene; R0C..R9C are the same with a compensatory
# mutation.
states      <- c(paste0("R", 0:9), paste0("R", 0:9, "C"))
copy_number <- c(0:9, 0:9)

# Initial state ----

#' Initial state of a run
#'
#' @param amax Antibiotic concentration at time 0, ug per mL.
#' @return A list: `pop`, the named densities of the 20 states in cells per mL
#'   (the founder population in R0..R9, zero in every compensated state);
#'   `r`, the resource, 1000 ug per mL; `A`, the antibiotic, equal to `amax`;
#'   and `B`, beta-lactamase, 0.
base_y0 <- function(amax) {
  pop0 <- c(99, 34600, 2.01e7, 55300, 493, 211, 4.13, 2.44, 0.872, 0.434,
            rep(0, 10))
  names(pop0) <- states
  list(pop = pop0, r = 1000, A = amax, B = 0)
}

# Condition builder ----

#' Parameter list and initial state of one condition
#'
#' @param scenario One of "S1", "S2" or "omega0".
#'   "S1": Scenario 1, every compensated state has an MIC of 25 ug per mL.
#'   "S2": Scenario 2, the MIC of each compensated state is 3 times that of
#'   the non-compensated state with the same copy number.
#'   "omega0": control without compensation; omega is set to 0 whatever the
#'   `omega` argument, so the compensated states are never populated. Its
#'   compensated MICs are those of Scenario 1.
#' @param cost_per_copy Reduction of the maximum growth rate per gene copy,
#'   per cell per hour. Default 0.0013.
#' @param omega Compensation rate, per cell per hour (Supplementary Table 2,
#'   Compensation Rate). Default 1e-8.
#' @param amax Antibiotic dose added at each transfer, and the concentration
#'   at time 0, ug per mL (Supplementary Table 2, Antibiotic Dose at Each
#'   Transfer). Default 16.
#' @return A list with `p`, the parameter list read by run_model() and
#'   hrsim::run_sim(), and `y0`, the initial state from base_y0().
make_condition <- function(scenario = c("S1", "S2", "omega0"),
                           cost_per_copy = 0.0013,
                           omega = 1e-8,
                           amax  = 16) {
  scenario <- match.arg(scenario)

  # Scalar parameters. The name in brackets is the row of Supplementary
  # Table 2; the symbol is the paper's where it differs from the name here.
  # Where a row is named, the unit given is the one listed in that row.
  #   DT        Euler step, hours
  #   transfer  time between serial transfers, hours
  #   C         resource after each transfer, ug per mL (Maximum Resource)
  #   e         conversion efficiency, ug per cell (Conversion Efficiency)
  #   kmon      Monod constant k, ug per mL (Monod Constant); the R
  #             reference engine and the compiled engine both add it to the
  #             resource r, which is also in ug per mL
  #   AMAX      antibiotic dose, ug per mL (Antibiotic Dose at Each Transfer)
  #   dil       dilution factor at transfer, dimensionless (Serial Transfer
  #             Rate)
  #   d         antibiotic decay rate, per hour (Antibiotic Decay Rate); in
  #             both engines d * A is a term of the rate of change of the
  #             antibiotic A
  #   KB        beta-lactamase association constant, mL per unit per hour
  #             (beta-lactamase Association Constant); in both engines
  #             KB * B * A is a term of the rate of change of A, where B is
  #             the beta-lactamase concentration in units per mL
  #   mu        recombination rate of a change from i copies to another copy
  #             number from 1 to min(2 * i, 9), doubling excepted, per cell
  #             per hour; multiplied by the per-transition count that
  #             build_transitions() gives (Recombination Transition Rates)
  #   gamma     recombination rate of deletion to zero copies and of
  #             doubling, per cell per hour (Recombination Transition Rates)
  #   omega     compensation rate, per cell per hour (Compensation Rate)
  p <- list(
    DT = 1e-5, transfer = 24,
    C  = 1000, e = 5e-7, kmon = 1,
    AMAX = amax, dil = 0.01, d = 0, KB = 1,
    mu = 4e-4, gamma = 4e-7,
    omega = if (scenario == "omega0") 0 else omega
  )

  # Maximum growth rates v, per cell per hour (Supplementary Table 2, Maximum
  # Growth Rates for the non-compensated states and Maximum Growth Rates of
  # Compensated States for the compensated ones): 2 without the gene, less
  # cost_per_copy for each copy. A compensated state recovers half of that
  # cost, which puts its rate midway between the non-compensated rate and 2.
  # vm is the minimum growth rate vMIN, per cell per hour (Minimum
  # Growth/Death Rate), and kH the Hill coefficient kappa, dimensionless
  # (Hill Coefficient); both are the same for every state.
  v_R  <- 2 - cost_per_copy * (0:9)
  v_C  <- (v_R + 2) / 2
  p$v  <- c(v_R, v_C)
  p$vm <- rep(-2, 20)
  p$kH <- rep(1, 20)

  # MIC of each state, M, in ug per mL: 0.375 without the gene plus 1.5 per
  # copy for the non-compensated states (Supplementary Table 2, MIC of Each
  # State). The row for the compensated states is MIC of Each Compensated
  # State, Scenario 2 when scenario is "S2" and MIC of Each Compensated State,
  # Scenario 1 otherwise.
  M_R  <- 0.375 + 1.5 * (0:9)
  M_C  <- if (scenario == "S2") 3 * M_R else rep(25, 10)
  p$M  <- c(M_R, M_C)

  # Beta-lactamase production rate of each state, p in the paper and P here,
  # in units per ug per hour (Supplementary Table 2, beta-lactamase
  # Production Rate for the non-compensated states and beta-lactamase
  # Production Rate of Compensated States for the compensated ones): 0.01 per
  # gene copy, the same with and without compensation.
  P_R  <- 0.01 * (0:9)
  p$P  <- c(P_R, P_R)

  # Labels carried with the parameters; the engines do not read them.
  p$scenario      <- scenario
  p$cost_per_copy <- cost_per_copy

  list(p = p, y0 = base_y0(amax))
}
