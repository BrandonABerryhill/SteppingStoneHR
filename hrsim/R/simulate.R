# simulate.R
#
# R wrapper around the compiled engine of the hrsim package, the simulation
# engine for Berryhill and Andersson (2026), "Unstable Heteroresistance as a
# Stepping Stone to Stable Resistance via Compensatory Mutations".
# Defines run_sim().

#' Run one replicate of the heteroresistance model with the compiled engine
#'
#' Takes the same named arguments as the R reference engine
#' \code{run_model()}; the order of \code{stochastic} and \code{n_out} and
#' the default of \code{n_out} differ between the two.
#' The function unpacks the parameter list, the transition table and the
#' initial state into the plain vectors and scalars that
#' \code{\link{simulate_hr}} takes, and returns its result as a data frame.
#'
#' Simulated time is in hours. The Euler step is \code{p$DT}.
#'
#' @param p Parameter list, as returned in \code{make_condition()$p} by the R
#'   reference implementation. Elements used: \code{DT} (Euler step, hours),
#'   \code{transfer} (hours between serial transfers), \code{C} (resource
#'   after each transfer), \code{e} (conversion efficiency), \code{kmon}
#'   (Monod constant), \code{AMAX} (antibiotic dose at each transfer),
#'   \code{dil} (dilution factor at each transfer), \code{d} (antibiotic decay
#'   rate), \code{KB} (beta-lactamase association constant), \code{mu},
#'   \code{gamma} and \code{omega} (transition rates, per cell per hour), and
#'   the length-20 vectors \code{v}, \code{vm}, \code{kH}, \code{M} and
#'   \code{P} (maximum growth rate, minimum growth rate, Hill coefficient, MIC
#'   and beta-lactamase production rate of each state).
#' @param transitions Data frame of transitions, as returned by
#'   \code{build_transitions()} in the R reference implementation, with
#'   columns \code{from} and \code{to} (state indices, 1 to 20), \code{type}
#'   (\code{"gamma"}, \code{"mu"} or \code{"omega"}) and \code{mult} (the
#'   number by which the base rate of that type is multiplied).
#' @param y0 Initial state: a list with \code{pop} (numeric vector of the 20
#'   state densities in cells per mL, in the order R0 to R9 then R0C to R9C)
#'   and the scalars \code{r}, \code{A} and \code{B} (resource, antibiotic and
#'   beta-lactamase).
#' @param t_end Simulated time in hours.
#' @param seed Seed of the engine's own random number generator: a
#'   non-negative whole number. It governs every random draw of the run, and
#'   a run's random stream depends only on its own seed. It has no effect when
#'   \code{stochastic = FALSE}.
#' @param stochastic Logical. If \code{FALSE} the Monte Carlo transition layer
#'   is disabled, so no transition fires and the run is a deterministic Euler
#'   integration.
#' @param n_out Approximate number of evenly spaced rows to record in addition
#'   to the transfer rows; must not exceed the number of steps,
#'   \code{round(t_end / p$DT)}. A value of 0 or less, the default, records
#'   only time 0 and the transfer rows. The R reference engine's default is
#'   400.
#' @return A data frame with 25 columns: \code{time} (hours), \code{R0} to
#'   \code{R9} and \code{R0C} to \code{R9C} (cells per mL), \code{r},
#'   \code{A}, \code{B} and \code{cycle_end}. \code{cycle_end} is 1 on the row
#'   recorded just before each serial transfer and 0 on every other row.
#' @export
run_sim <- function(p, transitions, y0, t_end, seed = 1L,
                    stochastic = TRUE, n_out = -1L) {
  # A grid finer than the step count would record every step; refuse it here,
  # before the compiled code is called.
  n_steps <- round(t_end / p$DT)
  if (n_out > 0 && n_out > n_steps)
    stop(sprintf("n_out (%d) must be <= n_steps (%.0f), or <= 0 for transfer-rows-only; R uses 400.",
                 as.integer(n_out), n_steps))
  base     <- c(gamma = p$gamma, mu = p$mu, omega = p$omega)
  rate_vec <- as.numeric(transitions$mult * base[transitions$type])
  out <- simulate_hr(
    init = as.numeric(y0$pop), r0 = y0$r, A0 = y0$A, B0 = y0$B,
    v = p$v, vm = p$vm, kH = p$kH, M = p$M, Pb = p$P,
    from = transitions$from, to = transitions$to, rate_vec = rate_vec,
    dt = p$DT, horizon_h = t_end, transfer_interval_h = p$transfer,
    kmon = p$kmon, eR = p$e, d0 = p$d, KB = p$KB,
    dil = p$dil, AMAX = p$AMAX, Cres = p$C,
    n_out = as.integer(n_out),
    seed = as.numeric(seed), stochastic = stochastic)
  as.data.frame(out)
}
