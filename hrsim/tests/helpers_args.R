# helpers_args.R
#
# Part of the tests of the hrsim package, the simulation engine for Berryhill
# and Andersson (2026), "Unstable Heteroresistance as a Stepping Stone to
# Stable Resistance via Compensatory Mutations".
# Defines debug_args().

#' Build the argument list of debug_step1() and debug_step1_fired()
#'
#' Unpacks the parameter list, transition table and initial state in the same
#' way as hrsim::run_sim(), so the debug functions receive exactly what
#' simulate_hr() receives. The rate vector is built with the expression used
#' by run_sim() and by the R reference engine.
#'
#' @param p Parameter list from make_condition()$p.
#' @param transitions Data frame from build_transitions().
#' @param y0 Initial state: list(pop, r, A, B).
#' @return A named list with the 17 arguments of hrsim::debug_step1(). Append
#'   `seed` to it to call hrsim::debug_step1_fired().
debug_args <- function(p, transitions, y0) {
  base     <- c(gamma = p$gamma, mu = p$mu, omega = p$omega)
  rate_vec <- as.numeric(transitions$mult * base[transitions$type])
  list(
    init = as.numeric(y0$pop), r0 = y0$r, A0 = y0$A, B0 = y0$B,
    v = p$v, vm = p$vm, kH = p$kH, M = p$M, Pb = p$P,
    from = transitions$from, to = transitions$to, rate_vec = rate_vec,
    dt = p$DT, kmon = p$kmon, eR = p$e, d0 = p$d, KB = p$KB)
}
