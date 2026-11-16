\name{run_sim}
\alias{run_sim}
\title{Run One Replicate of the Heteroresistance Model with the Compiled Engine}
\description{
Takes the same named arguments as the R reference engine \code{run_model()};
the order of \code{stochastic} and \code{n_out} and the default of
\code{n_out} differ between the two. The function unpacks the
parameter list, the transition table and the initial state into the plain
vectors and scalars that \code{\link{simulate_hr}} takes, and returns its
result as a data frame.
}
\usage{
run_sim(p, transitions, y0, t_end, seed = 1L, stochastic = TRUE,
        n_out = -1L)
}
\arguments{
  \item{p}{Parameter list, as returned in \code{make_condition()$p} by the R
    reference implementation. Elements used: \code{DT} (Euler step, hours),
    \code{transfer} (hours between serial transfers), \code{C} (resource
    after each transfer), \code{e} (conversion efficiency), \code{kmon}
    (Monod constant), \code{AMAX} (antibiotic dose at each transfer),
    \code{dil} (dilution factor at each transfer), \code{d} (antibiotic decay
    rate), \code{KB} (beta-lactamase association constant), \code{mu},
    \code{gamma} and \code{omega} (transition rates, per cell per hour), and
    the length-20 vectors \code{v}, \code{vm}, \code{kH}, \code{M} and
    \code{P} (maximum growth rate, minimum growth rate, Hill coefficient, MIC
    and beta-lactamase production rate of each state).}
  \item{transitions}{Data frame of transitions, as returned by
    \code{build_transitions()} in the R reference implementation, with
    columns \code{from} and \code{to} (state indices, 1 to 20), \code{type}
    (\code{"gamma"}, \code{"mu"} or \code{"omega"}) and \code{mult} (the
    number by which the base rate of that type is multiplied).}
  \item{y0}{Initial state: a list with \code{pop} (numeric vector of the 20
    state densities in cells per mL, in the order R0 to R9 then R0C to R9C)
    and the scalars \code{r}, \code{A} and \code{B} (resource, antibiotic and
    beta-lactamase).}
  \item{t_end}{Simulated time in hours.}
  \item{seed}{Seed of the engine's own random number generator: a
    non-negative whole number. It governs every random draw of the run, and
    a run's random stream depends only on its own seed. It has no effect when
    \code{stochastic = FALSE}.}
  \item{stochastic}{Logical. If \code{FALSE} the Monte Carlo transition layer
    is disabled, so no transition fires and the run is a deterministic Euler
    integration.}
  \item{n_out}{Approximate number of evenly spaced rows to record in addition
    to the transfer rows; must not exceed the number of steps,
    \code{round(t_end / p$DT)}. A value of 0 or less, the default, records
    only time 0 and the transfer rows. The R reference engine's default is
    400.}
}
\details{
Simulated time is in hours. The Euler step is \code{p$DT}.
}
\value{
A data frame with 25 columns: \code{time} (hours), \code{R0} to \code{R9} and
\code{R0C} to \code{R9C} (cells per mL), \code{r}, \code{A}, \code{B} and
\code{cycle_end}. \code{cycle_end} is 1 on the row recorded just before each
serial transfer and 0 on every other row.
}
\seealso{\code{\link{simulate_hr}}}
\references{
Berryhill, B. A. and Andersson, D. I. (2026). Unstable Heteroresistance as a
Stepping Stone to Stable Resistance via Compensatory Mutations. bioRxiv.
doi:10.64898/2026.10.08.756632
}
