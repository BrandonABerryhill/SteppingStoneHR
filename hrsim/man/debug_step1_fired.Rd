\name{debug_step1_fired}
\alias{debug_step1_fired}
\title{Transitions Fired in a Single Stochastic Step}
\description{
Computes the transition probabilities for the state supplied, as
\code{\link{debug_step1}} does, and applies the engine's Monte Carlo rule to
them once, with a generator seeded from \code{seed}. The rule and the order
of the random draws are those of \code{\link{simulate_hr}}, so for the same
inputs and seed the result is the set of transitions \code{simulate_hr} fires
in its first step. The function exists for the package's tests.
}
\usage{
debug_step1_fired(init, r0, A0, B0, v, vm, kH, M, Pb, from, to, rate_vec,
                  dt, kmon, eR, d0, KB, seed)
}
\arguments{
  \item{init}{Numeric vector of length 20: initial densities of the states in
    cells per mL, in the order R0 to R9 then R0C to R9C.}
  \item{r0, A0, B0}{Initial resource (ug per mL), antibiotic (ug per mL) and
    beta-lactamase (units per mL).}
  \item{v}{Numeric vector of length 20: maximum growth rate of each state,
    per cell per hour.}
  \item{vm}{Numeric vector of length 20: minimum growth rate of each state,
    per cell per hour.}
  \item{kH}{Numeric vector of length 20: Hill coefficient of each state,
    dimensionless.}
  \item{M}{Numeric vector of length 20: MIC of each state, ug per mL.}
  \item{Pb}{Numeric vector of length 20: beta-lactamase production rate of
    each state, units per ug per hour.}
  \item{from, to}{Integer vectors, one element per transition: the source and
    target state of each transition as 1-based indices (1 to 20) into the
    state order.}
  \item{rate_vec}{Numeric vector, one element per transition: its rate
    constant, per cell per hour. The number of transitions is taken from the
    length of this vector.}
  \item{dt}{Euler step, hours.}
  \item{kmon}{Monod constant of the resource limitation term.}
  \item{eR}{Conversion efficiency.}
  \item{d0}{Antibiotic decay rate.}
  \item{KB}{Beta-lactamase association constant.}
  \item{seed}{Seed of the random stream: a non-negative whole number, passed
    as a double.}
}
\details{
A transition whose probability is 0 or less never fires and one whose
probability is 1 or more always fires; neither consumes a random number.
Every other transition consumes one uniform draw, in the order of the
transition vectors, and fires when the draw is below its probability.

The lengths of the vectors and the range of the indices are not checked.
\code{to} does not enter the calculation.
}
\value{
A logical vector with one element per transition, in the order of
\code{from}: \code{TRUE} where the transition fired.
}
\seealso{\code{\link{debug_step1}}, \code{\link{simulate_hr}}}
\references{
Berryhill, B. A. and Andersson, D. I. (2026). Unstable Heteroresistance as a
Stepping Stone to Stable Resistance via Compensatory Mutations. bioRxiv.
doi:10.64898/2026.10.08.756632
}
