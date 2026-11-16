\name{simulate_hr}
\alias{simulate_hr}
\title{Run One Replicate with the Compiled Engine}
\description{
Integrates the model from time 0 to \code{horizon_h} and returns the recorded
rows. Every input is a separate vector or scalar; \code{\link{run_sim}}
builds this argument list from a parameter list, a transition table and an
initial state, and is the more convenient entry point.
}
\usage{
simulate_hr(init, r0, A0, B0, v, vm, kH, M, Pb, from, to, rate_vec, dt,
            horizon_h, transfer_interval_h, kmon, eR, d0, KB, dil, AMAX,
            Cres, n_out, seed, stochastic = TRUE)
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
  \item{horizon_h}{Simulated time, hours. The number of steps is
    \code{horizon_h / dt}, rounded to the nearest whole number.}
  \item{transfer_interval_h}{Time between serial transfers, hours.}
  \item{kmon}{Monod constant of the resource limitation term.}
  \item{eR}{Conversion efficiency.}
  \item{d0}{Antibiotic decay rate.}
  \item{KB}{Beta-lactamase association constant.}
  \item{dil}{Factor by which the populations and the antibiotic are
    multiplied at each transfer, dimensionless.}
  \item{AMAX}{Antibiotic added at each transfer, ug per mL.}
  \item{Cres}{Resource concentration set at each transfer, ug per mL.}
  \item{n_out}{Integer. Approximate number of evenly spaced grid rows to
    record in addition to the transfer rows; a grid row is recorded every
    \code{max(1, floor(n_steps / n_out))} steps. It must not exceed the
    number of steps. A value of 0 or less records only time 0 and the
    transfer rows.}
  \item{seed}{Seed of the run's random stream: a non-negative whole number,
    passed as a double. The stream depends on this value alone.}
  \item{stochastic}{Logical. If \code{FALSE} the Monte Carlo transition layer
    is disabled: no transition fires and no random number is drawn.}
}
\details{
The lengths of the vectors and the range of \code{from} and \code{to} are not
checked, and \code{run_sim} passes them on from its own arguments without
checking them either. The caller must ensure that the six per-state vectors
have length 20 and that \code{from} and \code{to} have the length of
\code{rate_vec} and values from 1 to 20.

The function stops with an error if \code{n_out} exceeds the number of steps
or if the number of rows it allocates for the output exceeds 1e8. That number
is an upper bound on the number of rows returned: the grid rows, plus
\code{ceiling(horizon_h / transfer_interval_h)}, plus 5.

The clock that triggers a transfer is advanced by adding \code{dt} at every
step. With \code{dt = 1e-5} and an interval of 24 the first transfer falls
one step after 24 h, at a recorded time of 24.00001, so a run of exactly
24 h contains no transfer row.
}
\value{
A numeric matrix with 25 named columns, in this order: \code{time} (hours),
\code{R0} to \code{R9} and \code{R0C} to \code{R9C} (cells per mL), \code{r},
\code{A}, \code{B} and \code{cycle_end}. Rows are in time order. The first
row is the initial state at time 0. Each transfer adds one row holding the
state just before dilution, with \code{cycle_end = 1}. Grid rows have
\code{cycle_end = 0}; a grid row that falls on a transfer step follows the
transfer row and holds the state after the transfer.

The matrix has two attributes: \code{"final_state"}, a numeric vector of
length 23 with the 20 densities, \code{r}, \code{A} and \code{B} after the
last step, and \code{"n_steps"}, the number of steps taken.
}
\seealso{\code{\link{run_sim}}}
\references{
Berryhill, B. A. and Andersson, D. I. (2026). Unstable Heteroresistance as a
Stepping Stone to Stable Resistance via Compensatory Mutations. bioRxiv.
doi:10.64898/2026.10.08.756632
}
