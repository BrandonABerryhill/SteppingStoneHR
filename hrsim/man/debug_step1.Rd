\name{debug_step1}
\alias{debug_step1}
\title{Rates and Transition Probabilities of a Single Step}
\description{
Evaluates, for the state supplied, the quantities the engine computes at the
start of a step, using the same compiled code as \code{\link{simulate_hr}}.
The state is not advanced and no random number is drawn. The function exists
so that the package's tests can compare these quantities with the R
reference.
}
\usage{
debug_step1(init, r0, A0, B0, v, vm, kH, M, Pb, from, to, rate_vec, dt,
            kmon, eR, d0, KB)
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
}
\details{
The arguments have the same meaning as in \code{\link{simulate_hr}}, and as
there the lengths of the vectors and the range of the indices are not
checked. \code{to} does not enter the calculation and is only returned.
}
\value{
A named list:
\describe{
  \item{\code{psi}}{the resource limitation term, \code{r0 / (r0 + kmon)}.}
  \item{\code{growth}}{numeric vector of length 20: the Euler increment of
    each state over one step, cells per mL.}
  \item{\code{dr}, \code{dB}, \code{dA}}{the time derivatives of the
    resource, beta-lactamase and the antibiotic.}
  \item{\code{probs}}{numeric vector, one element per transition: the
    probability that the transition fires in the step,
    \code{rate_vec * dt * init[from] * psi}.}
  \item{\code{from_1based}, \code{to_1based}}{the \code{from} and \code{to}
    vectors as supplied.}
}
}
\seealso{\code{\link{debug_step1_fired}}, \code{\link{simulate_hr}}}
\references{
Berryhill, B. A. and Andersson, D. I. (2026). Unstable Heteroresistance as a
Stepping Stone to Stable Resistance via Compensatory Mutations. bioRxiv.
doi:10.64898/2026.10.08.756632
}
