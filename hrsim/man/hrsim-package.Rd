\name{hrsim-package}
\alias{hrsim-package}
\alias{hrsim}
\docType{package}
\title{Simulation Engine for the Heteroresistance Model with Compensatory Mutations}
\description{
Compiled simulation engine for the model in Berryhill and Andersson (2026),
"Unstable Heteroresistance as a Stepping Stone to Stable Resistance via
Compensatory Mutations".
}
\details{
The model has twenty bacterial states. R0 to R9 carry 0 to 9 copies of a
beta-lactamase gene, and R0C to R9C are the same with a compensatory
mutation. The states grow on a limiting resource in the presence of an
antibiotic that the beta-lactamase they produce degrades, and the culture is
serially transferred at a fixed interval.

The engine advances the model by fixed-step explicit Euler integration.
At every step each transition between states (recombination between copy
numbers, or compensation) can fire, with a probability computed from the
state at the start of the step; a transition that fires moves one cell per mL
from its source state to its target state. At each transfer the populations
and the antibiotic are diluted, a fresh antibiotic dose is added, the
resource is restored and beta-lactamase is set to zero. Simulated time is in
hours.

Random numbers come from a generator inside the compiled code
(xoshiro256++, seeded through SplitMix64 from the seed of the run). No draw
is taken from R's random number generator, and \code{set.seed()} has no
effect on a run. A run's random stream depends only on its own seed.

Functions:
\describe{
  \item{\code{\link{run_sim}}}{runs one replicate from a parameter list, a
    transition table and an initial state, and returns a data frame. This is
    the function to call from R code.}
  \item{\code{\link{simulate_hr}}}{the compiled function that
    \code{run_sim} calls, with every input as a separate vector or scalar.}
  \item{\code{\link{debug_step1}}, \code{\link{debug_step1_fired}}}{return
    the rates, transition probabilities and fired transitions of a single
    step. They exist for the package's tests.}
}

Each call runs one replicate on one thread; replicates can be run in
parallel by the caller.

The source package contains, in its \code{tests} folder, the R reference
implementation the engine was verified against and the scripts that make the
comparison. \code{citation("hrsim")} gives the reference of the paper.
}
\references{
Berryhill, B. A. and Andersson, D. I. (2026). Unstable Heteroresistance as a
Stepping Stone to Stable Resistance via Compensatory Mutations. bioRxiv.
doi:10.64898/2026.10.08.756632
}
\keyword{package}
