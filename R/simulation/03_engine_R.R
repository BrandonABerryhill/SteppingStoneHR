# 03_engine_R.R
#
# The R reference engine: one simulation of the model from an initial state to
# a given time. The compiled engine, hrsim::run_sim(), takes the same named
# arguments and returns the same columns. Its last two arguments are in the
# other order (stochastic, then n_out), which matters to a call that passes
# them by position, and its default for n_out differs: with it the compiled
# engine records the initial row and the transfer rows alone. 07_engines.R
# loads this file into an environment of its own and exposes its function as
# run_model_R, whichever engine is active. Sourcing this file directly defines
# run_model as the R reference.
#
# Numerical scheme. Supplementary Equations 3 to 25 are integrated with a
# fixed-step explicit Euler method (step p$DT, in hours). At every step a Monte
# Carlo layer draws one uniform random number per transition and moves one
# cell per mL for each transition that fires. Every p$transfer hours a serial
# transfer dilutes densities and antibiotic, adds the dose and resets r and B.
#
# Part of the code for "Unstable Heteroresistance as a Stepping Stone to Stable
# Resistance via Compensatory Mutations" (Berryhill and Andersson).
#
# Defines: run_model().

#' Run one simulation with the R reference engine
#'
#' @param p Parameter list, the element `p` of `make_condition()`: the scalars
#'   `DT`, `transfer`, `C`, `e`, `kmon`, `AMAX`, `dil`, `d`, `KB`, `mu`,
#'   `gamma` and `omega`, and the length-20 vectors `v`, `vm`, `kH`, `M` and
#'   `P`.
#' @param transitions Data frame from `build_transitions()`: `from` and `to`
#'   (positions 1 to 20), `type` ("gamma", "mu" or "omega") and `mult`.
#' @param y0 Initial state, the element `y0` of `make_condition()`: `pop`
#'   (named numeric vector of 20 densities, cells per mL), `r`, `A` and `B`.
#' @param t_end Numeric scalar, simulated time in hours.
#' @param seed Integer scalar passed to `set.seed()`. It governs every random
#'   draw of the run, so with a given kind of R generator a replicate's stream
#'   depends only on its own seed. It has no effect on the result when
#'   `stochastic = FALSE`.
#' @param n_out Positive integer, the approximate number of regularly spaced
#'   rows to record. Rows at transfers are recorded in addition.
#' @param stochastic Logical scalar. `FALSE` switches the Monte Carlo layer
#'   off, so that no transition fires and the run is a deterministic Euler
#'   integration; this is used to compare the two engines.
#' @return A data frame with 25 columns: `time` (h), `R0` to `R9` and `R0C` to
#'   `R9C` (cells per mL), `r` and `A` (µg per mL), `B` (units per mL) and
#'   `cycle_end`. The first row is the initial state. `cycle_end` is 1 on the
#'   row recorded at a transfer, which holds the state just before dilution,
#'   and 0 on every other row. The classifier uses the rows with
#'   `cycle_end == 1`.
run_model <- function(p, transitions, y0, t_end, seed = 1, n_out = 400,
                      stochastic = TRUE) {
  set.seed(seed)
  DT     <- p$DT
  nsteps <- round(t_end / DT)
  pop <- as.numeric(y0$pop); r <- y0$r; A <- y0$A; B <- y0$B

  v <- p$v; vm <- p$vm; kH <- p$kH; M <- p$M; Pb <- p$P
  # Rate of each transition per cell per hour: the base rate of its type times
  # its multiplier.
  base     <- c(gamma = p$gamma, mu = p$mu, omega = p$omega)
  rate_vec <- as.numeric(transitions$mult * base[transitions$type])
  from <- transitions$from; to <- transitions$to; ntr <- nrow(transitions)
  kmon <- p$kmon; eR <- p$e; d0 <- p$d; KB <- p$KB
  cyc  <- p$transfer; dil <- p$dil; AMAX <- p$AMAX; Cres <- p$C

  # A regular row is recorded every out_every steps. The output matrix has
  # room for those rows, one row per transfer and a small margin.
  out_every <- max(1, floor(nsteps / n_out))
  n_rec <- floor(nsteps / out_every) + ceiling(t_end / cyc) + 5
  out <- matrix(NA_real_, nrow = n_rec, ncol = 25)
  colnames(out) <- c("time", names(y0$pop), "r", "A", "B", "cycle_end")
  out[1, ] <- c(0, pop, r, A, B, 0); rec <- 2

  # tq is the time since the last transfer, in hours.
  tq <- 0
  for (step in 1:nsteps) {
    # psi is the resource limitation of Supplementary Equation 2, and Q the
    # net growth rate of each state, Supplementary Equation 1.
    psi <- r / (r + kmon)
    AM  <- (A / M)^kH
    HH  <- ((v - vm) * AM) / (AM - vm / v)
    Q   <- (v - HH) * psi

    # Rates of change, all evaluated on the state at the start of the step.
    # growth is the first term of Supplementary Equations 6 to 25, already
    # multiplied by DT; dr, dB and dA are Supplementary Equations 3 (resource),
    # 5 (β-lactamase) and 4 (antibiotic).
    growth <- DT * pop * Q
    dr <- -eR * psi * sum(v  * pop)
    dB <-  eR * psi * sum(Pb * pop)
    dA <- -A * d0 - KB * B * A

    # Monte Carlo layer. A transition fires when its uniform draw is below
    # rate * DT * density of its source state * psi. One number is drawn for
    # each row of the transition table at every step, whatever the
    # probability. A transition that fires moves one cell per mL from its
    # source state to its destination state.
    probs <- rate_vec * DT * pop[from] * psi
    if (stochastic) {
      fired <- runif(ntr) < probs
      net <- if (any(fired)) tabulate(to[fired], 20) - tabulate(from[fired], 20) else 0
    } else {
      net <- 0
    }

    # The antibiotic and the densities are floored at zero. The resource and
    # the β-lactamase are not.
    pop <- pop + growth + net
    r   <- r + DT * dr
    A   <- max(0, A + DT * dA)
    B   <- B + DT * dB
    pop[pop < 0] <- 0
    now <- step * DT

    # Serial transfer. tq is a running floating-point sum of DT, and with the
    # default step and interval it first reaches 24 one step after 24 h, so
    # the first transfer row has time 24.00001. The row is recorded before the
    # transfer. Then the densities and the antibiotic are diluted by dil, the
    # dose AMAX is added to the diluted antibiotic, the resource is reset to
    # its maximum and the β-lactamase to zero.
    tq <- tq + DT
    if (tq >= cyc) {
      if (rec <= n_rec) { out[rec, ] <- c(now, pop, r, A, B, 1); rec <- rec + 1 }
      pop <- pop * dil; A <- A * dil + AMAX; r <- Cres; B <- 0; tq <- tq - cyc
    }
    # Regular row. On a step that also made a transfer it holds the state
    # after the transfer.
    if (step %% out_every == 0 && rec <= n_rec) {
      out[rec, ] <- c(now, pop, r, A, B, 0); rec <- rec + 1
    }
  }
  as.data.frame(out[1:(rec - 1), , drop = FALSE])
}
