# 03_engine_R.R
#
# R reference engine: the R implementation of the model that the compiled
# engine in src/simulate_hr.cpp was verified against. Belongs to Berryhill
# and Andersson (2026), "Unstable Heteroresistance as a Stepping Stone to
# Stable Resistance via Compensatory Mutations"; Supplementary Equations
# (S. Eq) are those of that paper.
#
# Defines run_model().
#
# Numerical scheme: fixed-step explicit Euler for the populations, the
# resource, the antibiotic and beta-lactamase; a Monte Carlo layer in which
# every transition draws one uniform random number per step; and a serial
# transfer every p$transfer hours.

#' Run one replicate with the R reference engine
#'
#' @param p Parameter list from make_condition()$p: DT (Euler step, hours),
#'   transfer (hours between transfers), C, e, kmon, AMAX, dil, d, KB, mu,
#'   gamma, omega, and the length-20 vectors v, vm, kH, M and P.
#' @param transitions Data frame from build_transitions(): from and to
#'   (1-based state indices, 1 to 20), type ("gamma", "mu" or "omega") and
#'   mult.
#' @param y0 Initial state: list(pop, r, A, B), with pop the named densities
#'   of the 20 states in cells per mL.
#' @param t_end Simulated time, hours.
#' @param seed Seed passed to set.seed() at the start of the run. It governs
#'   every random draw of the run, so a replicate's random stream depends
#'   only on its own seed and on the kind of generator R is using
#'   (RNGkind()). It has no effect on the result when stochastic = FALSE.
#' @param n_out Approximate number of evenly spaced rows to record in
#'   addition to the transfer rows; a grid row is recorded every
#'   max(1, floor(nsteps / n_out)) steps.
#' @param stochastic If FALSE the Monte Carlo layer is disabled and no
#'   transition fires.
#' @return A data frame with 25 columns: time (hours), R0..R9 and R0C..R9C
#'   (cells per mL), r, A, B and cycle_end. The first row is the initial
#'   state. cycle_end is 1 on the row recorded just before each serial
#'   transfer, before dilution, and 0 on every other row.
run_model <- function(p, transitions, y0, t_end, seed = 1, n_out = 400,
                      stochastic = TRUE) {
  set.seed(seed)
  DT     <- p$DT
  nsteps <- round(t_end / DT)
  pop <- as.numeric(y0$pop); r <- y0$r; A <- y0$A; B <- y0$B

  v <- p$v; vm <- p$vm; kH <- p$kH; M <- p$M; Pb <- p$P
  base     <- c(gamma = p$gamma, mu = p$mu, omega = p$omega)
  rate_vec <- as.numeric(transitions$mult * base[transitions$type])
  from <- transitions$from; to <- transitions$to; ntr <- nrow(transitions)
  kmon <- p$kmon; eR <- p$e; d0 <- p$d; KB <- p$KB
  cyc  <- p$transfer; dil <- p$dil; AMAX <- p$AMAX; Cres <- p$C

  # n_rec is the number of rows allocated, an upper bound on the number
  # returned: the grid rows, ceiling(t_end / cyc) for the transfer rows, and a
  # margin of 5 that covers the row at time 0.
  out_every <- max(1, floor(nsteps / n_out))
  n_rec <- floor(nsteps / out_every) + ceiling(t_end / cyc) + 5
  out <- matrix(NA_real_, nrow = n_rec, ncol = 25)
  colnames(out) <- c("time", names(y0$pop), "r", "A", "B", "cycle_end")
  out[1, ] <- c(0, pop, r, A, B, 0); rec <- 2

  # tq is the transfer clock, in hours: the sum of DT over the steps taken so
  # far, less the transfer interval (cyc) for each transfer carried out.
  tq <- 0
  for (step in 1:nsteps) {
    # Everything from here to the update is computed from the state at the
    # start of the step. psi is S. Eq 2 and Q the net growth rate of S. Eq 1.
    psi <- r / (r + kmon)
    AM  <- (A / M)^kH
    HH  <- ((v - vm) * AM) / (AM - vm / v)
    Q   <- (v - HH) * psi

    # Euler increment of each population, and the time derivatives of the
    # resource (S. Eq 3), beta-lactamase (S. Eq 5) and antibiotic (S. Eq 4).
    growth <- DT * pop * Q
    dr <- -eR * psi * sum(v  * pop)
    dB <-  eR * psi * sum(Pb * pop)
    dA <- -A * d0 - KB * B * A

    # Monte Carlo layer. A transition fires when its uniform draw is below
    # rate * DT * (density of the source state) * psi. When the layer is
    # enabled, one draw is made for every transition at every step. Each
    # firing moves one cell per mL from the source state to the target state.
    probs <- rate_vec * DT * pop[from] * psi
    if (stochastic) {
      fired <- runif(ntr) < probs
      net <- if (any(fired)) tabulate(to[fired], 20) - tabulate(from[fired], 20) else 0
    } else {
      net <- 0
    }

    # The antibiotic and the populations are floored at zero; the resource
    # and beta-lactamase are not.
    pop <- pop + growth + net
    r   <- r + DT * dr
    A   <- max(0, A + DT * dA)
    B   <- B + DT * dB
    pop[pop < 0] <- 0
    now <- step * DT

    # Serial transfer. The clock is advanced by adding DT at every step. With
    # DT = 1e-5 h and a 24 h interval the sum after 2,400,000 additions is
    # just below 24, so the first transfer falls at step 2,400,001, at a
    # recorded time of 24.00001 h. The row is recorded before dilution. The
    # interval is subtracted from the clock, so the remainder carries over.
    tq <- tq + DT
    if (tq >= cyc) {
      if (rec <= n_rec) { out[rec, ] <- c(now, pop, r, A, B, 1); rec <- rec + 1 }
      pop <- pop * dil; A <- A * dil + AMAX; r <- Cres; B <- 0; tq <- tq - cyc
    }
    # Grid row; on a transfer step it holds the state after the transfer.
    if (step %% out_every == 0 && rec <= n_rec) {
      out[rec, ] <- c(now, pop, r, A, B, 0); rec <- rec + 1
    }
  }
  as.data.frame(out[1:(rec - 1), , drop = FALSE])
}
