# 04_classify.R
#
# Classification of a replicate into one of three outcomes, and the outcome
# frequencies of a batch with their confidence intervals. The rule is the one
# given under "Route Classification" in the Methods: extinction if the total
# density on the last row is below EXTINCT_BELOW; otherwise climb if the peak
# mean copy number reached CLIMB_THRESHOLD; otherwise low-copy. Extinction is
# tested first because the mean copy number of an empty population is
# undefined. Both constants are set in 00_config.R.
#
# Part of the code for "Unstable Heteroresistance as a Stepping Stone to Stable
# Resistance via Compensatory Mutations" (Berryhill and Andersson).
#
# Defines: mean_copy_number(), classify_outcome(), summarise_frequencies(),
# subspace_summary(), cp_ci(), classify_batch() and freq_table().

# Outcome of one replicate ----

#' Mean copy number of each row of a simulation
#'
#' @param d Data frame with the 20 state columns named in `states` (cells per
#'   mL).
#' @return Numeric vector with one value per row: the sum of copy number times
#'   density divided by the total density, or `NA` where the total density is
#'   zero.
mean_copy_number <- function(d) {
  cells <- as.matrix(d[, states])
  N     <- rowSums(cells)
  ifelse(N > 0, as.numeric(cells %*% copy_number) / N, NA_real_)
}

#' Classify one replicate
#'
#' The peak is taken over every row of `d`, and the final density from the row
#' with the largest time; a final density that is missing counts as
#' extinction. To classify on the transfer rows, as in the paper, pass only
#' those rows; `classify_batch()` does this.
#'
#' @param d Data frame of one replicate with columns `time` (h), the 20 state
#'   columns (cells per mL) and `seed`.
#' @param thresh Numeric scalar, the peak mean copy number from which a
#'   surviving replicate counts as a climb.
#' @param extinct_below Numeric scalar, cells per mL. A final total density
#'   below it is an extinction.
#' @return A data frame of one row with columns `seed`, `peak_mACN` (peak mean
#'   copy number), `final_N` (total density on the last row, cells per mL),
#'   `outcome` ("extinction", "climb" or "low-copy") and `nonfinite` (`TRUE` if
#'   any state value is not finite).
classify_outcome <- function(d,
                             thresh        = CLIMB_THRESHOLD,
                             extinct_below = EXTINCT_BELOW) {
  cells  <- as.matrix(d[, states])
  N      <- rowSums(cells)
  macn   <- ifelse(N > 0, as.numeric(cells %*% copy_number) / N, NA_real_)
  finalN <- N[which.max(d$time)]
  peak   <- if (all(is.na(macn))) NA_real_ else max(macn, na.rm = TRUE)

  outcome <- if (is.na(finalN) || finalN < extinct_below) {
    "extinction"
  } else if (!is.na(peak) && peak >= thresh) {
    "climb"
  } else {
    "low-copy"
  }
  data.frame(seed = d$seed[1], peak_mACN = peak, final_N = finalN,
             outcome = outcome, nonfinite = any(!is.finite(cells)),
             stringsAsFactors = FALSE)
}

# Outcome frequencies of a batch ----

#' Outcome frequencies with exact binomial confidence intervals
#'
#' @param cls Data frame of classifications, one row per replicate, as
#'   returned by `classify_batch()`.
#' @return A data frame of three rows, in the order climb, low-copy,
#'   extinction, with columns `outcome`, `n` (count), `percent`, `ci_low` and
#'   `ci_high` (95% Clopper-Pearson bounds, in percent). Nothing is rounded.
#'   The attribute `n_total` holds the number of replicates and `nonfinite`
#'   the number of replicates with a non-finite state value.
summarise_frequencies <- function(cls) {
  n   <- nrow(cls)
  lev <- c("climb", "low-copy", "extinction")
  out <- lapply(lev, function(o) {
    k  <- sum(cls$outcome == o)
    ci <- if (n > 0) binom.test(k, n)$conf.int else c(NA, NA)
    data.frame(outcome = o, n = k,
               percent = if (n > 0) 100 * k / n else NA_real_,
               ci_low  = 100 * ci[1], ci_high = 100 * ci[2],
               stringsAsFactors = FALSE)
  })
  res <- do.call(rbind, out)
  attr(res, "n_total")   <- n
  attr(res, "nonfinite") <- sum(cls$nonfinite)
  res
}

# Non-compensated and compensated states separately ----

#' Density and mean copy number of the two groups of states
#'
#' A single mean copy number over all twenty states cannot show the population
#' moving from the non-compensated to the compensated states, so the two groups
#' are summarised separately.
#'
#' @param d Data frame of one replicate with columns `time`, `R0` to `R9` and
#'   `R0C` to `R9C`.
#' @return A data frame with one row per row of `d` and columns `time`, `N_R`
#'   and `N_C` (total density of the non-compensated and of the compensated
#'   states, cells per mL), `mACN_R` and `mACN_C` (mean copy number within each
#'   group, `NA` where the group is empty) and `frac_C` (compensated fraction
#'   of the total density, `NA` where the total is zero).
subspace_summary <- function(d) {
  R <- as.matrix(d[, paste0("R", 0:9)])
  C <- as.matrix(d[, paste0("R", 0:9, "C")])
  cn <- 0:9
  N_R <- rowSums(R); N_C <- rowSums(C)
  data.frame(
    time = d$time,
    N_R = N_R, N_C = N_C,
    mACN_R = ifelse(N_R > 0, as.numeric(R %*% cn) / N_R, NA_real_),
    mACN_C = ifelse(N_C > 0, as.numeric(C %*% cn) / N_C, NA_real_),
    frac_C = ifelse((N_R + N_C) > 0, N_C / (N_R + N_C), NA_real_)
  )
}

# Batch helpers ----

#' Clopper-Pearson 95% confidence interval of a proportion
#'
#' @param k Integer scalar, number of replicates with the outcome.
#' @param n Integer scalar, number of replicates; at least 1.
#' @return Named numeric vector `c(low, high)`, in percent.
cp_ci <- function(k, n) {
  ci <- stats::binom.test(k, n)$conf.int
  c(low = 100*ci[1], high = 100*ci[2])
}

#' Classify every replicate of a batch on its transfer rows
#'
#' @param sims List of data frames, one per replicate, as in the element `sims`
#'   of `run_batch()`. Each has the engine's 25 columns and `seed`.
#' @return A data frame with one row per replicate and the columns of
#'   `classify_outcome()`.
classify_batch <- function(sims) {
  do.call(rbind, lapply(sims, function(d) {
    ce <- d[d$cycle_end == 1, , drop = FALSE]
    # A run that ends before its first transfer has no transfer row and is
    # classified on all of its rows.
    if (nrow(ce) == 0) ce <- d
    ce$seed <- d$seed[1]
    classify_outcome(ce)
  }))
}

#' Outcome frequencies rounded for display
#'
#' @param cls Data frame of classifications, one row per replicate, as
#'   returned by `classify_batch()`; at least one row.
#' @return A data frame of three rows, in the order climb, low-copy,
#'   extinction, with columns `outcome`, `n` (count), `pct` (percent) and
#'   `ci_low` and `ci_high` (95% Clopper-Pearson bounds, in percent). The last
#'   three are rounded to three decimals.
freq_table <- function(cls) {
  n <- nrow(cls)
  outs <- c("climb","low-copy","extinction")
  do.call(rbind, lapply(outs, function(o) {
    k <- sum(cls$outcome == o); ci <- cp_ci(k, n)
    data.frame(outcome = o, n = k, pct = round(100*k/n, 3),
               ci_low = round(ci["low"],3), ci_high = round(ci["high"],3))
  }))
}
