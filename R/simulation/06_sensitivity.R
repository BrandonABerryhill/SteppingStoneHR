# 06_sensitivity.R
#
# The sensitivity analysis: a Latin hypercube sample of parameter sets, the
# outcome frequencies of each set, and PAWN indices. The quantity analysed is
# the frequency of each outcome, not the mean copy number: every parameter set
# is run for many seeds, each replicate is classified, and the set is reduced
# to its three frequencies. PAWN then measures how strongly each parameter
# shifts the distribution of one of those frequencies.
#
# The functions take whatever parameters the caller supplies, provided
# make_condition() accepts them by name. The paper's sweep varies six: omega,
# amax, cost_per_copy, mu, dil and comp_frac.
#
# Part of the code for "Unstable Heteroresistance as a Stepping Stone to Stable
# Resistance via Compensatory Mutations" (Berryhill and Andersson).
#
# Defines: lhs_sample(), freq_for_set(), run_sweep() and pawn_given_data().

# Sampling design ----

#' Latin hypercube sample of parameter sets
#'
#' The sample is drawn by `lhs::randomLHS()` from R's random number generator.
#' The caller sets the seed immediately before the call, and that seed governs
#' the whole design. The design also depends on the kind of R's generator in
#' force; the design of the paper was drawn with R's default kind,
#' Mersenne-Twister.
#'
#' @param n Positive integer, number of parameter sets.
#' @param ranges Named list with one element per parameter, each `c(min, max)`
#'   in the parameter's own unit.
#' @param log Character vector naming the parameters sampled uniformly on a
#'   log10 scale; they need positive bounds. The others are sampled uniformly
#'   on a linear scale.
#' @return A data frame of `n` rows with one numeric column per parameter,
#'   named as in `ranges`.
lhs_sample <- function(n, ranges, log = c("omega")) {
  X <- lhs::randomLHS(n, length(ranges))
  colnames(X) <- names(ranges)
  for (j in names(ranges)) {
    lo <- ranges[[j]][1]; hi <- ranges[[j]][2]
    if (j %in% log) X[, j] <- 10^(log10(lo) + X[, j] * (log10(hi) - log10(lo)))
    else            X[, j] <- lo + X[, j] * (hi - lo)
  }
  as.data.frame(X)
}

# Outcome frequencies of the sampled sets ----

#' Outcome frequencies of one parameter set
#'
#' @param scenario Character scalar passed to `make_condition()`.
#' @param pars Named list, or one-row data frame, of parameter values. The
#'   names are passed to `make_condition()` as they are.
#' @param seeds Integer vector of replicate seeds; each replicate's random
#'   stream depends only on its own seed.
#' @param t_end Numeric scalar, simulated time in hours; by default
#'   `CLASSIFY_HORIZON`.
#' @param n_workers Positive integer, parallel workers; does not affect
#'   results.
#' @param set_tag Character scalar, the batch tag of this set. It must differ
#'   between sets, because `run_batch()` resumes from the batch stored under a
#'   tag.
#' @return Named numeric vector `c(f_climb, f_lowcopy, f_extinct)`, the
#'   fraction (0 to 1) of the set's successful replicates with each outcome;
#'   all `NA` if no replicate succeeded.
freq_for_set <- function(scenario, pars, seeds,
                         t_end = CLASSIFY_HORIZON, n_workers = 4,
                         set_tag = "sweep") {
  cond <- do.call(make_condition, c(list(scenario = scenario), as.list(pars)))
  b   <- run_batch(cond, t_end = t_end, seeds = seeds,
                   n_workers = n_workers, tag = set_tag, manifest = FALSE)
  if (length(b$sims) == 0)
    return(setNames(c(NA, NA, NA), c("f_climb", "f_lowcopy", "f_extinct")))
  cls <- do.call(rbind, lapply(b$sims, classify_outcome))
  s   <- summarise_frequencies(cls)
  setNames(s$percent / 100, c("f_climb", "f_lowcopy", "f_extinct"))
}

#' Outcome frequencies of every parameter set
#'
#' Set i is stored as the batch tagged "sweep_set" followed by i padded to four
#' digits. Finished sets are also kept in `ckpt_path`, keyed by set number, and
#' are not run again. A resumed sweep therefore needs the same `X`, row for
#' row, which the caller ensures by drawing it from a fixed seed.
#'
#' @param scenario Character scalar passed to `make_condition()`.
#' @param X Data frame of parameter sets, one row per set and one column per
#'   parameter, as returned by `lhs_sample()`.
#' @param seeds_per_set Positive integer, replicates per set. The seeds are 1
#'   to `seeds_per_set` in every set, and a replicate's random stream depends
#'   only on its own seed.
#' @param n_workers Positive integer, parallel workers; does not affect
#'   results.
#' @param ckpt_path Path of the `.rds` file that holds the frequencies of the
#'   finished sets, or `NULL` to keep no such file.
#' @return `X` with three columns added: `f_climb`, `f_lowcopy` and
#'   `f_extinct`, the outcome frequencies (0 to 1) of each set.
run_sweep <- function(scenario, X, seeds_per_set, n_workers = 4,
                      ckpt_path = NULL) {
  done <- list()
  if (!is.null(ckpt_path) && file.exists(ckpt_path)) {
    done <- readRDS(ckpt_path)
    message(sprintf("  resuming sweep: %d/%d sets already on disk",
                    length(done), nrow(X)))
  }
  fr <- vector("list", nrow(X))
  for (i in seq_len(nrow(X))) {
    key <- sprintf("set%04d", i)
    if (!is.null(done[[key]])) { fr[[i]] <- done[[key]]; next }
    pars <- as.list(X[i, , drop = FALSE])
    fr[[i]] <- freq_for_set(scenario, pars, seeds = seq_len(seeds_per_set),
                            n_workers = n_workers, set_tag = paste0("sweep_", key))
    done[[key]] <- fr[[i]]
    if (!is.null(ckpt_path)) saveRDS(done, ckpt_path)
    message(sprintf("set %d/%d  %s -> f_lowcopy=%.3f", i, nrow(X),
                    paste(sprintf("%s=%.3g", names(pars), unlist(pars)), collapse = " "),
                    fr[[i]]["f_lowcopy"]))
  }
  cbind(X, do.call(rbind, fr[seq_len(nrow(X))]))
}

# PAWN indices ----

#' PAWN sensitivity indices, given-data estimator
#'
#' For each input the sample is split into bins at the input's quantiles, and
#' the Kolmogorov-Smirnov distance is taken between the distribution of the
#' output within each bin and its distribution over the whole sample. The
#' distance is evaluated at the observed output values. A larger index means a
#' more influential input. The estimator is that of Pianosi and Wagener,
#' Environmental Modelling & Software 108, 197-207 (2018).
#'
#' @param X Data frame of inputs, one row per parameter set and one column per
#'   parameter.
#' @param y Numeric vector of the output, one value per row of `X`.
#' @param n_bins Positive integer, number of quantile bins per input.
#' @return A data frame with one row per input, named by the columns of `X`,
#'   and columns `median` and `max`: the median and the maximum of the
#'   distance over the bins that hold at least three sets. Both are `NA` for
#'   an input with fewer than three distinct bin edges.
pawn_given_data <- function(X, y, n_bins = 8) {
  Fall <- ecdf(y)
  res <- lapply(colnames(X), function(j) {
    x  <- X[[j]]
    br <- unique(quantile(x, probs = seq(0, 1, length.out = n_bins + 1),
                          na.rm = TRUE))
    if (length(br) < 3) return(c(median = NA, max = NA))
    bin <- cut(x, breaks = br, include.lowest = TRUE)
    ks  <- tapply(seq_along(y), bin, function(idx) {
      if (length(idx) < 3) return(NA_real_)
      max(abs(ecdf(y[idx])(sort(y)) - Fall(sort(y))))
    })
    ks <- ks[is.finite(ks)]
    c(median = stats::median(ks), max = max(ks))
  })
  out <- do.call(rbind, res); rownames(out) <- colnames(X)
  as.data.frame(out)
}
