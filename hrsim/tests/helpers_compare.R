# helpers_compare.R
#
# Part of the tests of the hrsim package, the simulation engine for Berryhill
# and Andersson (2026), "Unstable Heteroresistance as a Stepping Stone to
# Stable Resistance via Compensatory Mutations".
# Defines compare_trajectories().

#' Compare two recorded trajectories cell by cell
#'
#' A cell passes when |test - ref| <= atol + rtol * |ref|. The absolute term
#' keeps cells whose reference value is zero or nearly zero from failing on a
#' relative criterion alone. The relative error that is reported is
#' |test - ref| / max(|ref|, atol).
#'
#' @param ref,test Numeric matrices with identical dimensions and column
#'   names, including a `time` column. The two time columns must agree to
#'   within 1e-6 h; otherwise the function stops, because the rows would not
#'   be the same time points.
#' @param rtol Relative tolerance, dimensionless.
#' @param atol Absolute tolerance, in the unit of each column.
#' @return A list: `pass` (TRUE when every cell other than time passes),
#'   `n_fail` (number of failing cells), `max_relerr`, and for the cell with
#'   the largest relative error `worst_time`, `worst_state` (its column name)
#'   and `worst_pair` (the reference and test values).
compare_trajectories <- function(ref, test, rtol = 1e-9, atol = 1e-12) {
  stopifnot(identical(dim(ref), dim(test)),
            identical(colnames(ref), colnames(test)))
  if (max(abs(ref[, "time"] - test[, "time"])) > 1e-6)
    stop("time grids differ — align n_out before diffing")

  state_cols <- setdiff(colnames(ref), "time")
  r <- ref[,  state_cols, drop = FALSE]
  t <- test[, state_cols, drop = FALSE]

  abserr <- abs(t - r)
  tol    <- atol + rtol * abs(r)
  ok     <- abserr <= tol
  relerr <- abserr / pmax(abs(r), atol)

  worst <- which.max(relerr)
  idx   <- arrayInd(worst, dim(r))
  list(
    pass        = all(ok),
    n_fail      = sum(!ok),
    max_relerr  = max(relerr),
    worst_time  = ref[idx[1], "time"],
    worst_state = state_cols[idx[2]],
    worst_pair  = c(ref = r[worst], test = t[worst])
  )
}
