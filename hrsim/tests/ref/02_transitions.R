# 02_transitions.R
#
# Table of the transitions between states: recombination, which changes the
# gene copy number, and compensation. Part of the R reference implementation
# the compiled engine was verified against. Belongs to Berryhill and
# Andersson (2026), "Unstable Heteroresistance as a Stepping Stone to Stable
# Resistance via Compensatory Mutations".
#
# Defines build_transitions().

#' Build the transition table
#'
#' Each transition has a type, which names its base rate, and a multiplier;
#' its rate constant is the multiplier times that base rate (p$gamma, p$mu or
#' p$omega). For a state with i copies (i = 1 to 9) the table holds:
#'   deletion to zero copies, type "gamma";
#'   doubling to 2 * i copies where 2 * i <= 9, type "gamma";
#'   a change to any other copy number j from 1 to min(2 * i, 9), type "mu",
#'   with multiplier i - |i - j|.
#' The same 65 recombination transitions are built for the non-compensated
#' states (indices 1 to 10) and for the compensated states (indices 11 to
#' 20). Each state is linked to its compensated counterpart by a pair of
#' "omega" transitions, one in each direction. That gives 150 rows.
#'
#' @return A data frame with one row per transition and columns `from` and
#'   `to` (1-based state indices, 1 to 20, in the order of `states`), `type`
#'   ("gamma", "mu" or "omega") and `mult` (the multiplier). The row order is
#'   the order in which the engines visit the transitions.
build_transitions <- function() {
  rows <- list()
  add <- function(from, to, type, mult)
    rows[[length(rows) + 1]] <<- data.frame(from = from, to = to,
                                            type = type, mult = mult,
                                            stringsAsFactors = FALSE)
  # offset is 0 for the non-compensated states and 10 for the compensated
  # ones; copy number i has state index i + 1 + offset.
  recomb <- function(offset) {
    for (i in 1:9) {
      add(i + 1 + offset, 1 + offset, "gamma", 1)
      hi <- min(2 * i, 9)
      for (j in setdiff(1:hi, i)) {
        if (j == 2 * i) add(i + 1 + offset, j + 1 + offset, "gamma", 1)
        else            add(i + 1 + offset, j + 1 + offset, "mu", i - abs(i - j))
      }
    }
  }
  recomb(0); recomb(10)
  for (i in 0:9) {
    add(i + 1, i + 11, "omega", 1)
    add(i + 11, i + 1, "omega", 1)
  }
  do.call(rbind, rows)
}
