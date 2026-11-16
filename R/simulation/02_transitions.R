# 02_transitions.R
#
# The table of transitions between the twenty states. Each row is one
# transition that the Monte Carlo layer of the engines can fire: recombination
# between copy numbers within the non-compensated states and, identically,
# within the compensated states, and compensation between each state and its
# compensated counterpart. These are the terms M of Supplementary Equations 6
# to 25. The paper shows the recombination events and their relative rates in
# Supplementary Figure 1.
#
# Part of the code for "Unstable Heteroresistance as a Stepping Stone to Stable
# Resistance via Compensatory Mutations" (Berryhill and Andersson).
#
# Defines: build_transitions().

#' Build the transition table
#'
#' A state with i copies (i from 1 to 9) can lose all copies, double its copy
#' number when 2i is at most 9, or change stepwise to any copy number from 1 to
#' min(2i, 9) other than i and 2i. Every state can gain or lose the
#' compensatory mutation.
#'
#' @return A data frame of 150 rows, 65 for recombination among R0 to R9, 65
#'   for recombination among R0C to R9C and 20 for compensation, with columns
#'   `from` and `to` (positions in `states`, 1 to 20), `type` ("gamma", "mu" or
#'   "omega", the base rate that applies) and `mult` (the factor by which that
#'   base rate is multiplied). The engines go through the transitions in the
#'   row order of this table when they draw random numbers.
build_transitions <- function() {
  rows <- list()
  add <- function(from, to, type, mult)
    rows[[length(rows) + 1]] <<- data.frame(from = from, to = to,
                                            type = type, mult = mult,
                                            stringsAsFactors = FALSE)
  # offset is 0 for the non-compensated states and 10 for the compensated
  # ones, so i + 1 + offset is the position of the state with i copies.
  recomb <- function(offset) {
    for (i in 1:9) {
      # Deletion of every copy, at rate gamma.
      add(i + 1 + offset, 1 + offset, "gamma", 1)
      hi <- min(2 * i, 9)
      for (j in setdiff(1:hi, i)) {
        # Doubling, at rate gamma. Any other change from i to j copies is
        # stepwise, at rate mu times i - |i - j|.
        if (j == 2 * i) add(i + 1 + offset, j + 1 + offset, "gamma", 1)
        else            add(i + 1 + offset, j + 1 + offset, "mu", i - abs(i - j))
      }
    }
  }
  recomb(0); recomb(10)
  # Compensation and its reversal, both at rate omega.
  for (i in 0:9) {
    add(i + 1, i + 11, "omega", 1)
    add(i + 11, i + 1, "omega", 1)
  }
  do.call(rbind, rows)
}
