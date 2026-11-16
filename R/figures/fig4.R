# fig4.R
#
# Draws Figure 4, the global sensitivity of the three outcome frequencies to
# the six swept parameters, as PAWN indices, and saves it as figures/Fig4.png.
#
# Panels and the panel data frames they are drawn from:
#   A  Fig4A  PAWN indices for the frequency of the climb outcome
#   B  Fig4B  the same for the low-copy outcome
#   C  Fig4C  the same for extinction
#
# For each parameter the 300 parameter sets of the sweep are split into up to
# eight quantile bins of that parameter (eight for every parameter of the
# shipped sweep). In each bin the Kolmogorov-Smirnov distance is taken between
# the distribution of the outcome frequency among its sets and that among all
# sets (pawn_given_data() in R/simulation/06_sensitivity.R). The bar is the
# median of that distance over the bins, and the tick is its maximum over them.
#
# The frames are read from source_data/panel_data.rds, which 01_panel_data.R
# writes; if that file is absent, source 01_panel_data.R first. The working
# directory must be the repository root. The panels carry no titles, because
# the figure legend of the paper says which outcome each one shows.
#
# Part of the code for "Unstable Heteroresistance as a Stepping Stone to Stable
# Resistance via Compensatory Mutations" (Berryhill and Andersson).
#
# Defines: PD, pawn_panel(), the panels pA, pB and pC, and the assembled
# figure fig4.

if (!exists("theme_hr")) source("R/figures/00_style.R")
PD <- readRDS(PANEL_DATA_RDS)$data

#' PAWN indices of one outcome as horizontal bars
#'
#' The parameters are drawn from top to bottom in the order of the rows of
#' `df`. The tick is plotting symbol 124, a vertical bar. The limits of the
#' x scale are 0 and 0.6, or 5% beyond the largest tick if that is further,
#' and the axis is drawn 2% of that range past the upper limit. The upper
#' limit is set panel by panel, so the three panels of the figure need not
#' share it.
#'
#' @param df One of the frames Fig4A to Fig4C: one row per parameter with
#'   `label` (axis label), `median` and `max` (median and maximum of the PAWN
#'   index over the quantile bins of the parameter, dimensionless, between 0
#'   and 1).
#' @param show_y Logical scalar; `FALSE` leaves out the parameter labels and
#'   the ticks of the y axis, for a panel placed to the right of one that
#'   shows them.
#' @return A ggplot object.
pawn_panel <- function(df, show_y = TRUE) {
  df$label <- factor(df$label, levels = rev(df$label))
  xmax <- max(0.6, max(df$max) * 1.05)
  p <- ggplot(df, aes(median, label)) +
    geom_col(width = 0.7, fill = unname(HR_TOL["indigo"])) +
    geom_point(aes(x = max), shape = 124, size = 2.2, colour = HR_PAL$ink) +
    scale_x_continuous(limits = c(0, xmax), expand = expansion(c(0, 0.02))) +
    labs(x = "PAWN index", y = NULL) +
    theme_hr()
  if (!show_y) p <- p + theme(axis.text.y = element_blank(),
                              axis.ticks.y = element_blank())
  p
}

# The three frames list the parameters in the same order, so the labels of
# panel A also serve panels B and C.
pA <- pawn_panel(PD$Fig4A, show_y = TRUE)
pB <- pawn_panel(PD$Fig4B, show_y = FALSE)
pC <- pawn_panel(PD$Fig4C, show_y = FALSE)

# Assemble and save ----

# The three panels sit in one row, and patchwork adds the panel tags A to C.
fig4 <- (pA | pB | pC) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text())

save_fig(fig4, "Fig4", width_mm = WIDTH_DOUBLE, height_mm = 62)

# The figure object is the last value of the script, which hands it to the
# caller and displays it when the script is run line by line.
fig4
