# figS5.R
#
# Draws Supplementary Figure 5, the raw parameter sweep, and saves it as
# figures/FigS5.png.
#
# Six panels, all drawn from the panel data frame FigS5, which has one row per
# parameter set of the sweep (300 sets): the six sampled parameters and the
# frequency of each outcome among the replicates of the set. Each panel puts
# one parameter on the x axis:
#   A  amax           antibiotic dose
#   B  omega          compensation rate
#   C  cost_per_copy  per-copy cost
#   D  mu             recombination rate
#   E  dil            serial transfer rate
#   F  comp_frac      fraction of the cost recovered by compensation
#
# The x axis follows the way the parameter was sampled: logarithmic for the
# compensation rate, the recombination rate and the serial transfer rate, and
# linear for the other three.
#
# The frame is read from source_data/panel_data.rds, which 01_panel_data.R
# writes; if that file is absent, source 01_panel_data.R first. The working
# directory must be the repository root. The panels carry no titles; the
# parameter is named by the x-axis label.
#
# Part of the code for "Unstable Heteroresistance as a Stepping Stone to Stable
# Resistance via Compensatory Mutations" (Berryhill and Andersson).
#
# Defines: PD, sweep_panel(), d, the panels pDose, pOmega, pCost, pMu, pDil
# and pFrac, and the assembled figure figS5.

if (!exists("theme_hr")) source("R/figures/00_style.R")
suppressPackageStartupMessages(library(tidyr))
PD <- readRDS(PANEL_DATA_RDS)$data

#' Outcome frequencies of the sweep against one parameter
#'
#' Each parameter set gives three points, one per outcome, and each outcome
#' has a loess curve fitted to its 300 points (`span = 0.9`, no confidence
#' band). On a logarithmic x axis the curve is fitted to log10 of the
#' parameter.
#'
#' The y axis is limited to the range of a frequency, 0 to 1. Every point
#' lies within it, so no point is left out. A fitted curve is not confined
#' to that range and can pass below 0 or above 1. The curve is evaluated at
#' 80 values of x, ggplot2 removes the fitted values below 0 or above 1 with
#' a warning, and the line is drawn through the values that remain. A curve
#' that leaves the range at one end of the x axis therefore stops short of
#' that end, and one that leaves it in mid-course is closed by a straight
#' segment between the fitted values either side of the gap.
#'
#' @param df The frame FigS5: one row per parameter set with the six
#'   parameter columns and `f_climb`, `f_lowcopy` and `f_extinct`, the
#'   frequency of each outcome among the replicates of the set (0 to 1).
#' @param param Character scalar, name of the parameter column put on the x
#'   axis.
#' @param xlab Character scalar, x-axis label.
#' @param logx Logical scalar; `TRUE` draws the x axis on a log10 scale.
#' @return A ggplot object.
sweep_panel <- function(df, param, xlab, logx = FALSE) {
  long <- pivot_longer(df, c(f_climb, f_lowcopy, f_extinct),
                       names_to = "outcome", values_to = "freq")
  long$outcome <- factor(long$outcome,
                         levels = c("f_climb", "f_lowcopy", "f_extinct"),
                         labels = c("Climb", "Low-copy", "Extinction"))
  long$x <- long[[param]]
  p <- ggplot(long, aes(x, freq, colour = outcome)) +
    geom_point(alpha = 0.22, size = 0.5) +
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE,
                linewidth = 0.5, span = 0.9) +
    scale_colour_manual(values = outcome_colors, name = NULL) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    labs(x = xlab, y = "Frequency") +
    theme_hr()
  if (logx) p <- p + scale_x_log10()
  p
}

# With the shipped sweep, fitted values are removed in two of the six panels,
# and drawing the figure gives one warning for each.
#
# Panel A (dose), 71 of the 240 fitted values:
#   Climb       31 below 0 (lowest -0.057), doses 35.1 to 57.8; straight
#               segment across the gap
#   Low-copy    13 below 0 (lowest -0.020), dose 54.8 to the upper end; the
#               line stops at 54.0
#   Extinction  11 below 0 (lowest -0.065), lower end to dose 11.6; the line
#               starts at 12.4
#               16 above 1 (highest 1.028), doses 48.7 to 60.1; straight
#               segment across the gap
#
# Panel B (compensation rate), 23 of the 240 fitted values:
#   Climb       8 below 0 (lowest -0.046), 4.3e-7 to the upper end; the line
#               stops at 3.8e-7
#   Low-copy    15 below 0 (lowest -0.011), 2.9e-10 to 1.5e-9; straight
#               segment across the gap
d <- PD$FigS5
pDose  <- sweep_panel(d, "amax",          "Dose (AMAX)",                 logx = FALSE)
pOmega <- sweep_panel(d, "omega",         "Compensation rate (\u03c9)",  logx = TRUE)
pCost  <- sweep_panel(d, "cost_per_copy", "Per-copy cost",               logx = FALSE)
pMu    <- sweep_panel(d, "mu",            "Recombination rate (\u03bc)", logx = TRUE)
pDil   <- sweep_panel(d, "dil",           "Serial transfer rate",        logx = TRUE)
pFrac  <- sweep_panel(d, "comp_frac",     "Recovered fraction",          logx = FALSE)

# Assemble and save ----

# Two rows of three panels. The six identical outcome legends are collected
# into one, placed below the figure, and patchwork adds the panel tags A to F.
figS5 <- (pDose | pOmega | pCost) / (pMu | pDil | pFrac) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(), legend.position = "bottom")

save_fig(figS5, "FigS5", width_mm = WIDTH_DOUBLE, height_mm = 112)

# The figure object is the last value of the script, which hands it to the
# caller and displays it when the script is run line by line.
figS5
