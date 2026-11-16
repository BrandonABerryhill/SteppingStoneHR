# fig2.R
#
# Draws Figure 2, the long-term dynamics of a heteroresistant population under
# Scenario 1, and saves it as figures/Fig2.png. Nothing shown is a single
# run: every curve is the mean density of a state over the 100 replicates of
# the 3,000 h Scenario 1 batch, and every heatmap cell is the replicate-mean
# density of a compensated state divided by the replicate-mean total density.
#
# Panels and the panel data frames they are drawn from:
#   A  Fig2A  mean density of the states, 0 to 72 h
#   B  Fig2B  the same, 0 to 1,000 h
#   C  Fig2C  the same, 0 to 3,000 h
#   D  Fig2D  heatmap of each compensated state as a fraction of the total
#             population, 0 to 3,000 h
#
# Panels A to C do not draw all 20 states: a state is drawn in a panel only if
# it reaches 0.5% of the peak total population within that panel's window (see
# keep_states() below). The legend of the three panels is taken from panel A.
#
# The frames are read from source_data/panel_data.rds, which 01_panel_data.R
# writes; if that file is absent, source 01_panel_data.R first. The working
# directory must be the repository root. The panels carry no titles; each
# window is read from its x axis and from the figure legend of the paper.
#
# Part of the code for "Unstable Heteroresistance as a Stepping Stone to Stable
# Resistance via Compensatory Mutations" (Berryhill and Andersson).
#
# Defines: PD, keep_states(), copy_levels, dyn_panel(), the panels pA, pB, pC
# and pD, and the assembled figure fig2.

if (!exists("theme_hr")) source("R/figures/00_style.R")
suppressPackageStartupMessages(library(dplyr))
PD <- readRDS(PANEL_DATA_RDS)$data

# Panels A, B and C: mean densities over a window ----

#' States drawn in one window
#'
#' A state is kept if its largest mean density within the window is at least
#' `thresh_frac` of the largest total density, summed over the 20 states,
#' found at any time in the window. States that stay below it are left out
#' so that the lines of the abundant states remain legible.
#'
#' @param df One of the frames Fig2A to Fig2C, already cut to its window: `t`
#'   (h), `state`, `density` (mean density, cells per mL), `copy` and
#'   `compartment`, one row per time and state.
#' @param thresh_frac Numeric scalar between 0 and 1, the threshold as a
#'   fraction of the peak total density.
#' @return A tibble with one row per kept state: `state`, `mx` (its largest
#'   mean density in the window, cells per mL) and `copy`.
keep_states <- function(df, thresh_frac = 0.005) {
  maxtot <- max(tapply(df$density, df$t, sum), na.rm = TRUE)
  df |>
    group_by(state) |>
    summarise(mx = max(density, na.rm = TRUE), copy = first(copy), .groups = "drop") |>
    filter(mx >= thresh_frac * maxtot)
}

# With the shipped results the peak total density is 1.95e9 cells per mL in
# all three windows, which puts the threshold at 9.7e6 cells per mL. Panel A
# then draws R2 to R9 and R2C, and panels B and C draw R2 to R9 and R1C to
# R9C. R0, R1 and R0C stay below the threshold in every window and are drawn
# in no panel.

# Copy numbers given to the colour scale as its limits: those of the states
# kept in at least one of the three windows. With the shipped results these
# are 1 to 9, so the legend has no entry for copy number 0.
copy_levels <- as.character(sort(unique(unlist(lapply(
  list(PD$Fig2A, PD$Fig2B, PD$Fig2C),
  function(d) keep_states(d)$copy)))))

#' Mean densities of the kept states over one window
#'
#' Rows with a mean density of zero are left out because the y axis is
#' logarithmic. With the shipped results these are rows of compensated states
#' at 0 h and, for R6C to R9C, at 12 h, so the line of a compensated state
#' starts at 12 h or at 24 h.
#'
#' Only a panel drawn with `legend = TRUE` has guides, and the figure gives
#' them to panel A alone; `plot_layout(guides = "collect")` then places them
#' beside the figure, which leaves one copy-number legend and one line-type
#' legend for the three panels. Every panel uses the same colour scale with
#' the same limits, so a copy number has the same colour in each. A legend
#' key is drawn only for a copy number that panel A itself draws: with the
#' shipped results copy number 1 passes the threshold only as R1C in panels B
#' and C, so its legend entry has a label and an empty key. The colour keys
#' are set to solid lines so that line type is read from its own legend.
#'
#' @param df One of the frames Fig2A to Fig2C (see `keep_states()`).
#' @param thresh_frac Numeric scalar passed to `keep_states()`.
#' @param legend Logical scalar; `TRUE` gives the panel the colour and
#'   line-type legends, `FALSE` gives it none.
#' @return A ggplot object with one line per kept state, coloured by copy
#'   number, solid for non-compensated and dashed for compensated states.
dyn_panel <- function(df, thresh_frac = 0.005, legend = FALSE) {
  keep <- keep_states(df, thresh_frac)$state
  d <- df[df$state %in% keep & df$density > 0, ]
  g <- if (legend) "legend" else "none"
  p <- ggplot(d, aes(t, density, colour = factor(copy), group = state,
                     linetype = compartment)) +
    geom_line(linewidth = 0.4) +
    scale_y_log10() +
    scale_copy_colour(limits = copy_levels, drop = FALSE, guide = g) +
    scale_linetype_manual(values = compartment_linetypes, name = NULL, guide = g) +
    labs(x = "Time (h)", y = "Density (cells/mL)") +
    theme_hr()
  if (legend)
    p <- p + guides(colour = guide_legend(override.aes = list(linetype = "solid")))
  p
}
pA <- dyn_panel(PD$Fig2A, legend = TRUE)
pB <- dyn_panel(PD$Fig2B)
pC <- dyn_panel(PD$Fig2C)

# Panel D: heatmap of the compensated states ----

# Each cell is the mean density of the compensated state with that copy number
# divided by the total mean density of all 20 states, at one time of the 12 h
# output grid. Non-compensated states have no row, so the ten cells at one
# time add up to the compensated fraction of the population and not to one.
#
# The colour scale ends at 0.3, and oob = scales::squish gives every value
# above 0.3 the colour of 0.3, which keeps the smaller values apart. With the
# shipped results 35 cells lie above 0.3, all at copy number 2: one of 0.76
# at 12 h, where the total mean density is at its lowest, and 34 between
# 0.301 and 0.315 from 2,604 h onwards.
pD <- ggplot(PD$Fig2D, aes(t, copy, fill = freq)) +
  geom_raster() +
  scale_freq_fill(limits = c(0, 0.3), breaks = c(0, 0.1, 0.2, 0.3),
                  oob = scales::squish) +
  scale_y_continuous(breaks = c(0, 3, 6, 9), expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Time (h)", y = "Copy number") +
  theme_hr()

# Assemble and save ----

# A, B and C share the top row and D spans the bottom row. The collected
# legends, panel A's two and the colour bar of D, are placed to the right of
# the figure. patchwork adds the panel tags A to D.
fig2 <- ((pA | pB | pC) / pD) +
  plot_layout(guides = "collect", heights = c(1, 1.05)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(), legend.position = "right")

save_fig(fig2, "Fig2", width_mm = WIDTH_DOUBLE, height_mm = 120)

# The figure object is the last value of the script, which hands it to the
# caller and displays it when the script is run line by line.
fig2
