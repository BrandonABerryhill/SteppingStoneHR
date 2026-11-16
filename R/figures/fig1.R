# fig1.R
#
# Draws Figure 1, the outcomes of selection on a heteroresistant population
# under Scenario 1, and saves it as figures/Fig1.png.
#
# Panels and the panel data frames they are drawn from:
#   A  Fig1A  mean copy number of one replicate that climbs to high copy number
#   B  Fig1B  mean copy number of one replicate that stays at low copy number
#   C  Fig1C  frequency of each outcome among the 4,000 compiled replicates,
#             with 95% Clopper-Pearson intervals
#   D  Fig1D  control without compensation: share of each non-compensated
#             state in the total population at 288 h
#
# The frames are read from source_data/panel_data.rds, which 01_panel_data.R
# writes; if that file is absent, source 01_panel_data.R first. The working
# directory must be the repository root. The panels carry no titles, because
# the figure legend of the paper describes them.
#
# Part of the code for "Unstable Heteroresistance as a Stepping Stone to Stable
# Resistance via Compensatory Mutations" (Berryhill and Andersson).
#
# Defines: PD, mc_panel(), d1c, d1d, floorY, the panels pA, pB, pC and pD,
# bottom and the assembled figure fig1.

if (!exists("theme_hr")) source("R/figures/00_style.R")
PD <- readRDS(PANEL_DATA_RDS)$data

# Panels A and B: mean copy number of single replicates ----

#' Mean copy number of one replicate against time
#'
#' The replicate is recorded at 0 h and at each of the 12 transfers, so the
#' points mark the recorded rows and the line only joins them.
#'
#' @param df Data frame of one replicate with `time` (h) and `mean_copy`
#'   (mean copy number of the population), one row per recorded time.
#' @return A ggplot object.
mc_panel <- function(df) {
  ggplot(df, aes(time, mean_copy)) +
    geom_line(colour = HR_PAL$line, linewidth = 0.5) +
    geom_point(colour = HR_PAL$line, size = 0.8) +
    scale_y_continuous(limits = c(0, 9.5), breaks = c(0, 3, 6, 9)) +
    scale_x_continuous(breaks = c(0, 96, 192, 288)) +
    labs(x = "Time (h)", y = "Mean copy number") +
    theme_hr()
}
pA <- mc_panel(PD$Fig1A)
pB <- mc_panel(PD$Fig1B)

# Panel C: outcome frequencies ----

# Bars and error bars are in percent of the 4,000 replicates and are drawn
# from the stored summary, whose percentages and interval bounds are rounded
# to three decimals. The factor levels fix the order of the bars. An outcome
# that no replicate reached has a bar of no height, and only its error bar
# shows.
d1c <- PD$Fig1C
d1c$outcome <- factor(d1c$outcome, levels = c("Climb", "Low-copy", "Extinction"))
pC <- ggplot(d1c, aes(outcome, pct, fill = outcome)) +
  geom_col(width = 0.65) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                width = 0.18, linewidth = 0.3, colour = HR_PAL$ink) +
  scale_fill_manual(values = outcome_colors, guide = "none") +
  scale_y_continuous(limits = c(0, 100), expand = expansion(c(0, 0.03))) +
  labs(x = NULL, y = "Frequency (%)") +
  theme_hr() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

# Panel D: control without compensation at 288 h ----

# Only the ten non-compensated states are drawn. Bar heights are the mean
# share of each state in the total population at 288 h over the 100 replicates
# of the control trajectory batch, and the error bars reach one standard
# deviation either side of the mean. At 288 h every compensated state is at
# zero in every replicate of that batch.
#
# A bar on a logarithmic axis needs a finite lower end, and geom_col() draws
# bars from zero. Each bar is therefore a rectangle from floorY up to the
# mean, and the mean and both ends of the error bar are raised to floorY where
# they lie below it. With the shipped results this applies to R0 alone: its
# mean share, 6.2e-7, and its whole error bar lie below the floor, so both
# are drawn with no height on the bottom axis. Fill colour repeats the copy
# number, which the x axis already gives, so the fill legend is switched off.
d1d    <- subset(PD$Fig1D, compartment == "Non-compensated")
floorY <- 1e-6
pD <- ggplot(d1d, aes(x = copy, fill = factor(copy))) +
  geom_rect(aes(xmin = copy - 0.4, xmax = copy + 0.4,
                ymin = floorY, ymax = pmax(mean_freq, floorY)),
            colour = "black", linewidth = 0.2) +
  geom_errorbar(aes(ymin = pmax(mean_freq - sd_freq, floorY),
                    ymax = pmax(mean_freq + sd_freq, floorY)),
                width = 0.3, linewidth = 0.25, colour = HR_PAL$ink) +
  scale_copy_fill(guide = "none") +
  scale_x_continuous(breaks = 0:9, labels = paste0("R", 0:9)) +
  scale_y_log10(breaks = c(1e-6, 1e-4, 1e-2, 1e0),
                labels = expression(10^-6, 10^-4, 10^-2, 10^0),
                limits = c(floorY, 1.2),
                expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Copy number", y = "Frequency") +
  theme_hr()

# Assemble and save ----

# A and B fill the top row; C and D share the bottom row, with D 2.2 times as
# wide as C. patchwork adds the panel tags A to D.
bottom <- (pC | pD) + plot_layout(widths = c(1, 2.2))
fig1 <- (pA | pB) / bottom +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text())

save_fig(fig1, "Fig1", width_mm = WIDTH_DOUBLE, height_mm = 120)

# The figure object is the last value of the script, which hands it to the
# caller and displays it when the script is run line by line.
fig1
