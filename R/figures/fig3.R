# fig3.R
#
# Draws Figure 3, Scenario 2, in which the resistance of a compensated state
# rises with its copy number, and saves it as figures/Fig3.png.
#
# Panels and the panel data frames they are drawn from:
#   A  Fig3A  mean copy number against time for Scenario 1 and Scenario 2
#   B  Fig3B  Scenario 2 heatmap of each compensated state as a fraction of
#             the total population, 0 to 3,000 h
#   C  Fig3C  net growth rate of each compensated state under Scenario 2,
#             computed from the parameters
#
# Panels A and B are computed from the mean densities over the 100 replicates
# of the 3,000 h batches, not from single runs: panel A is the mean copy
# number of the replicate-mean densities, and each cell of panel B is the
# replicate-mean density of a compensated state divided by the replicate-mean
# total density. Panel C involves no simulation.
#
# The frames are read from source_data/panel_data.rds, which 01_panel_data.R
# writes; if that file is absent, source 01_panel_data.R first. The working
# directory must be the repository root. The panels carry no titles, because
# the figure legend of the paper describes them.
#
# Part of the code for "Unstable Heteroresistance as a Stepping Stone to Stable
# Resistance via Compensatory Mutations" (Berryhill and Andersson).
#
# Defines: PD, the panels pA, pB and pC, and the assembled figure fig3.

if (!exists("theme_hr")) source("R/figures/00_style.R")
PD <- readRDS(PANEL_DATA_RDS)$data

# Panel A: mean copy number over time, both scenarios ----

# Each line is the mean copy number of the mean densities of the 100
# replicates, that is of the replicates pooled, at each time of the 12 h
# output grid. It is not the average of the replicates' own mean copy
# numbers.
pA <- ggplot(PD$Fig3A, aes(t, mean_copy, colour = scenario)) +
  geom_line(linewidth = 0.55) +
  scale_colour_manual(values = scenario_colors, name = NULL) +
  scale_y_continuous(limits = c(0, 9.5), breaks = c(0, 3, 6, 9)) +
  labs(x = "Time (h)", y = "Mean copy number") +
  theme_hr() +
  theme(legend.position = "top")

# Panel B: Scenario 2 heatmap of the compensated states ----

# As in Fig 2D, each cell is the mean density of the compensated state with
# that copy number divided by the total mean density of all 20 states, so the
# ten cells at one time add up to the compensated fraction of the population.
# Unlike Fig 2D the colour scale is given no limits and spans the range of
# the values.
pB <- ggplot(PD$Fig3B, aes(t, copy, fill = freq)) +
  geom_raster() +
  scale_freq_fill() +
  scale_y_continuous(breaks = c(0, 3, 6, 9), expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Time (h)", y = "Copy number") +
  theme_hr()

# Panel C: net growth rate against copy number under Scenario 2 ----

# The rates are those of the ten compensated states at an antibiotic
# concentration of 16 micrograms per mL, the dose added at each transfer, as
# 01_panel_data.R computes them from the parameters. The dashed horizontal
# line marks zero net growth.
pC <- ggplot(PD$Fig3C, aes(copy, net_growth)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = HR_PAL$ref,
             linetype = "22") +
  geom_line(colour = HR_PAL$line, linewidth = 0.55) +
  geom_point(colour = HR_PAL$line, size = 1.3) +
  scale_x_continuous(breaks = 0:9) +
  labs(x = "Copy number", y = "Net per-capita growth (/h)") +
  theme_hr()

# Assemble and save ----

# The three panels sit in one row, and patchwork adds the panel tags A to C.
fig3 <- (pA | pB | pC) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text())

save_fig(fig3, "Fig3", width_mm = WIDTH_DOUBLE, height_mm = 62)

# The figure object is the last value of the script, which hands it to the
# caller and displays it when the script is run line by line.
fig3
