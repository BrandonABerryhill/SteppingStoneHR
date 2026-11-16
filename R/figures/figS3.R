# figS3.R
#
# Draws Supplementary Figure 3, the dynamics of the control without
# compensation (compensation rate set to zero) to 3,000 h, and saves it as
# figures/FigS3.png.
#
# One panel, drawn from the panel data frame FigS3: the mean density of each
# of the 20 states over the 100 replicates of the 3,000 h control batch, at
# each time of the 12 h output grid. Every curve is that mean, not a single
# run.
#
# The frame is read from source_data/panel_data.rds, which 01_panel_data.R
# writes; if that file is absent, source 01_panel_data.R first. The working
# directory must be the repository root.
#
# Part of the code for "Unstable Heteroresistance as a Stepping Stone to Stable
# Resistance via Compensatory Mutations" (Berryhill and Andersson).
#
# Defines: PD, d and the figure figS3.

if (!exists("theme_hr")) source("R/figures/00_style.R")
PD <- readRDS(PANEL_DATA_RDS)$data

# The y axis is logarithmic and cannot show zero, so mean densities below
# 1 cell per mL are drawn at 1 cell per mL. The ten compensated states are
# at zero at every time in this batch, and their dashed lines lie on top of
# one another along that floor. With the shipped results three values of
# non-compensated states are raised as well: R8 and R9 at 0 h, whose
# founder densities are below 1, and R0 at 12 h.
d <- PD$FigS3
d$dens <- pmax(d$density, 1)

# All 20 states are drawn, one line each: colour gives the copy number and
# line type the compartment, solid for non-compensated and dashed for
# compensated states.
figS3 <- ggplot(d, aes(t, dens, colour = factor(copy), group = state,
                       linetype = compartment)) +
  geom_line(linewidth = 0.4) +
  scale_y_log10() +
  scale_copy_colour() +
  scale_linetype_manual(values = compartment_linetypes, name = NULL) +
  labs(x = "Time (h)", y = "Density (cells/mL)") +
  theme_hr()

save_fig(figS3, "FigS3", width_mm = WIDTH_1_5, height_mm = 78)

# The figure object is the last value of the script, which hands it to the
# caller and displays it when the script is run line by line.
figS3
