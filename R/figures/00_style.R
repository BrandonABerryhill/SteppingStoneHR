# 00_style.R
#
# Everything the figures share: the state space as the figure scripts use it,
# the font, the palettes, the ggplot theme, the figure widths, the output
# folders and the helper that saves a figure. 01_panel_data.R and every figure
# script source this file if it has not been sourced already, so a change made
# here reaches every panel.
#
# Working directory. This file takes the working directory at the moment it is
# sourced as the repository root, the folder that holds R/, results/ and the
# notebook, and builds the output paths from it. Sourcing it creates the
# folders figures/ and source_data/ there if they do not exist.
#
# Part of the code for "Unstable Heteroresistance as a Stepping Stone to Stable
# Resistance via Compensatory Mutations" (Berryhill and Andersson).
#
# Defines: STATES, COPY, state_info, HR_FONT, HR_TOL, HR_PAL, outcome_colors,
# scenario_colors, compartment_linetypes, HR_COPY, scale_copy_colour(),
# scale_copy_fill(), scale_freq_fill(), BASE_SIZE, theme_hr(), WIDTH_SINGLE,
# WIDTH_1_5, WIDTH_DOUBLE, PROJ, FIG_DIR, DATA_DIR, PANEL_DATA_RDS and
# save_fig().

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

# State space ----

# The twenty states in the order the engines use: R0 to R9 without a
# compensatory mutation, then R0C to R9C with one. COPY is the copy number of
# each state in that order.
STATES <- c(paste0("R", 0:9), paste0("R", 0:9, "C"))
COPY   <- c(0:9, 0:9)
state_info <- data.frame(
  state       = STATES,
  copy        = COPY,
  compartment = rep(c("Non-compensated", "Compensated"), each = 10),
  stringsAsFactors = FALSE
)

# Fonts ----

# All text is set in Arial, axis text included. The PNG device finds it
# through systemfonts. If systemfonts does not list Arial, HR_FONT is set to
# "", which selects the device's default sans-serif font, and one message
# says so.
HR_FONT <- "Arial"
if (requireNamespace("systemfonts", quietly = TRUE)) {
  .fams <- systemfonts::system_fonts()$family
  if (!(HR_FONT %in% .fams)) {
    message("00_style.R: '", HR_FONT,
            "' not found by systemfonts; falling back to default sans. ",
            "Install Arial or change HR_FONT.")
    HR_FONT <- ""
  }
}

# Palette ----

# Categories are coloured with Paul Tol's muted qualitative scheme.
# Non-compensated and compensated states are told apart by line type, not by
# colour.
HR_TOL <- c(
  indigo = "#332288", cyan = "#88CCEE", teal   = "#44AA99",
  green  = "#117733", olive = "#999933", sand   = "#DDCC77",
  rose   = "#CC6677", wine  = "#882255", purple = "#AA4499",
  grey   = "#DDDDDD"
)

# Colours that carry no category: ink for text and axes, line for a single
# data series, ref for faint reference lines such as a zero line.
HR_PAL <- list(
  ink  = "#1a1a1a",
  line = "#332288",
  ref  = "#999999"
)

# The three outcomes (Fig 1C, Supplementary Fig 5).
outcome_colors <- c(
  "Climb"      = unname(HR_TOL["indigo"]),
  "Low-copy"   = unname(HR_TOL["rose"]),
  "Extinction" = unname(HR_TOL["sand"])
)

# The two scenarios (Fig 3A).
scenario_colors <- c(
  "Scenario 1" = unname(HR_TOL["rose"]),
  "Scenario 2" = unname(HR_TOL["indigo"])
)

# Line type by compartment (Fig 2, Supplementary Fig 3): solid for
# non-compensated states and short dashes for compensated ones.
compartment_linetypes <- c(
  "Non-compensated" = "solid",
  "Compensated"     = "22"
)

# One colour per copy number, 0 to 9: charcoal for copy 0 and the nine colours
# of the muted scheme for copies 1 to 9. The scales below are discrete, so the
# colour or fill aesthetic has to be mapped to factor(copy).
HR_COPY <- c(
  "0" = "#444444", "1" = "#332288", "2" = "#88CCEE", "3" = "#44AA99",
  "4" = "#117733", "5" = "#999933", "6" = "#DDCC77", "7" = "#CC6677",
  "8" = "#882255", "9" = "#AA4499"
)
scale_copy_colour <- function(...) {
  scale_colour_manual(values = HR_COPY, name = "Copy\nnumber", ...)
}
scale_copy_fill <- function(...) {
  scale_fill_manual(values = HR_COPY, name = "Copy\nnumber", ...)
}
# Continuous fill for the heatmap panels: option "C" of the viridis scales.
scale_freq_fill <- function(...) {
  scale_fill_viridis_c(name = "Frequency", option = "C", ...)
}

# Type sizes ----

# Point size of axis text. The theme sets every other text size relative to
# it.
BASE_SIZE <- 8

# Theme ----

#' Theme shared by every panel
#'
#' @param base_size Numeric scalar, point size of axis text. Axis titles and
#'   plot titles are one point larger, panel tags four points larger and
#'   legend text one point smaller.
#' @param base_family Character scalar, font family; "" selects the device's
#'   default sans-serif font.
#' @return A complete ggplot2 theme, built on `theme_classic()`.
theme_hr <- function(base_size = BASE_SIZE, base_family = HR_FONT) {
  theme_classic(base_size = base_size, base_family = base_family) %+replace%
    theme(
      axis.line   = element_line(linewidth = 0.3, colour = HR_PAL$ink),
      axis.ticks  = element_line(linewidth = 0.3, colour = HR_PAL$ink),
      axis.text   = element_text(size = base_size, colour = HR_PAL$ink),
      axis.title  = element_text(size = base_size + 1, colour = HR_PAL$ink),
      axis.title.y = element_text(angle = 90, margin = margin(r = 2)),
      plot.title  = element_text(size = base_size + 1, hjust = 0,
                                 margin = margin(b = 3), colour = HR_PAL$ink),
      plot.tag    = element_text(size = base_size + 4, face = "bold",
                                 colour = HR_PAL$ink),
      legend.key.size = unit(3.5, "mm"),
      legend.text  = element_text(size = base_size - 1, colour = HR_PAL$ink),
      legend.title = element_text(size = base_size,     colour = HR_PAL$ink),
      legend.position = "right",
      legend.background = element_blank(),
      legend.key = element_blank(),
      strip.background = element_blank(),
      strip.text  = element_text(size = base_size, face = "bold",
                                 colour = HR_PAL$ink),
      plot.margin = margin(2, 2, 2, 2, "mm")
    )
}

# Text drawn by a geom does not take its font from the theme, so the default
# family of the text and label geoms is set as well.
if (nzchar(HR_FONT)) {
  update_geom_defaults("text",  list(family = HR_FONT))
  update_geom_defaults("label", list(family = HR_FONT))
}

# Figure dimensions ----

# Figure widths in millimetres, passed to save_fig(): single column, one and a
# half columns, and double column (full page width).
WIDTH_SINGLE <- 89
WIDTH_1_5    <- 120
WIDTH_DOUBLE <- 183

# Output folders ----

# PROJ is the working directory when this file is sourced, which has to be the
# repository root. Figures are written to FIG_DIR, and 01_panel_data.R writes
# the panel data to PANEL_DATA_RDS, from which the figure scripts read it.
PROJ      <- getwd()
FIG_DIR   <- file.path(PROJ, "figures")
DATA_DIR  <- file.path(PROJ, "source_data")
for (d in c(FIG_DIR, DATA_DIR))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

PANEL_DATA_RDS <- file.path(DATA_DIR, "panel_data.rds")

# Save helper ----

#' Save a figure as PDF and PNG
#'
#' Writes `<name>.pdf` and `<name>.png` to FIG_DIR with the same dimensions.
#' The PNG is the reference output. The PDF is an extra: if the PDF device
#' cannot be opened, or the export stops with an error, the PDF is skipped
#' with a message and the PNG is still written.
#'
#' @param plot A ggplot or patchwork object.
#' @param name Character scalar, file name without extension.
#' @param width_mm,height_mm Numeric scalars, figure size in millimetres.
#' @param dpi Numeric scalar, resolution of the PNG in dots per inch.
#' @return `plot`, invisibly.
save_fig <- function(plot, name, width_mm, height_mm, dpi = 600) {
  pdf_dev <- if (requireNamespace("grDevices", quietly = TRUE))
    grDevices::cairo_pdf else "pdf"
  pdf_file <- file.path(FIG_DIR, paste0(name, ".pdf"))
  pdf_skipped <- function(cond) {
    message("  PDF export of ", name, " was skipped: ",
            strsplit(conditionMessage(cond), "\n")[[1]][1])
    FALSE
  }
  # A PDF device that cannot be opened signals a warning, not an error, and
  # ggsave() would then draw the figure on whatever device the session opens
  # by default. The device is therefore opened and closed once on a temporary
  # file first, and a warning or an error at that point skips the PDF.
  pdf_written <- tryCatch({
    if (is.function(pdf_dev)) {
      current <- grDevices::dev.cur()
      probe   <- tempfile(fileext = ".pdf")
      pdf_dev(filename = probe)
      grDevices::dev.off()
      unlink(probe)
      if (current > 1) grDevices::dev.set(current)
    }
    TRUE
  }, warning = pdf_skipped, error = pdf_skipped)
  if (pdf_written) {
    pdf_written <- tryCatch({
      ggsave(pdf_file, plot = plot,
             width = width_mm, height = height_mm, units = "mm",
             device = pdf_dev)
      TRUE
    }, error = function(cond) {
      # Whatever the stopped export has written is deleted, so that no PDF is
      # kept from an export that did not finish.
      unlink(pdf_file)
      pdf_skipped(cond)
    })
  }
  png_dev <- if (requireNamespace("ragg", quietly = TRUE))
    ragg::agg_png else "png"
  ggsave(file.path(FIG_DIR, paste0(name, ".png")), plot = plot,
         width = width_mm, height = height_mm, units = "mm",
         dpi = dpi, device = png_dev)
  message("  wrote ", if (pdf_written) paste0(name, ".pdf and "),
          name, ".png  (", width_mm, " x ", height_mm, " mm)")
  invisible(plot)
}
