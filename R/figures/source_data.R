# source_data.R
#
# Writes the source data workbook, source_data/Source_Data.xlsx: an Index tab
# followed by one tab for each of the 24 data tabs listed in SOURCE_DATA_TABS,
# in manuscript order. Each data tab holds the display item's title, the date
# the code behind it was run, the question the panel answers and its data:
# the plotted values and, for the tabs to which source_data_tab() gives one
# (Fig 1C, Fig 1D, Fig 2A to 2C, Fig 3A and Supp Fig 3), a table of
# per-replicate values from the batch behind the panel. A per-replicate table
# with more than MAX_BLOCK_ROWS rows is left out, and the Index names the
# files it is computed from; with the shipped results that is the table of
# Fig 3A. Fig 2D and Fig 3B have no per-replicate table. The tables come from
# source_data_tab() of 01_panel_data.R, the same frames the figures are drawn
# from, and are written as they are returned: nothing is filtered, averaged,
# rounded or converted here. This file only lays the tables out, gives each
# column a header with its unit, and looks up the run dates.
#
# Use. The working directory must be the repository root. Sourcing this file
# defines the functions below and then writes the workbook once with
# write_source_data(). It needs the panel data layer and, for the
# per-replicate table of Fig 1C, the two supplementary tables and the note on
# the outcomes in the Index, 00_config.R, 01_parameters.R and 04_classify.R of
# R/simulation/. If source_data_tab() is not defined yet, this file sources
# R/figures/01_panel_data.R first; the simulation modules have to be sourced
# by the caller. The run dates are read from SOURCE_DATA_RUN_DATES if a script
# has assigned it before sourcing this file (a data frame, or the path of a
# CSV file, with the columns batch and run_date), and otherwise from
# results/run_dates.csv. A re-simulation passes the dates of its own run this
# way. The number of seeds per parameter set of the sweep, shown on the Index,
# is taken from SWEEP_SEEDS if that is defined.
#
# Part of the code for "Unstable Heteroresistance as a Stepping Stone to Stable
# Resistance via Compensatory Mutations" (Berryhill and Andersson).
#
# Defines: SOURCE_DATA_XLSX, MAX_BLOCK_ROWS, TAB_QUESTIONS, TAB_BATCHES,
# COLUMN_HEADERS, display_item(), column_headers(), read_run_dates(),
# tab_run_date(), index_rows(), write_data_tab(), write_index_tab() and
# write_source_data().

if (!exists("source_data_tab")) source("R/figures/01_panel_data.R")

suppressPackageStartupMessages(library(openxlsx))

# Settings ----

# Output file, in the folder that 00_style.R creates for the source data.
SOURCE_DATA_XLSX <- file.path(DATA_DIR, "Source_Data.xlsx")

# A per-replicate table with more rows than this is left out of the workbook,
# and the Index names the files it is computed from.
MAX_BLOCK_ROWS <- 50000

# Questions ----

# The question each tab answers, written in cell B3. For Fig 1A and 1B the
# seed of the replicate shown is part of the question line; these are the two
# seeds 01_panel_data.R takes from the Scenario 1 trajectory batch.
TAB_QUESTIONS <- c(
  "Fig 1A" = "How does mean copy number change over a treatment course in a population that climbs to high copy number? The replicate shown has seed 2.",
  "Fig 1B" = "How does mean copy number change in a population that survives at low copy number? The replicate shown has seed 7.",
  "Fig 1C" = "How often does each outcome occur across 4,000 replicate populations under Scenario 1?",
  "Fig 1D" = "Without compensation, what is the distribution of states after 288 h?",
  "Fig 2A" = "How do the states change over the first 72 h under Scenario 1?",
  "Fig 2B" = "How do the states change over the first 1,000 h?",
  "Fig 2C" = "How do the states change over 3,000 h?",
  "Fig 2D" = "How does the share of each compensated state change over 3,000 h under Scenario 1?",
  "Fig 3A" = "Does mean copy number descend over 3,000 h under Scenario 1 and under Scenario 2?",
  "Fig 3B" = "How does the share of each compensated state change over 3,000 h under Scenario 2?",
  "Fig 3C" = "How does net growth at the treatment dose depend on copy number in compensated cells under Scenario 2?",
  "Fig 4A" = "How sensitive is the frequency of the climb outcome to each of six parameters?",
  "Fig 4B" = "How sensitive is the frequency of the low-copy outcome to each parameter?",
  "Fig 4C" = "How sensitive is the frequency of extinction to each parameter?",
  "Supp Fig 3" = "What are the dynamics over 3,000 h when compensation is absent?",
  "Supp Fig 4" = "When does the zero-copy compensated state overtake the single-copy state?",
  "Supp Fig 5A" = "How do the outcome frequencies vary with antibiotic dose across the 300 sampled parameter sets?",
  "Supp Fig 5B" = "How do the outcome frequencies vary with compensation rate across the 300 sampled parameter sets?",
  "Supp Fig 5C" = "How do the outcome frequencies vary with per-copy cost across the 300 sampled parameter sets?",
  "Supp Fig 5D" = "How do the outcome frequencies vary with recombination rate across the 300 sampled parameter sets?",
  "Supp Fig 5E" = "How do the outcome frequencies vary with serial transfer rate across the 300 sampled parameter sets?",
  "Supp Fig 5F" = "How do the outcome frequencies vary with recovered fraction across the 300 sampled parameter sets?",
  "Supp Table 1" = "What are the states of the model?",
  "Supp Table 2" = "What parameter values does the model use?"
)

# Run dates ----

# The entry of the run-dates table that dates each tab. A tab drawn from two
# batches (Fig 3A) lists both, named by the scenario each one holds. Fig 3C and
# the two supplementary tables have entries of their own: Fig 3C is computed
# from parameters, and the two tables take the date of the Scenario 1
# frequency batch whose parameters they report.
TAB_BATCHES <- list(
  "Fig 1A" = "S1_traj300", "Fig 1B" = "S1_traj300",
  "Fig 1C" = "S1_freq_cpp", "Fig 1D" = "omega0_traj300",
  "Fig 2A" = "S1_long", "Fig 2B" = "S1_long", "Fig 2C" = "S1_long",
  "Fig 2D" = "S1_long",
  "Fig 3A" = c("Scenario 1" = "S1_long", "Scenario 2" = "S2_long"),
  "Fig 3B" = "S2_long", "Fig 3C" = "Fig 3C",
  "Fig 4A" = "sweep_S1", "Fig 4B" = "sweep_S1", "Fig 4C" = "sweep_S1",
  "Supp Fig 3" = "omega0_long", "Supp Fig 4" = "long_run",
  "Supp Fig 5A" = "sweep_S1", "Supp Fig 5B" = "sweep_S1",
  "Supp Fig 5C" = "sweep_S1", "Supp Fig 5D" = "sweep_S1",
  "Supp Fig 5E" = "sweep_S1", "Supp Fig 5F" = "sweep_S1",
  "Supp Table 1" = "Supp Table 1", "Supp Table 2" = "Supp Table 2"
)

#' Read the run-dates table
#'
#' @param run_dates A data frame, or the path of a CSV file, with at least the
#'   columns `batch` and `run_date` (ISO date, YYYY-MM-DD), one row per batch
#'   and per computed tab.
#' @return A named character vector of run dates, named by batch. Stops if an
#'   entry is duplicated or a date is not an ISO date.
read_run_dates <- function(run_dates) {
  if (is.character(run_dates))
    run_dates <- read.csv(run_dates, colClasses = "character",
                          stringsAsFactors = FALSE)
  if (!all(c("batch", "run_date") %in% names(run_dates)))
    stop("The run-dates table needs the columns batch and run_date.",
         call. = FALSE)
  if (anyDuplicated(run_dates$batch))
    stop("The run-dates table lists a batch more than once.", call. = FALSE)
  dates <- as.character(run_dates$run_date)
  if (!all(grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", dates)))
    stop("Every run date must be an ISO date, YYYY-MM-DD.", call. = FALSE)
  setNames(dates, run_dates$batch)
}

#' Run date of one tab, as written in cell B2
#'
#' @param tab Character scalar, one of `SOURCE_DATA_TABS`.
#' @param dates Named character vector from `read_run_dates()`.
#' @return Character scalar: the ISO date of the tab's batch, or, for a tab
#'   drawn from two batches, each date after the name of its scenario
#'   ("Scenario 1: YYYY-MM-DD; Scenario 2: YYYY-MM-DD"). Stops if a batch has
#'   no entry, so no tab is ever written without a date from the table.
tab_run_date <- function(tab, dates) {
  batches <- TAB_BATCHES[[tab]]
  absent  <- setdiff(batches, names(dates))
  if (length(absent) > 0)
    stop("No run date for ", paste(absent, collapse = ", "), " (tab ", tab,
         ") in the run-dates table.", call. = FALSE)
  if (length(batches) == 1) return(unname(dates[batches]))
  paste(paste0(names(batches), ": ", dates[batches]), collapse = "; ")
}

# Column headers ----

# Header of every column that source_data_tab() returns, with its unit. A row
# with `tabs` empty is the header of that column in every tab and block; a
# row that names tabs (separated by "|") and a block ("plotted" or
# "per_replicate") replaces it there. Units: times are in hours and densities
# in cells per mL (R/simulation/01_parameters.R and the documentation of the
# functions in 01_panel_data.R that build each table); the units of the
# swept parameters are those of the arguments of make_condition(); a PAWN
# index is a Kolmogorov-Smirnov distance and has none.
COLUMN_HEADERS <- local({
  h <- function(column, header, tabs = "", block = "")
    data.frame(column = column, header = header, tabs = tabs, block = block,
               stringsAsFactors = FALSE)
  mean_panels <- "Fig 2A|Fig 2B|Fig 2C|Fig 2D|Fig 3A|Fig 3B|Supp Fig 3"
  st <- c(paste0("R", 0:9), paste0("R", 0:9, "C"))
  rbind(
    # Time.
    h("time", "Time (h)"),
    h("time", "Time, rounded to the hour (h)", mean_panels, "plotted"),
    h("time", "Time as recorded (h)", mean_panels, "per_replicate"),
    # States, totals and copy number.
    h(st, paste0(st, " (cells per mL)")),
    h(st, paste0(st, " at 288 h (fraction of population)"), "Fig 1D",
      "per_replicate"),
    h("total", "Total density (cells per mL)"),
    h("mean_copy", "Mean copy number"),
    h("mean_copy", "Mean copy number of the replicate-mean densities",
      "Fig 3A", "plotted"),
    h("mean_copy", "Mean copy number of the replicate", "Fig 3A",
      "per_replicate"),
    h("frac_compensated", "Compensated cells (fraction of population)"),
    h("state", "State"),
    h("copy", "Copy number"),
    h("copy", "Copy number of the compensated state", "Fig 2D|Fig 3B|Fig 3C",
      "plotted"),
    h("compartment", "Compartment"),
    h("density", "Mean density over the replicates (cells per mL)"),
    h("freq", paste("Frequency of the compensated state (fraction of the",
                    "population, from the replicate-mean densities)")),
    h("scenario", "Scenario"),
    # Outcomes.
    h("outcome", "Outcome"),
    h("n", "Number of replicates"),
    h("pct", paste("Percent of replicates (%, as stored by the run, which",
                   "rounded it to three decimals)")),
    h("ci_low", paste("Lower 95% Clopper-Pearson bound (%, as stored by the",
                      "run, which rounded it to three decimals)")),
    h("ci_high", paste("Upper 95% Clopper-Pearson bound (%, as stored by the",
                       "run, which rounded it to three decimals)")),
    h("seed", "Seed"),
    h("peak_mACN", "Peak mean copy number over the transfer rows"),
    h("final_N", "Total density at 288 h (cells per mL)"),
    h("mean_freq", "Mean frequency at 288 h (fraction of population)"),
    h("sd_freq",
      "Standard deviation of the frequency at 288 h (fraction of population)"),
    # Net growth.
    h("net_growth", "Net per-capita growth rate (per h)"),
    # PAWN indices.
    h("parameter", "Parameter"),
    h("label", "Parameter label"),
    h("median", "PAWN index, median over the quantile bins (dimensionless)"),
    h("max", "PAWN index, maximum over the quantile bins (dimensionless)"),
    # The sweep.
    h("set", "Parameter set"),
    h("amax", "Dose, AMAX (µg per mL)"),
    h("omega", "Compensation rate, ω (per cell per hour)"),
    h("cost_per_copy", "Per-copy cost (per cell per hour)"),
    h("mu", "Recombination rate, μ (per cell per hour)"),
    h("dil", "Serial transfer rate (dimensionless)"),
    h("comp_frac", "Recovered fraction (dimensionless)"),
    h("f_climb", "Frequency of climb (fraction of the set's replicates)"),
    h("f_lowcopy", "Frequency of low-copy (fraction of the set's replicates)"),
    h("f_extinct",
      "Frequency of extinction (fraction of the set's replicates)"),
    # Supplementary tables.
    h("dimensions", "Dimensions"),
    h("definition", "Definition"),
    h("value_S1", "Value in Scenario 1 (in the unit under Dimensions)"),
    h("value_S2", "Value in Scenario 2 (in the unit under Dimensions)"),
    h("source", paste("Source (reference number in the manuscript, or This",
                      "Study, or Assumed)"))
  )
})

#' Headers of the columns of one table
#'
#' @param columns Character vector, column names of a table returned by
#'   `source_data_tab()`.
#' @param tab Character scalar, the tab.
#' @param block Character scalar, "plotted" or "per_replicate".
#' @return Character vector of headers, one per column. Stops if a column has
#'   no header in `COLUMN_HEADERS`.
column_headers <- function(columns, tab, block) {
  vapply(columns, function(column) {
    rows <- COLUMN_HEADERS[COLUMN_HEADERS$column == column, ]
    here <- rows[rows$block == block &
                 vapply(strsplit(rows$tabs, "|", fixed = TRUE),
                        function(t) tab %in% t, logical(1)), ]
    if (nrow(here) == 0) here <- rows[rows$tabs == "", ]
    if (nrow(here) != 1)
      stop("No header for column ", column, " of tab ", tab, ".",
           call. = FALSE)
    here$header
  }, character(1), USE.NAMES = FALSE)
}

# Index ----

#' Title of the display item of a tab, as written in cell A1
#'
#' @param tab Character scalar, for example "Fig 1A" or "Supp Table 2".
#' @return Character scalar, for example "Figure 1A" or "Supplementary Table
#'   2".
display_item <- function(tab) {
  tab <- sub("^Supp Fig ", "Supplementary Figure ", tab)
  tab <- sub("^Supp Table ", "Supplementary Table ", tab)
  sub("^Fig ", "Figure ", tab)
}

#' One line of the Index for every data tab
#'
#' @param tables Named list, one element per tab, each the list returned by
#'   `source_data_tab()`.
#' @param dates Named character vector from `read_run_dates()`.
#' @param results Character scalar, the results folder as a path relative to
#'   the repository root, used to name files.
#' @param sweep_seeds Integer scalar, seeds per parameter set of the sweep, or
#'   `NULL` if unknown.
#' @return A list with one element per tab, each a list of the eight cells of
#'   its Index line: tab, display item, question, run date, source batch,
#'   engine, number of replicates and number of data rows. The last two are
#'   numbers where they are a single count and text otherwise.
index_rows <- function(tables, dates, results, sweep_seeds) {
  rp <- function(...) file.path(results, ...)
  sims <- function(scenario, tag) rp(scenario, tag, "sims.rds")
  seeds_in <- function(d, scenario = NULL) {
    if (!is.null(scenario)) d <- d[d$scenario == scenario, ]
    length(unique(d$seed))
  }
  n_sets <- nrow(tables[["Supp Fig 5A"]]$plotted)
  sweep_reps <- if (is.null(sweep_seeds))
    paste(n_sets, "parameter sets") else
    paste(n_sets, "parameter sets of", sweep_seeds, "replicates")
  s1_long <- seeds_in(tables[["Fig 2C"]]$per_replicate)
  pr3A    <- tables[["Fig 3A"]]$per_replicate

  compiled <- "Compiled (hrsim)"
  from_par <- "Computed from parameters"
  info <- list(
    "Fig 1A" = list(paste0("S1_traj300 (", sims("S1", "S1_traj300"),
                           "), seed 2"), compiled, 1),
    "Fig 1B" = list(paste0("S1_traj300 (", sims("S1", "S1_traj300"),
                           "), seed 7"), compiled, 1),
    "Fig 1C" = list(paste0("S1_freq_cpp: stored frequencies in ",
                           rp("frequency_summary.rds"),
                           "; per-replicate values classified from ",
                           sims("S1", "S1_freq_cpp")), compiled,
                    nrow(tables[["Fig 1C"]]$per_replicate)),
    "Fig 1D" = list(paste0("omega0_traj300 (",
                           sims("omega0", "omega0_traj300"), ")"), compiled,
                    seeds_in(tables[["Fig 1D"]]$per_replicate)),
    "Fig 2A" = list(paste0("S1_long (", sims("S1", "S1_long"), ")"),
                    compiled, seeds_in(tables[["Fig 2A"]]$per_replicate)),
    "Fig 2B" = list(paste0("S1_long (", sims("S1", "S1_long"), ")"),
                    compiled, seeds_in(tables[["Fig 2B"]]$per_replicate)),
    "Fig 2C" = list(paste0("S1_long (", sims("S1", "S1_long"), ")"),
                    compiled, s1_long),
    "Fig 2D" = list(paste0("S1_long (", sims("S1", "S1_long"), ")"),
                    compiled, s1_long),
    "Fig 3A" = list(paste0("S1_long (", sims("S1", "S1_long"),
                           ") for Scenario 1 and S2_long (",
                           sims("S2", "S2_long"), ") for Scenario 2"),
                    compiled,
                    paste0(seeds_in(pr3A, "Scenario 1"), " (Scenario 1) and ",
                           seeds_in(pr3A, "Scenario 2"), " (Scenario 2)")),
    "Fig 3B" = list(paste0("S2_long (", sims("S2", "S2_long"), ")"),
                    compiled, seeds_in(pr3A, "Scenario 2")),
    "Fig 3C" = list(paste("None: computed from the Scenario 2 parameters of",
                          "make_condition() at an antibiotic concentration",
                          "of 16 µg per mL"),
                    from_par, "None (no simulation)"),
    "Fig 4A" = list(paste0("sweep_S1 (", rp("sensitivity", "sweep_S1",
                                            "pawn.rds"), ")"),
                    compiled, sweep_reps),
    "Fig 4B" = list(paste0("sweep_S1 (", rp("sensitivity", "sweep_S1",
                                            "pawn.rds"), ")"),
                    compiled, sweep_reps),
    "Fig 4C" = list(paste0("sweep_S1 (", rp("sensitivity", "sweep_S1",
                                            "pawn.rds"), ")"),
                    compiled, sweep_reps),
    "Supp Fig 3" = list(paste0("omega0_long (",
                               sims("omega0", "omega0_long"), ")"),
                        compiled,
                        seeds_in(tables[["Supp Fig 3"]]$per_replicate)),
    "Supp Fig 4" = list(paste0("long_run (",
                               "results/long_run/long_run_percycle.csv, ",
                               "reduced from ",
                               "results/long_run/long_run_raw.csv)"),
                        paste("None of these: the data come from an old",
                              "iteration of the model which is not",
                              "regenerated by this code"), 1),
    "Supp Table 1" = list(paste("The states of the model, checked against",
                                "the states vector of",
                                "R/simulation/01_parameters.R; dated by",
                                "S1_freq_cpp, the Scenario 1 frequency",
                                "batch whose parameters it reports"),
                          "None (no simulation)", "None (no simulation)"),
    "Supp Table 2" = list(paste("make_condition() for S1 and S2",
                                "(R/simulation/01_parameters.R); dated by",
                                "S1_freq_cpp, the Scenario 1 frequency",
                                "batch whose parameters it reports"),
                          from_par, "None (no simulation)")
  )
  for (tab in sprintf("Supp Fig 5%s", LETTERS[1:6]))
    info[[tab]] <- list(paste0("sweep_S1 (", rp("sensitivity", "sweep_S1",
                                                "sweep.csv"), ")"),
                        compiled, sweep_reps)

  # Files a left-out per-replicate table is computed from.
  left_out_files <- list(
    "Fig 1C" = sims("S1", "S1_freq_cpp"),
    "Fig 1D" = sims("omega0", "omega0_traj300"),
    "Fig 2A" = sims("S1", "S1_long"), "Fig 2B" = sims("S1", "S1_long"),
    "Fig 2C" = sims("S1", "S1_long"),
    "Fig 3A" = c(sims("S1", "S1_long"), sims("S2", "S2_long")),
    "Supp Fig 3" = sims("omega0", "omega0_long")
  )

  rows <- lapply(SOURCE_DATA_TABS, function(tab) {
    x  <- tables[[tab]]
    np <- nrow(x$plotted)
    nr <- if (is.null(x$per_replicate)) 0L else nrow(x$per_replicate)
    n_rows <- if (nr == 0) {
      np
    } else if (nr <= MAX_BLOCK_ROWS) {
      paste0(np, " plotted values; ", nr, " per-replicate values")
    } else {
      paste0(np, " plotted values; the ", nr, " per-replicate values are ",
             "more than ", MAX_BLOCK_ROWS, " rows and are not in this ",
             "workbook: they are computed from ",
             paste(left_out_files[[tab]], collapse = " and "),
             " by the function source_data_tab() of ",
             "R/figures/01_panel_data.R")
    }
    if (tab == "Supp Fig 4")
      n_rows <- paste0(np, " plotted values, one per 24 h cycle; the raw ",
                       "hourly series is not in this workbook and is ",
                       "results/long_run/long_run_raw.csv")
    list(tab, display_item(tab), TAB_QUESTIONS[[tab]],
         tab_run_date(tab, dates), info[[tab]][[1]], info[[tab]][[2]],
         info[[tab]][[3]], n_rows)
  })
  rows
}

# Writing ----

#' Write one data tab
#'
#' @param wb An openxlsx Workbook with a worksheet named `tab`.
#' @param tab Character scalar, the tab.
#' @param x List returned by `source_data_tab(tab)`.
#' @param date Character scalar, the run date from `tab_run_date()`.
#' @param bold An openxlsx style that sets bold text and nothing else.
#' @return Logical scalar, `TRUE` if a per-replicate table was left out for
#'   length, invisibly.
write_data_tab <- function(wb, tab, x, date, bold) {
  put <- function(value, row, col = 1)
    writeData(wb, tab, value, startCol = col, startRow = row,
              colNames = FALSE)
  block <- function(d, kind, header_row) {
    names(d) <- column_headers(names(d), tab, kind)
    writeData(wb, tab, d, startRow = header_row, colNames = TRUE,
              keepNA = FALSE)
    addStyle(wb, tab, bold, rows = header_row, cols = seq_len(ncol(d)),
             stack = TRUE)
  }
  put(display_item(tab), 1)
  addStyle(wb, tab, bold, rows = 1, cols = 1, stack = TRUE)
  put("Date code was run", 2); put(date, 2, 2)
  put("Question", 3); put(TAB_QUESTIONS[[tab]], 3, 2)
  put("Data", 5)
  addStyle(wb, tab, bold, rows = 5, cols = 1, stack = TRUE)
  put("Plotted values", 6)
  block(x$plotted, "plotted", 7)
  freezePane(wb, tab, firstActiveRow = 8)

  pr <- x$per_replicate
  left_out <- !is.null(pr) && nrow(pr) > MAX_BLOCK_ROWS
  if (!is.null(pr) && !left_out) {
    # One blank row below the plotted values, then the label and the table.
    label_row <- 7 + nrow(x$plotted) + 2
    put("Per-replicate values", label_row)
    block(pr, "per_replicate", label_row + 1)
  }
  invisible(left_out)
}

#' Write the Index tab
#'
#' @param wb An openxlsx Workbook with a worksheet named "Index".
#' @param idx List from `index_rows()`.
#' @param bold An openxlsx style that sets bold text and nothing else.
#' @return `NULL`, invisibly.
write_index_tab <- function(wb, idx, bold) {
  put <- function(value, row, col = 1)
    writeData(wb, "Index", value, startCol = col, startRow = row,
              colNames = FALSE)
  put(paste("Unstable Heteroresistance as a Stepping Stone to Stable",
            "Resistance via Compensatory Mutations"), 1)
  addStyle(wb, "Index", bold, rows = 1, cols = 1, stack = TRUE)
  put("Brandon A. Berryhill and Dan I. Andersson", 2)

  header <- c("Tab", "Display item", "Question", "Date code was run",
              "Source batch", "Engine", "Number of replicates",
              "Number of data rows")
  put(t(header), 4)
  addStyle(wb, "Index", bold, rows = 4, cols = seq_along(header),
           stack = TRUE)
  for (i in seq_along(idx))
    for (j in seq_along(header)) put(idx[[i]][[j]], 4 + i, j)

  # Notes. The outcome rule is read from the classification settings of
  # 00_config.R and the transfer interval from make_condition(), so the
  # sentence states what the classifier does.
  transfer <- make_condition("S1")$p$transfer
  last_row <- floor(CLASSIFY_HORIZON / transfer) * transfer
  unit     <- if (EXTINCT_BELOW == 1) "cell per mL" else "cells per mL"
  notes <- c(
    paste0("Outcomes. Each replicate of a ", CLASSIFY_HORIZON, " h run is ",
           "classified on the rows recorded at its serial transfers, every ",
           transfer, " h up to ", last_row, " h: extinction if the total ",
           "density on the last of these rows is below ", EXTINCT_BELOW, " ",
           unit, "; otherwise climb if the mean copy number on at least one ",
           "of these rows is ", CLIMB_THRESHOLD, " or more; otherwise ",
           "low-copy."),
    "Units. Times are in hours and densities in cells per mL.",
    paste("Supplementary Figure 5. The smoothed lines are drawn by the",
          "plotting code from the points given in tabs Supp Fig 5A to",
          "Supp Fig 5F."),
    paste("Supplementary Figure 4. The data come from an old iteration of",
          "the model which is not regenerated by this code.")
  )
  first <- 4 + length(idx) + 2
  put("Notes", first)
  addStyle(wb, "Index", bold, rows = first, cols = 1, stack = TRUE)
  for (k in seq_along(notes)) put(notes[k], first + k)
  setColWidths(wb, "Index", cols = 1:8,
               widths = c(13, 24, 60, 18, 60, 30, 22, 60))
  invisible(NULL)
}

#' Write the source data workbook
#'
#' @param run_dates A data frame, or the path of a CSV file, with the columns
#'   `batch` and `run_date`; by default the shipped results/run_dates.csv. A
#'   re-simulation passes the dates of its own run.
#' @param path Character scalar, the file to write.
#' @param pd Named list of the 17 panel data frames, passed to
#'   `source_data_tab()`.
#' @param sweep_seeds Integer scalar, seeds per parameter set of the sweep,
#'   shown on the Index; `NULL` leaves it out.
#' @return `path`, invisibly. The workbook holds the Index and the 24 data
#'   tabs of `SOURCE_DATA_TABS`.
write_source_data <- function(run_dates = file.path(PROJ, "results",
                                                    "run_dates.csv"),
                              path = SOURCE_DATA_XLSX,
                              pd = panel_data,
                              sweep_seeds = if (exists("SWEEP_SEEDS"))
                                SWEEP_SEEDS else NULL) {
  if (!identical(sort(names(TAB_QUESTIONS)), sort(SOURCE_DATA_TABS)))
    stop("TAB_QUESTIONS does not list the tabs of SOURCE_DATA_TABS.",
         call. = FALSE)
  dates <- read_run_dates(run_dates)

  # The results folder is named on the Index relative to the repository root,
  # so that no path of the machine that wrote the workbook appears in it.
  results <- HR_RESULTS
  if (startsWith(results, paste0(PROJ, "/")))
    results <- substring(results, nchar(PROJ) + 2)
  if (grepl("^(/|~|[A-Za-z]:)", results))
    stop("HR_RESULTS must lie inside the repository root.", call. = FALSE)

  message("source_data.R: reading the tables of ", length(SOURCE_DATA_TABS),
          " tabs ...")
  tables <- setNames(lapply(SOURCE_DATA_TABS, source_data_tab, pd = pd),
                     SOURCE_DATA_TABS)

  wb   <- createWorkbook(creator = "")
  bold <- createStyle(textDecoration = "bold")
  addWorksheet(wb, "Index")
  for (tab in SOURCE_DATA_TABS) {
    addWorksheet(wb, tab)
    left_out <- write_data_tab(wb, tab, tables[[tab]],
                               tab_run_date(tab, dates), bold)
    if (left_out)
      message("  ", tab, ": per-replicate values left out (",
              nrow(tables[[tab]]$per_replicate), " rows)")
  }
  write_index_tab(wb, index_rows(tables, dates, results, sweep_seeds), bold)

  # openxlsx declares a drawing part for every worksheet but writes it only
  # for a sheet that holds an image, so a workbook without images refers to
  # parts that are not in the file. Readers that follow every reference, among
  # them Python's openpyxl, stop on them. openxlsx also gives every worksheet
  # a printer-settings part. That part is a template shipped with the library
  # and holds printer settings that are not this workbook's, so the workbook
  # is written without it. The references to both kinds of part are removed
  # here: the relationship entries of each sheet, the content-type entries
  # and, for the printer settings, the r:id attribute of each sheet's page
  # setup, whose other attributes stay. A sheet with no relationship left gets
  # no relationships file. No cell, style or other sheet setting is affected.
  stopifnot(all(vapply(wb$drawings, function(d) all(d == ""), logical(1))))
  for (i in seq_along(wb$worksheets_rels))
    wb$worksheets_rels[[i]] <- grep(
      "relationships/(drawing|vmlDrawing|printerSettings)\"",
      wb$worksheets_rels[[i]], invert = TRUE, value = TRUE)
  wb$Content_Types <- grep("/xl/drawings/|printerSettings", wb$Content_Types,
                           invert = TRUE, value = TRUE)
  for (i in seq_along(wb$worksheets))
    wb$worksheets[[i]]$pageSetup <- sub(" r:id=\"[^\"]*\"", "",
                                        wb$worksheets[[i]]$pageSetup)
  # openxlsx also records the used range of every sheet as A1. Readers that
  # trust that record, such as openpyxl in read-only mode, would see a single
  # cell, so the range actually written is recorded instead.
  for (i in seq_along(wb$worksheets)) {
    sd <- wb$worksheets[[i]]$sheet_data
    wb$worksheets[[i]]$dimension <- sprintf(
      "<dimension ref=\"A1:%s%d\"/>", int2col(max(sd$cols)), max(sd$rows))
  }
  saveWorkbook(wb, path, overwrite = TRUE)
  # saveWorkbook() writes the printer-settings parts themselves whether or not
  # anything refers to them, and has no setting that leaves them out. The
  # saved file is therefore unpacked, those parts are deleted and the other
  # parts are packed again as they are, as saveWorkbook() packs them: files
  # only, at compression level 6.
  unpacked <- tempfile("source_data_")
  repacked <- tempfile(fileext = ".xlsx")
  zip::unzip(path, exdir = unpacked)
  unlink(file.path(unpacked, "xl", "printerSettings"), recursive = TRUE)
  zip::zipr(repacked, list.files(unpacked, full.names = TRUE),
            include_directories = FALSE, compression_level = 6)
  stopifnot(!any(grepl("printerSettings", zip::zip_list(repacked)$filename)))
  stopifnot(file.copy(repacked, path, overwrite = TRUE))
  unlink(c(unpacked, repacked), recursive = TRUE)
  message("source_data.R: wrote ", file.path(basename(dirname(path)),
                                             basename(path)))
  invisible(path)
}

# Write ----

write_source_data(
  run_dates = if (exists("SOURCE_DATA_RUN_DATES")) SOURCE_DATA_RUN_DATES else
    file.path(PROJ, "results", "run_dates.csv")
)
