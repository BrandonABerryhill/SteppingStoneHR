# 08_notebook.R
#
# The steps of the notebook SteppingStone.Rmd that are more than one call to
# another module: the table of model parameters, the comparison of the two
# engines, the helper that loads a batch or runs it, the frequency summary, the
# sensitivity sweep with its PAWN indices, and the table of run dates. The
# notebook holds the settings (horizons, seed counts, sweep ranges, seeds) and
# passes them to these functions, so that a script can call the same functions
# with fewer seeds.
#
# Two modes. Functions with a `resimulate` argument read stored results when it
# is FALSE and run simulations when it is TRUE. In both cases the folder is the
# global RESULTS_DIR of 00_config.R, which run_batch() also uses: the folder of
# shipped results when reading, and a separate folder when simulating, because
# run_batch() resumes from any batch it finds. Nothing is written when
# `resimulate` is FALSE.
#
# This file is sourced after 00_config.R to 07_engines.R. Times are in hours.
#
# Part of the code for "Unstable Heteroresistance as a Stepping Stone to Stable
# Resistance via Compensatory Mutations" (Berryhill and Andersson).
#
# Defines: condition_tables(), engine_gate(), load_or_run_batch(),
# frequency_summary(), interval_overlap(), same_frequency_tables(),
# sweep_design(), pawn_tables(), sensitivity_sweep() and resim_run_dates().

# Model parameters ----

#' Parameter values of several conditions, side by side
#'
#' @param conds Named list of conditions from `make_condition()`, named by
#'   scenario.
#' @return A list of two data frames. `scalars` has one row per condition, with
#'   the column `condition` and one column per scalar parameter of `p` (`DT`,
#'   `transfer`, `C`, `e`, `kmon`, `AMAX`, `dil`, `d`, `KB`, `mu`, `gamma`,
#'   `omega`). `by_state` has one row per state, in the order of `states`,
#'   with the columns `state`, `copy_number` and `y0` (initial density, cells
#'   per mL) and the per-state parameters `v`, `vm`, `kH`, `M` and `P`: a
#'   parameter or initial density that is the same in every condition has one
#'   column, and one that differs has a column per condition, named with the
#'   scenario as suffix. Units are those given in 01_parameters.R.
condition_tables <- function(conds) {
  scalar_names <- c("DT", "transfer", "C", "e", "kmon", "AMAX", "dil", "d",
                    "KB", "mu", "gamma", "omega")
  scalars <- data.frame(condition = names(conds), stringsAsFactors = FALSE)
  for (name in scalar_names)
    scalars[[name]] <- unname(vapply(conds, function(cond) cond$p[[name]],
                                     numeric(1)))

  per_state <- list(y0 = lapply(conds, function(cond) unname(cond$y0$pop)))
  for (name in c("v", "vm", "kH", "M", "P"))
    per_state[[name]] <- lapply(conds, function(cond) cond$p[[name]])
  columns <- list(state = states, copy_number = copy_number)
  for (name in names(per_state)) {
    values <- per_state[[name]]
    same   <- all(vapply(values, function(x) identical(x, values[[1]]),
                         logical(1)))
    if (same) {
      columns[[name]] <- values[[1]]
    } else {
      for (sc in names(values))
        columns[[paste(name, sc, sep = "_")]] <- values[[sc]]
    }
  }
  list(scalars = scalars,
       by_state = as.data.frame(columns, stringsAsFactors = FALSE))
}

# Engine verification ----

#' Compare the R reference and the compiled engine without the random layer
#'
#' Runs Scenario 1 once on each engine with `stochastic = FALSE`, on the same
#' output grid, and takes the largest relative difference over every numeric
#' column the two results share, which is all 25, and every row. The relative
#' difference of two values is their absolute difference divided by the
#' absolute value from the R reference, or by 1e-12 where that value is below
#' 1e-12. The seed passed to the engines is 1; with the random layer off it
#' has no effect on either result. The R reference engine still calls
#' `set.seed(1)`, which resets R's random number generator in the calling
#' session.
#'
#' @param t_end Numeric scalar, simulated time in hours.
#' @param n_out Positive integer, the approximate number of regularly spaced
#'   rows both engines record.
#' @param tol Numeric scalar, the largest relative difference accepted.
#' @param halt Logical scalar; `TRUE` stops with an error if the difference
#'   exceeds `tol`.
#' @param record Path of an `.rds` file in which `list(relerr, tol)` is saved,
#'   or `NULL` to save nothing.
#' @return A list: `relerr` (the largest relative difference), `tol`, `pass`
#'   (`relerr <= tol`), `t_end`, `n_out`, `n_columns` and `n_rows` (the number
#'   of columns and of rows compared). With `halt = TRUE` the function returns
#'   only if `pass` is `TRUE`.
engine_gate <- function(t_end, n_out, tol, halt = TRUE, record = NULL) {
  if (!exists("run_model_cpp"))
    stop("The compiled engine is not loaded: the comparison needs the hrsim ",
         "package.", call. = FALSE)
  cond <- make_condition("S1"); tr <- build_transitions()
  dR <- run_model_R  (cond$p, tr, cond$y0, t_end = t_end, seed = 1,
                      n_out = n_out, stochastic = FALSE)
  dC <- run_model_cpp(cond$p, tr, cond$y0, t_end = t_end, seed = 1,
                      n_out = n_out, stochastic = FALSE)
  cols <- intersect(names(dR), names(dC))
  cols <- cols[vapply(cols, function(c) is.numeric(dR[[c]]), logical(1))]
  nrow_min <- min(nrow(dR), nrow(dC))
  relerr <- 0
  for (cc in cols) {
    a <- dR[[cc]][seq_len(nrow_min)]; b <- dC[[cc]][seq_len(nrow_min)]
    relerr <- max(relerr, max(abs(a - b) / pmax(abs(a), 1e-12)))
  }
  pass <- relerr <= tol
  if (!is.null(record)) saveRDS(list(relerr = relerr, tol = tol), record)
  if (halt && !isTRUE(pass))
    stop(sprintf(paste("The compiled engine differs from the R reference by",
                       "more than the tolerance (maximum relative error %.3e,",
                       "tolerance %.1e). Results of the compiled engine",
                       "should not be used."), relerr, tol), call. = FALSE)
  list(relerr = relerr, tol = tol, pass = pass, t_end = t_end, n_out = n_out,
       n_columns = length(cols), n_rows = nrow_min)
}

# Batches ----

#' Load a batch from the results folder, or run it
#'
#' With `resimulate = FALSE` the batch stored under RESULTS_DIR/<scenario>/<tag>
#' is read and nothing is written. With `resimulate = TRUE` the batch is run
#' there with `run_batch()`, which resumes from the replicates already on disk,
#' and its classification is saved beside it as `classification.rds`. In both
#' cases every replicate is classified on its transfer rows with
#' `classify_batch()`.
#'
#' @param scenario Character scalar passed to `make_condition()`: "S1", "S2"
#'   or "omega0". The condition is built with the default parameter values.
#' @param tag Character scalar, the batch's folder name within the scenario
#'   folder.
#' @param t_end Numeric scalar, simulated time in hours.
#' @param seeds Integer vector of replicate seeds. Each seed is one replicate,
#'   and a replicate's random stream depends only on its own seed.
#' @param keep "cycle_end" keeps the rows recorded at transfers only; "all"
#'   keeps every recorded row.
#' @param n_out `NULL` for the engine's default output grid, or a positive
#'   integer (see `run_single()`).
#' @param engine "cpp" for the compiled engine or "R" for the R reference.
#'   For "R" the global `run_model` is set to `run_model_R` while the batch
#'   runs and put back afterwards. Used only when `resimulate` is `TRUE`.
#' @param resimulate Logical scalar: `FALSE` loads the batch, `TRUE` runs it.
#' @param n_workers Positive integer, parallel workers of a run; does not
#'   affect results.
#' @return A list: `tag`, `batch_dir` (the batch's folder), `sims` (list of
#'   data frames, one per replicate among `seeds`, named by seed and ordered by
#'   increasing seed, as returned by `run_batch()`) and `cls` (data frame, one
#'   row per replicate, the columns of `classify_outcome()`).
load_or_run_batch <- function(scenario, tag, t_end, seeds,
                              keep = "cycle_end", n_out = NULL,
                              engine = c("cpp", "R"),
                              resimulate = FALSE, n_workers = N_WORKERS) {
  engine <- match.arg(engine)
  if (resimulate) {
    if (engine == "cpp" && !exists("run_model_cpp"))
      stop("The compiled engine is not loaded: batch '", tag,
           "' needs the hrsim package.", call. = FALSE)
    if (engine == "R") {
      rm_save <- run_model; run_model <<- run_model_R
      on.exit(run_model <<- rm_save, add = TRUE)
    }
    b <- run_batch(make_condition(scenario), t_end = t_end,
                   seeds = seeds, n_workers = n_workers,
                   tag = tag, keep = keep, n_out = n_out)
    cls <- classify_batch(b$sims)
    saveRDS(cls, file.path(b$batch_dir, "classification.rds"))
    return(list(tag = tag, batch_dir = b$batch_dir, sims = b$sims, cls = cls))
  }

  batch_dir <- file.path(RESULTS_DIR, scenario, tag)
  stored    <- file.path(batch_dir, "sims.rds")
  if (!file.exists(stored))
    stop("Batch '", tag, "' was not found: ", stored, " does not exist.",
         call. = FALSE)
  done <- readRDS(stored)
  # The stored file is run_batch()'s checkpoint. The replicates are selected
  # from it as run_batch() selects them: failed replicates are left out, and
  # the requested seeds are returned in increasing order.
  good <- done[!vapply(done, function(x) inherits(x, "failed_seed"),
                       logical(1))]
  want <- as.character(sort(as.integer(intersect(names(good),
                                                 as.character(seeds)))))
  if (length(want) < length(seeds))
    message(sprintf("  [%s] %d of the %d requested seeds are stored", tag,
                    length(want), length(seeds)))
  sims <- good[want]
  list(tag = tag, batch_dir = batch_dir, sims = sims,
       cls = classify_batch(sims))
}

# Outcome frequencies ----

#' Outcome frequencies of the frequency batches, in the form the run stored
#'
#' @param batches Named list with one element per condition, named by
#'   scenario. Each element is a list with `cpp`, the compiled frequency batch
#'   as returned by `load_or_run_batch()`, and optionally `R`, the R-reference
#'   batch.
#' @param save_path Path of the `.rds` file to write, or `NULL` to write
#'   nothing.
#' @return A named list with one element per condition, each a list of `cpp`
#'   (the table of `freq_table()` for the compiled batch, with percentages and
#'   bounds rounded to three decimals), `dir_cpp` (the batch's folder) and,
#'   where an R-reference batch was given, `R` and `dir_R`. This is the content
#'   of `frequency_summary.rds`, from which Fig 1C is drawn.
frequency_summary <- function(batches, save_path = NULL) {
  FREQ <- list()
  for (sc in names(batches)) {
    FREQ[[sc]]$cpp     <- freq_table(batches[[sc]]$cpp$cls)
    FREQ[[sc]]$dir_cpp <- batches[[sc]]$cpp$batch_dir
    if (!is.null(batches[[sc]]$R)) {
      FREQ[[sc]]$R     <- freq_table(batches[[sc]]$R$cls)
      FREQ[[sc]]$dir_R <- batches[[sc]]$R$batch_dir
    }
  }
  if (!is.null(save_path)) saveRDS(FREQ, save_path)
  FREQ
}

#' Do the confidence intervals of the two engines overlap?
#'
#' For each condition with an R-reference batch and each outcome, compares the
#' 95% Clopper-Pearson interval of the compiled batch with that of the
#' R-reference batch, using the rounded bounds of `frequency_summary()`.
#'
#' @param FREQ List returned by `frequency_summary()`.
#' @return A data frame with one row per condition and outcome: `condition`,
#'   `outcome`, `cpp_low`, `cpp_high`, `R_low`, `R_high` (bounds in percent)
#'   and `overlap` (`TRUE` if the two intervals share at least one value);
#'   `NULL` if no condition has an R-reference batch.
interval_overlap <- function(FREQ) {
  rows <- list()
  for (sc in names(FREQ)) {
    if (is.null(FREQ[[sc]]$R)) next
    for (o in c("climb","low-copy","extinction")) {
      a <- FREQ[[sc]]$cpp[FREQ[[sc]]$cpp$outcome==o,]
      r <- FREQ[[sc]]$R  [FREQ[[sc]]$R$outcome==o,]
      overlap <- !(a$ci_high < r$ci_low || r$ci_high < a$ci_low)
      rows[[length(rows) + 1]] <- data.frame(
        condition = sc, outcome = o,
        cpp_low = a$ci_low, cpp_high = a$ci_high,
        R_low = r$ci_low, R_high = r$ci_high,
        overlap = overlap, stringsAsFactors = FALSE)
    }
  }
  if (length(rows) == 0) return(NULL)
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Are two frequency summaries the same tables?
#'
#' @param FREQ,stored Lists in the form of `frequency_summary()`, for example
#'   the one computed in the session and the one read from
#'   `frequency_summary.rds`.
#' @return `TRUE` if both hold the same conditions and, for every condition,
#'   identical `cpp` and `R` tables; the folder names `dir_cpp` and `dir_R`
#'   are not compared.
same_frequency_tables <- function(FREQ, stored) {
  tables <- function(x) lapply(x, function(s) s[c("cpp", "R")])
  identical(tables(FREQ), tables(stored))
}

# Sensitivity sweep ----

#' Draw the Latin hypercube design of the sweep
#'
#' The seed is set immediately before the sample is drawn and governs the whole
#' design; no replicate depends on it. The design also depends on the kind of
#' R's random number generator, so the function stops unless the session is on
#' Mersenne-Twister, R's default and the kind with which the design of the
#' paper was drawn. A session leaves the default kind when
#' `furrr::future_map()` is called with `seed = TRUE`, as `run_batch()` does,
#' before any call to `set.seed()`; a script that sets a seed at its start
#' stays on the default.
#'
#' @param n_sets Positive integer, number of parameter sets.
#' @param ranges Named list of `c(min, max)` per parameter, as for
#'   `lhs_sample()`.
#' @param log Character vector, the parameters sampled on a log10 scale.
#' @param seed Integer scalar passed to `set.seed()`.
#' @return The data frame of `lhs_sample()`: `n_sets` rows, one column per
#'   parameter.
sweep_design <- function(n_sets, ranges, log, seed) {
  if (RNGkind()[1] != "Mersenne-Twister")
    stop("The sweep design is drawn with R's default generator, ",
         "Mersenne-Twister, but the session is on ", RNGkind()[1],
         ". Call set.seed() before the first batch of the session.",
         call. = FALSE)
  set.seed(seed)
  lhs_sample(n_sets, ranges, log = log)
}

#' PAWN indices of the three outcome frequencies
#'
#' Parameter sets without frequencies, which are sets in which no replicate
#' succeeded, are left out. `pawn_given_data()` is called with its default
#' number of quantile bins.
#'
#' @param sweep Data frame of the sweep: one row per parameter set with the
#'   parameter columns and `f_climb`, `f_lowcopy` and `f_extinct` (0 to 1).
#' @param parameters Character vector, the parameter columns of `sweep`.
#' @return A list of three data frames, `climb`, `lowcopy` and `extinct`, each
#'   as returned by `pawn_given_data()`: one row per parameter, columns
#'   `median` and `max`.
pawn_tables <- function(sweep, parameters) {
  ok  <- stats::complete.cases(sweep[, c("f_climb","f_lowcopy","f_extinct")])
  inp <- sweep[ok, parameters, drop = FALSE]
  list(
    climb   = pawn_given_data(inp, sweep$f_climb[ok]),
    lowcopy = pawn_given_data(inp, sweep$f_lowcopy[ok]),
    extinct = pawn_given_data(inp, sweep$f_extinct[ok]))
}

#' Load the sensitivity sweep of Scenario 1 and its PAWN indices, or run them
#'
#' The tables are kept in RESULTS_DIR/sensitivity/sweep_S1/ as `sweep.rds`,
#' `sweep.csv` and `pawn.rds`. With `resimulate = FALSE` the two `.rds` files
#' are read and nothing is written. With `resimulate = TRUE` the design is
#' drawn with `sweep_design()`, every set is run with `run_sweep()`, which
#' stores set i as the batch "sweep_set" followed by i in four digits and
#' keeps the frequencies of finished sets in
#' RESULTS_DIR/sensitivity/sweep_S1_progress.rds, and the tables are written.
#' The replicates of every set have the seeds 1 to `seeds_per_set`, and each
#' replicate's random stream depends only on its own seed.
#'
#' @param n_sets Positive integer, number of parameter sets of the design.
#' @param seeds_per_set Positive integer, replicates per set.
#' @param ranges,log Passed to `sweep_design()`; `names(ranges)` are also the
#'   inputs of the PAWN analysis.
#' @param design_seed Integer scalar, the seed of the design.
#' @param resimulate Logical scalar: `FALSE` loads the tables, `TRUE` runs the
#'   sweep.
#' @param n_workers Positive integer, parallel workers of a run; does not
#'   affect results.
#' @param first_sets Positive integer, at most `n_sets`. A value below `n_sets`
#'   runs only that many sets from the start of the full design, as a test
#'   does. The PAWN indices are then not computed, because they need the whole
#'   design, and `pawn.rds` is not written.
#' @return A list: `dir` (the folder of the tables), `sweep` (the design with
#'   the columns `f_climb`, `f_lowcopy` and `f_extinct` added, one row per set
#'   run), `pawn` (list from `pawn_tables()`, or `NULL` for a partial run) and
#'   `checkpoint` (path of the progress file of `run_sweep()`).
sensitivity_sweep <- function(n_sets, seeds_per_set, ranges, log, design_seed,
                              resimulate = FALSE, n_workers = N_WORKERS,
                              first_sets = n_sets) {
  sdir       <- file.path(RESULTS_DIR, "sensitivity", "sweep_S1")
  sweep_ckpt <- file.path(RESULTS_DIR, "sensitivity", "sweep_S1_progress.rds")
  if (!resimulate) {
    stored <- file.path(sdir, c("sweep.rds", "pawn.rds"))
    if (!all(file.exists(stored)))
      stop("The sweep tables were not found: ",
           paste(stored[!file.exists(stored)], collapse = ", "),
           " missing.", call. = FALSE)
    return(list(dir = sdir, sweep = readRDS(stored[1]),
                pawn = readRDS(stored[2]), checkpoint = sweep_ckpt))
  }

  # The design is drawn from a fixed seed, so a resumed sweep has the same
  # parameter sets, row for row, as the one that wrote the checkpoint.
  X <- sweep_design(n_sets, ranges, log, design_seed)
  X <- X[seq_len(first_sets), , drop = FALSE]
  dir.create(dirname(sweep_ckpt), recursive = TRUE, showWarnings = FALSE)
  sweep <- run_sweep("S1", X, seeds_per_set = seeds_per_set,
                     n_workers = n_workers, ckpt_path = sweep_ckpt)
  dir.create(sdir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(sweep, file.path(sdir, "sweep.rds"))
  write.csv(sweep, file.path(sdir, "sweep.csv"), row.names = FALSE)
  pawn <- NULL
  if (first_sets == n_sets) {
    pawn <- pawn_tables(sweep, names(ranges))
    saveRDS(pawn, file.path(sdir, "pawn.rds"))
  } else {
    message(sprintf("  sweep: %d of %d sets run; PAWN indices not computed",
                    first_sets, n_sets))
  }
  list(dir = sdir, sweep = sweep, pawn = pawn, checkpoint = sweep_ckpt)
}

# Run dates ----

#' Run dates of a re-simulation, in the form of results/run_dates.csv
#'
#' The source data workbook shows the date on which the code behind each tab
#' was run. For a re-simulation these dates are taken from its own records.
#' A batch is dated by the earliest row the manifest holds for it, because
#' `run_batch()` appends a row at the end of every call and a batch that was
#' resumed therefore has several. The sweep has no manifest row and is dated
#' by the time its progress file was last written, which is when its last set
#' finished. The two supplementary tables report the parameters of the
#' Scenario 1 frequency batch and take its date. Entries that no simulation
#' produces are copied from the table of the shipped results.
#'
#' @param manifest Path of the manifest of the re-simulation.
#' @param sweep_checkpoint Path of the progress file of the sweep, or `NULL`
#'   if no sweep was run.
#' @param shipped Path of the run-dates table of the shipped results, with the
#'   columns `batch`, `run_date` and `source`.
#' @param copied Character vector, the entries copied from `shipped`: by
#'   default the long run behind Supplementary Fig 4 and Fig 3C, which is
#'   computed from parameters.
#' @param path Path of the CSV file to write, or `NULL` to write nothing.
#' @return A data frame with the columns `batch`, `run_date` (YYYY-MM-DD) and
#'   `source`: a row per manifest batch, in order of first appearance, dated by
#'   its earliest manifest row; the sweep, dated by its progress file, if that
#'   exists; the `shipped` rows named in `copied`, dated and ordered as there;
#'   and, given a row "S1_freq_cpp", the two supplementary tables with its date.
resim_run_dates <- function(manifest = MANIFEST, sweep_checkpoint = NULL,
                            shipped, copied = c("long_run", "Fig 3C"),
                            path = NULL) {
  m <- read.csv(manifest, stringsAsFactors = FALSE)
  batches <- unique(m$batch_id)
  first   <- vapply(batches, function(b) min(m$datetime[m$batch_id == b]),
                    character(1))
  out <- data.frame(batch = batches, run_date = substr(first, 1, 10),
                    source = "manifest", stringsAsFactors = FALSE)
  if (!is.null(sweep_checkpoint) && file.exists(sweep_checkpoint))
    out <- rbind(out, data.frame(
      batch = "sweep_S1",
      run_date = format(file.mtime(sweep_checkpoint), "%Y-%m-%d"),
      source = "file timestamp of the sweep progress file",
      stringsAsFactors = FALSE))
  ship <- read.csv(shipped, colClasses = "character", stringsAsFactors = FALSE)
  out  <- rbind(out, ship[ship$batch %in% copied,
                          c("batch", "run_date", "source")])
  if ("S1_freq_cpp" %in% out$batch)
    out <- rbind(out, data.frame(
      batch = c("Supp Table 1", "Supp Table 2"),
      run_date = out$run_date[out$batch == "S1_freq_cpp"],
      source = "date of the Scenario 1 frequency batch",
      stringsAsFactors = FALSE))
  rownames(out) <- NULL
  if (!is.null(path)) write.csv(out, path, row.names = FALSE, quote = FALSE)
  out
}
