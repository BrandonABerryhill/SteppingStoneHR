# 05_run.R
#
# The batch runner. run_batch() runs many replicates of one condition in
# parallel and saves a checkpoint after every wave of replicates. When its
# manifest argument is TRUE, the default, it appends a row to the manifest at
# each call; freq_for_set() in 06_sensitivity.R passes FALSE, so the batches
# of a sweep add no row. A batch is written to RESULTS_DIR/<scenario>/<tag>/,
# so the tag identifies it: calling run_batch() again with the same scenario
# and tag loads the replicates already on disk and runs only the missing
# seeds. Each batch therefore needs a tag of its own. A replicate that stops
# with an error is recorded and left out of the result instead of ending the
# batch.
#
# Seeds and reproducibility. A replicate is identified by its integer seed,
# which is passed to the engine. The compiled engine seeds its own generator
# from it, and the R reference engine calls set.seed() with it, so a
# replicate's random stream depends only on its own seed and not on the number
# of workers, the size of a wave or the order in which replicates finish.
#
# Part of the code for "Unstable Heteroresistance as a Stepping Stone to Stable
# Resistance via Compensatory Mutations" (Berryhill and Andersson).
#
# Defines: run_single(), .transitions_cache(), run_batch() and
# write_manifest_row().

# One replicate ----

#' Run one replicate with the active engine
#'
#' @param cond Condition from `make_condition()`.
#' @param t_end Numeric scalar, simulated time in hours.
#' @param seed Integer scalar, the replicate's seed. It governs all of the
#'   replicate's random draws. With the compiled engine the replicate's
#'   random stream depends on nothing else. With the R reference engine it
#'   also depends on the kind of R's random number generator in force, which
#'   inside `run_batch()` is L'Ecuyer-CMRG (see the comment on furrr's seed
#'   option there).
#' @param keep "cycle_end" keeps only the rows recorded at transfers, which
#'   are the rows the classifier needs; "all" keeps every recorded row.
#' @param n_out `NULL`, or a positive integer giving the approximate number of
#'   regularly spaced rows the engine records. With `NULL` the engine is
#'   called without `n_out` and uses its own default, which differs between
#'   the engines: 400 regular rows for the R reference, the initial row and
#'   the transfer rows only for the compiled engine.
#' @return The engine's data frame (25 columns, see 03_engine_R.R) with a
#'   `seed` column added. With `keep = "cycle_end"` it is an error if the run
#'   has no transfer row.
run_single <- function(cond, t_end, seed, keep = c("cycle_end", "all"),
                       n_out = NULL) {
  keep <- match.arg(keep)
  d <- if (is.null(n_out)) {
    run_model(cond$p, .transitions_cache(), cond$y0, t_end = t_end, seed = seed)
  } else {
    run_model(cond$p, .transitions_cache(), cond$y0, t_end = t_end, seed = seed,
              n_out = n_out)
  }
  if (keep == "cycle_end") {
    d <- d[d$cycle_end == 1, ]
    if (nrow(d) == 0)
      stop(sprintf(paste("no end-of-cycle rows: t_end (%g h) must span at least one full",
                         "transfer interval (%g h) and should be the classification horizon."),
                   t_end, cond$p$transfer), call. = FALSE)
  }
  d$seed <- seed
  d
}

# The transition table does not depend on the condition, so it is built once
# and reused.
.transitions_cache <- local({
  tr <- NULL
  function() { if (is.null(tr)) tr <<- build_transitions(); tr }
})

# One batch ----

#' Run a batch of replicates with checkpoint and resume
#'
#' @param cond Condition from `make_condition()`. `cond$p$scenario` names the
#'   scenario folder.
#' @param t_end Numeric scalar, simulated time in hours.
#' @param seeds Integer vector of replicate seeds. Each seed is one replicate,
#'   and a replicate's random stream depends only on its own seed.
#' @param n_workers Positive integer, number of parallel R sessions; 1 runs
#'   the replicates in the calling session. It does not affect results.
#' @param batch_sz Positive integer, replicates per wave; a checkpoint is saved
#'   after each wave. `NULL` uses `n_workers`.
#' @param tag Character scalar, the batch's folder name within its scenario
#'   folder. It is the key for resuming.
#' @param keep Passed to `run_single()`: "cycle_end" or "all".
#' @param manifest Logical scalar; `TRUE` appends a row to the manifest.
#' @param n_out Passed to `run_single()`; `NULL` leaves the engine's default.
#' @return A list: `batch_id` (the tag), `batch_dir` (the folder), `sims` (a
#'   list of data frames, one per successful replicate among `seeds`, named by
#'   seed and ordered by increasing seed), `cond`, and `n_failed` (number of
#'   replicates on disk that stopped with an error). In `batch_dir` it writes
#'   `sims.rds`, the checkpoint holding every replicate run so far (a failed
#'   one as a marker), and `condition.rds`, which is rewritten at every call.
run_batch <- function(cond, t_end, seeds,
                      n_workers = 4, batch_sz = NULL,
                      tag = cond$p$scenario, keep = "cycle_end",
                      manifest = TRUE, n_out = NULL) {
  if (is.null(batch_sz)) batch_sz <- max(1L, n_workers)
  batch_id  <- tag
  batch_dir <- file.path(RESULTS_DIR, cond$p$scenario, tag)
  dir.create(batch_dir, showWarnings = FALSE, recursive = TRUE)
  ckpt <- file.path(batch_dir, "sims.rds")

  done <- list()
  if (file.exists(ckpt)) {
    done <- readRDS(ckpt)
    message(sprintf("  resuming '%s': %d seed(s) already on disk", tag, length(done)))
  }
  remaining <- setdiff(as.character(seeds), names(done))

  # A replicate that fails is stored as a marker in the checkpoint, so that a
  # resumed batch does not run it again, and is left out of the returned sims.
  run_one <- function(s) {
    tryCatch(run_single(cond, t_end, as.integer(s), keep = keep, n_out = n_out),
             error = function(e)
               structure(list(failed = TRUE, seed = as.integer(s),
                              msg = conditionMessage(e)), class = "failed_seed"))
  }

  if (length(remaining) > 0) {
    # Putting the plan back to sequential on exit, also after an error, closes
    # the worker sessions. A plan the caller had set is not restored.
    if (n_workers <= 1) future::plan(future::sequential)
    else                future::plan(future::multisession, workers = n_workers)
    on.exit(future::plan(future::sequential), add = TRUE)
    message(sprintf("  [%s] running on %d worker(s)", tag, future::nbrOfWorkers()))

    n_start <- length(done)
    n_tot   <- length(seeds)
    batches <- split(remaining, ceiling(seq_along(remaining) / batch_sz))
    t_start <- Sys.time()
    for (b in seq_along(batches)) {
      sb  <- batches[[b]]
      # seed = TRUE makes furrr hand each replicate a stream of R's
      # L'Ecuyer-CMRG generator. Neither engine keeps that stream: the compiled
      # engine has a generator of its own, and the R reference engine re-seeds
      # R's generator with set.seed(). Both use the replicate's own seed. The
      # option does set the kind of R's generator, so inside a batch the R
      # reference engine seeds L'Ecuyer-CMRG, and its draws differ from those
      # of a direct call in a session that has R's default kind of generator.
      res <- furrr::future_map(sb, run_one,
                               .options = furrr::furrr_options(seed = TRUE))
      names(res) <- sb
      done <- c(done, res)
      saveRDS(done, ckpt)
      n_done <- length(done)
      # The remaining time is estimated from the replicates run in this call.
      # Times are in minutes.
      el     <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
      per    <- el / max(1L, n_done - n_start)
      eta    <- per * (n_tot - n_done)
      message(sprintf("  [%s] wave %d/%d | %d/%d seeds (%.0f%%) | %.1f min elapsed | ETA %.1f min",
                      tag, b, length(batches), n_done, n_tot,
                      100 * n_done / n_tot, el, eta))
    }
  }

  failed <- done[vapply(done, function(x) inherits(x, "failed_seed"), logical(1))]
  if (length(failed))
    message(sprintf("  [%s] WARNING: %d seed(s) failed and were excluded (e.g. seed %s: %s)",
                    tag, length(failed), failed[[1]]$seed, failed[[1]]$msg))
  good <- done[!vapply(done, function(x) inherits(x, "failed_seed"), logical(1))]
  want <- as.character(sort(as.integer(intersect(names(good), as.character(seeds)))))
  sims <- good[want]

  saveRDS(cond, file.path(batch_dir, "condition.rds"))
  if (manifest) write_manifest_row(cond, batch_id, batch_dir, t_end, length(sims))
  list(batch_id = batch_id, batch_dir = batch_dir, sims = sims, cond = cond,
       n_failed = length(failed))
}

# Manifest ----

#' Append one row describing a batch to the manifest
#'
#' @param cond Condition from `make_condition()`.
#' @param batch_id Character scalar, the batch's tag.
#' @param batch_dir Character scalar, the batch's folder.
#' @param t_end Numeric scalar, simulated time in hours.
#' @param n_seeds Integer scalar, number of replicates returned.
#' @return Invisibly, the one-row data frame that was appended to `MANIFEST`,
#'   with columns `batch_id`, `datetime` (system clock at the time of the
#'   call), `scenario`, `horizon_h`, `n_seeds`, `omega`, `amax`,
#'   `cost_per_copy`, `classify_horizon`, `climb_threshold`, `code_version`
#'   and `results_path`. The header line is written when the file is created.
write_manifest_row <- function(cond, batch_id, batch_dir, t_end, n_seeds) {
  p <- cond$p
  row <- data.frame(
    batch_id      = batch_id,
    datetime      = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    scenario      = p$scenario,
    horizon_h     = t_end,
    n_seeds       = n_seeds,
    omega         = p$omega,
    amax          = p$AMAX,
    cost_per_copy = p$cost_per_copy,
    classify_horizon = CLASSIFY_HORIZON,
    climb_threshold = CLIMB_THRESHOLD,
    code_version  = CODE_VERSION,
    results_path  = batch_dir,
    stringsAsFactors = FALSE)
  hdr <- !file.exists(MANIFEST)
  write.table(row, MANIFEST, sep = ",", row.names = FALSE,
              col.names = hdr, append = !hdr)
  invisible(row)
}
