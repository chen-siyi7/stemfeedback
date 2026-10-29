simulate_competition <- function(model, initial, n = 100L, seed = 1L,
                                 workers = 1L, reference = NULL,
                                 enlarged_reference = NULL, max_events = 2e7,
                                 time_limit = Inf, observation_times = 0,
                                 keep_paths = FALSE) {
  stopifnot(length(n) == 1L, is.finite(n), n >= 1, n == round(n),
            length(seed) == 1L, is.finite(seed), seed >= 1, seed == round(seed),
            seed + n - 1 < .Machine$integer.max,
            length(workers) == 1L, is.finite(workers), workers >= 1,
            workers == round(workers))
  validate_state(initial, model)
  if (!is.null(reference)) validate_reference(reference, model)
  if (!is.null(enlarged_reference)) validate_reference(enlarged_reference, model)
  paths <- parallel_paths(seq_len(n), function(i)
    simulate_path(initial, model, as.integer(seed + i - 1L),
      max_events = max_events, time_limit = time_limit,
      observation_times = observation_times, reference = reference,
      enlarged_reference = enlarged_reference), workers = workers)
  if (any(vapply(paths, inherits, FALSE, what = "try-error")))
    stop("A simulation worker failed.")
  records <- bind_rows(lapply(seq_along(paths), function(i)
    path_record(paths[[i]], i, "custom")))
  case <- list(id = "custom", suite = "custom", model = model, initial = initial,
               stop_rule = "first_stem_loss", note = "User-specified competition")
  summary <- summarize_case(records, case, "custom")
  feedback_summary <- if (is.null(reference)) NULL else list(
    reference_initial = records$reference_initial[1],
    mean_integral = mean_interval(records$integral),
    mean_residual = mean_interval(records$residual),
    outer_exits = sum(records$status == "reference_exit"),
    censored = sum(records$censored))
  structure(list(summary = summary, records = records, feedback = feedback_summary,
                 paths = if (keep_paths) paths else NULL, model = model),
            class = "stemfeedback_result")
}

print.stemfeedback_result <- function(x, ...) {
  s <- x$summary
  cat("Stem-cell competition:", s$attempts, "attempts\n")
  if (s$unresolved > 0) {
    cat("Unresolved outcomes:", s$unresolved, "\n")
    cat("Leukemic fixation bounds:",
        format(c(s$outcome_lower_bound, s$outcome_upper_bound)), "\n")
  } else {
    cat("Leukemic fixation:", format(s$estimate), "\n")
    cat("95% Wilson interval:", format(c(s$lower, s$upper)), "\n")
  }
  invisible(x)
}

run_experiments <- function(output, suite = "all", profile = "smoke",
                            seed = 20261028L, workers = 1L, replicates = NULL,
                            max_events = 2e7, force = FALSE, plots = TRUE) {
  suite <- match.arg(suite, c("all", "baseline", "weak_selection", "clearance", "transit",
    "mechanism", "founders", "initialization", "entry", "staged", "noise", "occupation",
    "establishment", "latency", "deterministic", "audits"))
  profile <- match.arg(profile, c("smoke", "quick", "full"))
  stopifnot(length(output) == 1L, !is.na(output), nzchar(output),
            length(seed) == 1L, is.finite(seed), seed >= 1, seed == round(seed),
            seed + 500000000 < .Machine$integer.max,
            length(workers) == 1L, is.finite(workers), workers >= 1, workers == round(workers),
            length(max_events) == 1L, max_events >= 1)
  if (!is.null(replicates))
    stopifnot(length(replicates) == 1L, is.finite(replicates), replicates >= 1,
              replicates == round(replicates), replicates < 250000)
  seed <- as.integer(seed); workers <- as.integer(workers)
  if (!is.null(replicates)) replicates <- as.integer(replicates)
  output <- path.expand(output)
  preflight_output(output, suite, profile, seed, max_events, force)
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  output <- normalizePath(output, mustWork = TRUE)
  runtime <- runtime_record()
  write_csv(data.frame(date_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    R_version = runtime$R_version, Matrix_version = runtime$Matrix_version,
    deSolve_version = runtime$deSolve_version, code_fingerprint = code_fingerprint(),
    package_version = as.character(utils::packageVersion("stemfeedback")),
    profile = profile, suite = suite, base_seed = seed, workers = workers,
    replicates_override = if (is.null(replicates)) NA else replicates,
    max_events = max_events), file.path(output, "run_manifest.csv"))
  writeLines(capture.output(sessionInfo()), file.path(output, "session_info.txt"))
  if (suite %in% c("all", "deterministic")) run_deterministic(output, profile)
  if (suite %in% c("all", "audits")) run_audits(output, profile, seed)
  if (!suite %in% c("deterministic", "audits"))
    run_simulations(output, suite, profile, seed, workers, replicates, force, max_events)
  if (plots) plot_results(output, profile)
  writeLines(capture.output(sessionInfo()), file.path(output, "session_info.txt"))
  invisible(output)
}
