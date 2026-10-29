parallel_paths <- function(indices, fun, workers = 1L) {
  if (workers == 1L) return(lapply(indices, fun))
  if (.Platform$OS.type != "windows")
    return(parallel::mclapply(indices, fun, mc.cores = workers, mc.set.seed = FALSE))
  cluster <- parallel::makePSOCKcluster(workers); on.exit(parallel::stopCluster(cluster))
  parallel::clusterCall(cluster, function(libs) {
    .libPaths(libs); loadNamespace("stemfeedback"); NULL
  }, .libPaths())
  parallel::parLapply(cluster, indices, fun)
}

path_seed <- function(base_seed, case_index, attempt, continuation = 0L) {
  offset <- if (continuation > 0) 500000 + 2 * (attempt - 1) + continuation else attempt
  seed <- base_seed + case_index * 1000000 + offset
  if (seed >= .Machine$integer.max || offset >= 1000000) stop("Seed allocation exhausted.")
  as.integer(seed)
}

code_fingerprint <- function() {
  file <- system.file("SOURCE_SHA256.txt", package = "stemfeedback")
  if (!nzchar(file)) stop("The installed package source fingerprint is missing.")
  paste(as.character(utils::packageVersion("stemfeedback")),
        paste(readLines(file, warn = FALSE), collapse = ":"), sep = ":")
}

runtime_record <- function() {
  version <- function(name) tryCatch(as.character(utils::packageVersion(name)),
                                    error = function(e) NA_character_)
  list(R_version = R.version.string, Matrix_version = version("Matrix"),
       deSolve_version = version("deSolve"))
}

preflight_output <- function(output, suite, profile, seed, max_events, force = FALSE) {
  if (!dir.exists(output)) return(invisible(TRUE))
  files <- list.files(output, all.files = TRUE, no.. = TRUE)
  if (!length(files)) return(invisible(TRUE))
  manifest_file <- file.path(output, "run_manifest.csv")
  if (!file.exists(manifest_file))
    stop("The output folder is nonempty and has no run manifest. Choose a new folder.")
  manifest <- read.csv(manifest_file, stringsAsFactors = FALSE, colClasses = "character")
  expected <- c(list(suite = suite, profile = profile, base_seed = seed,
                    max_events = max_events, code_fingerprint = code_fingerprint()), runtime_record())
  required <- c("suite", "profile", "base_seed", "max_events", "code_fingerprint", "R_version")
  if (nrow(manifest) != 1L || !all(required %in% names(manifest)))
    stop("The existing run manifest is incomplete. Choose a new output folder.")
  mismatch <- names(expected)[!vapply(names(expected), function(name) {
    if (!name %in% names(manifest)) return(FALSE)
    actual <- manifest[[name]][1]
    if (is.numeric(expected[[name]])) actual <- as.numeric(actual)
    isTRUE(all.equal(actual, expected[[name]], check.attributes = FALSE))
  }, FALSE)]
  if (length(mismatch) && !force)
    stop("Existing run differs in ", paste(mismatch, collapse = ", "),
         ". No files were changed. Use a new folder or force = TRUE to replace the generated run.")
  if (force) {
    managed <- c("design", "paths", "summary", "theory", "deterministic", "audits",
                 "figures", "manuscript_figures", "continuations", "staged_scores", "run_manifest.csv", "session_info.txt")
    unlink(file.path(output, managed), recursive = TRUE)
  }
  invisible(TRUE)
}

case_signature <- function(case, base_seed, max_events) {
  f <- tempfile(); on.exit(unlink(f))
  writeBin(serialize(list(case = case, seed = base_seed, max_events = max_events,
                          code = code_fingerprint()), NULL, version = 2), f)
  unname(tools::md5sum(f))
}

case_observations <- function(case) {
  if (case$suite == "occupation") return(c(0, 2, 5, 10, 20, 40, 80, 160))
  if (case$suite %in% c("noise", "mechanism")) return(seq(0, 2, .25))
  0
}

run_case <- function(case, attempts, output, base_seed = 20261028L, workers = 1L,
                     force = FALSE, max_events = 2e7, reference = NULL,
                     enlarged_reference = NULL) {
  file <- file.path(output, "paths", paste0(case$id, ".csv"))
  signature <- case_signature(case, base_seed, max_events)
  old <- data.frame()
  if (file.exists(file) && !force) {
    old <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
    if (!nrow(old) || !all(old$signature == signature) || anyDuplicated(old$attempt))
      stop("Cached paths differ from the current implementation: ", file,
           ". Choose another output directory or use force = TRUE.")
    expected <- vapply(old$attempt, function(i) path_seed(base_seed, case$index, i), 0L)
    if (!all(old$seed == expected)) stop("Cached seeds do not match the design.")
  } else if (file.exists(file)) unlink(file)
  missing <- setdiff(seq_len(attempts), old$attempt)
  if (length(missing)) {
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    batches <- split(missing, ceiling(seq_along(missing) / 100L))
    for (batch in batches) {
      calculate <- function(i) {
        path <- simulate_path(case$initial, case$model, path_seed(base_seed, case$index, i),
          stop_rule = case$stop_rule, threshold = case$threshold, burden = case$burden,
          max_events = max_events, observation_times = case_observations(case),
          reference = reference, enlarged_reference = enlarged_reference)
        row <- path_record(path, i, case$id); row$signature <- signature; row
      }
      rows <- parallel_paths(batch, calculate, workers)
      if (any(vapply(rows, inherits, FALSE, what = "try-error"))) stop("A simulation worker failed.")
      rows <- bind_rows(rows)
      write.table(rows, file, sep = ",", row.names = FALSE, col.names = !file.exists(file),
                  append = file.exists(file), quote = TRUE, na = "")
    }
  }
  ans <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  ans <- ans[ans$attempt <= attempts, , drop = FALSE]
  ans[order(ans$attempt), , drop = FALSE]
}

summarize_case <- function(records, case, profile) {
  n <- nrow(records); censored <- sum(records$censored)
  unresolved <- censored + sum(records$status == "reference_exit")
  target <- if (case$stop_rule %in% c("establishment", "burden")) records$hit else records$healthy_loss
  w <- if (!unresolved) wilson(sum(target), n) else c(estimate = NA, lower = NA, upper = NA)
  m <- case$model; eq <- equilibrium(m)
  prediction <- if (case$stop_rule == "first_stem_loss") diffusion_probability(case$initial, m) else NA
  out <- data.frame(case_id = case$id, suite = case$suite, profile = profile,
    N = m$N, p = m$p, q = m$q, a = m$a, b = m$b, k = m$k,
    dx = m$healthy$clearance, dy = m$leukemic$clearance, noise = m$noise,
    K = eq$K, r = eq$r, attempts = n, censored = censored, unresolved = unresolved,
    healthy_losses = sum(records$healthy_loss), hits = sum(records$hit),
    estimate = unname(w[1]), lower = unname(w[2]), upper = unname(w[3]),
    outcome_lower_bound = sum(target) / n, outcome_upper_bound = (sum(target) + unresolved) / n,
    diffusion_prediction = prediction, mean_events = mean(records$events),
    mean_stopping_time = mean(records$time), note = case$note)
  if (case$suite %in% c("founders", "initialization")) {
    out$scaled_estimate <- m$N * out$estimate; out$scaled_lower <- m$N * out$lower
    out$scaled_upper <- m$N * out$upper
    out$founder_coefficient <- scale_function(m)$founder_coefficient
  }
  if (case$stop_rule %in% c("establishment", "burden")) {
    bp <- branching_prediction(m, if (is.null(case$threshold)) Inf else case$threshold)
    out$branching_survival <- bp["survival"]; out$branching_threshold <- bp["threshold_probability"]
    out$lambda <- bp["lambda"]
  }
  out
}

staged_founder <- function(case, entry_records, output, base_seed, workers = 1L,
                           max_events = 2e7, bootstrap = 2000L, force = FALSE) {
  if (any(entry_records$censored)) stop("Resolve all entry attempts before staged estimation.")
  hits <- which(entry_records$status == "threshold_hit")
  jobs <- if (length(hits)) expand.grid(row = hits, continuation = 1:2) else data.frame()
  continuations <- data.frame(case_id = character(), original_attempt = integer(),
    continuation = integer(), seed = integer(), healthy_loss = integer(),
    censored = logical(), status = character(), time = numeric(), events = integer())
  if (nrow(jobs)) {
    calculate <- function(j) {
      hit <- entry_records[jobs$row[j], ]; initial <- as.numeric(unlist(
        hit[paste0("final_", case$model$names)], use.names = FALSE))
      seed <- path_seed(base_seed, case$index, hit$attempt, jobs$continuation[j])
      p <- simulate_path(initial, case$model, seed, max_events = max_events, observation_times = 0)
      data.frame(case_id = case$id, original_attempt = hit$attempt,
                 continuation = jobs$continuation[j], seed = seed, healthy_loss = p$healthy_loss,
                 censored = p$censored, status = p$status, time = p$time, events = p$events)
    }
    continuations <- bind_rows(parallel_paths(seq_len(nrow(jobs)), calculate, workers))
    if (any(continuations$censored)) stop("A continuation was censored; no staged probability was reported.")
  }
  scores <- as.numeric(entry_records$status == "healthy_loss")
  for (i in hits) scores[i] <- mean(continuations$healthy_loss[
    continuations$original_attempt == entry_records$attempt[i]])

  ci <- cluster_score_interval(scores, bootstrap,
                               seed = base_seed + 300000000 + case$index)
  write_csv(continuations, file.path(output, "continuations", paste0(case$id, ".csv")))
  write_csv(data.frame(case_id = case$id, attempt = entry_records$attempt, score = scores),
            file.path(output, "staged_scores", paste0(case$id, ".csv")))
  data.frame(case_id = case$id, N = case$model$N, r = equilibrium(case$model)$r,
    original_attempts = length(scores), threshold_hits = length(hits),
    continuations = nrow(continuations), continuation_fixations = sum(continuations$healthy_loss),
    scaled_estimate = case$model$N * ci$estimate, scaled_lower = case$model$N * ci$lower,
    scaled_upper = case$model$N * ci$upper,
    founder_coefficient = scale_function(case$model)$founder_coefficient,
    bootstrap_resamples = ci$bootstrap_resamples, bootstrap_requested = bootstrap,
    interval_method = ci$method)
}

summarize_latency <- function(records, case, shift, bootstrap = 1500L, seed = 1L) {
  bp <- branching_prediction(case$model)
  lambda <- unname(bp["lambda"]); survival <- unname(bp["survival"])
  success <- records[records$status == "burden_hit", , drop = FALSE]
  centered <- success$time - log(case$model$N) / lambda
  complete <- !any(records$censored)
  values <- if (complete) centered else numeric()
  ci <- mean_interval(values)
  mean_ci <- if (length(values) >= 2L)
    bootstrap_interval(values, mean, B = bootstrap, seed = seed) else
      c(estimate = NA_real_, lower = NA_real_, upper = NA_real_)
  variance_ci <- if (length(values) >= 2L)
    bootstrap_interval(values, var, B = bootstrap, seed = seed + 1L) else
      c(estimate = NA_real_, lower = NA_real_, upper = NA_real_)
  moments <- data.frame(case_id = case$id, N = case$model$N,
    attempts = nrow(records), hits = nrow(success), censored = sum(records$censored),
    complete_run = complete,
    moment_status = if (!complete) "incomplete_censored" else if (!length(values)) "no_hits" else "complete",
    mean_crossing_time = if (complete && nrow(success)) mean(success$time) else NA_real_,
    empirical_mean = ci["mean"], mean_se = ci["se"],
    mean_lower = mean_ci["lower"], mean_upper = mean_ci["upper"],
    mean_bootstrap_resamples = if (length(values) >= 2L) bootstrap else 0L,
    mean_interval_method = "percentile_bootstrap",
    empirical_variance = variance_ci["estimate"],
    variance_lower = variance_ci["lower"], variance_upper = variance_ci["upper"],
    variance_bootstrap_resamples = if (length(values) >= 2L) bootstrap else 0L,
    empirical_skewness = sample_skewness(values),
    predicted_conditional_mean = shift + (0.5772156649015329 + log(survival)) / lambda,
    predicted_conditional_variance = pi^2 / (6 * lambda^2),
    predicted_conditional_skewness = 12 * sqrt(6) * 1.202056903159594 / pi^3,
    deterministic_shift = shift, branching_survival = survival)
  probabilities <- c(.1, .25, .5, .75, .9)
  quantiles <- data.frame(case_id = case$id, probability = probabilities,
    empirical = if (length(values)) unname(quantile(values, probabilities)) else NA_real_,
    predicted = shift + (log(survival) - log(-log(probabilities))) / lambda,
    complete_run = complete)
  paths <- data.frame(case_id = rep(case$id, nrow(success)), attempt = success$attempt,
    seed = success$seed, crossing_time = success$time, centered_time = centered)
  list(moments = moments, quantiles = quantiles, paths = paths)
}

summarize_entry <- function(records, case) {
  m <- case$threshold; eq <- equilibrium(case$model); N <- case$model$N
  if (any(records$censored)) stop("Resolve threshold runs before calculating entry statistics.")
  hit <- records$status == "threshold_hit"; R <- records$final_X1[hit] / (eq$K * N)
  weighted <- rep(0, nrow(records)); weighted[hit] <- m * R^(-case$model$q / case$model$p)
  ci <- wilson(sum(hit), nrow(records)); h <- mean_interval(weighted)
  data.frame(case_id = case$id, N = N, r = eq$r, threshold = m,
    attempts = nrow(records), hits = sum(hit), m_hit_probability = m * ci[1],
    m_hit_lower = m * ci[2], m_hit_upper = m * ci[3],
    mean_resident_ratio_at_hit = if (length(R)) mean(R) else NA,
    weighted_entry = h["mean"], weighted_entry_se = h["se"])
}

summarize_fluctuations <- function(records, case) {
  survivors <- records[is.finite(records$rms), , drop = FALSE]
  rms <- mean_interval(survivors$rms)
  times <- seq(0, 2, .25)
  tags <- vapply(times, function(t) gsub("[.]", "p", format(t, trim = TRUE, scientific = FALSE)), "")
  columns <- paste0("t", tags, "_signal")
  lag <- NA_real_
  if (nrow(survivors) >= 2 && all(columns %in% names(survivors))) {
    values <- as.matrix(survivors[columns]); centred <- sweep(values, 2, colMeans(values), "-")
    left <- as.vector(centred[, -ncol(centred), drop = FALSE])
    right <- as.vector(centred[, -1, drop = FALSE])
    if (sd(left) > 0 && sd(right) > 0) lag <- cor(left, right)
  }
  data.frame(case_id = case$id, suite = case$suite, N = case$model$N,
    speed = if (case$model$n == 3) case$model$healthy$division[2] else case$model$healthy$clearance,
    noise = case$model$noise, attempts = nrow(records), observed_at_two = nrow(survivors),
    signal_rms = rms["mean"], rms_se = rms["se"], rms_lower = rms["lower"], rms_upper = rms["upper"],
    lag_025 = lag, transverse_sd = if (nrow(survivors) > 1) sd(survivors$transverse) else NA)
}

summarize_occupation <- function(records, case, reference) {
  if (any(records$censored)) stop("Resolve occupation paths before reporting the complete-path identity.")
  I <- mean_interval(records$integral); D <- mean_interval(records$residual)
  w <- wilson(sum(records$healthy_loss), nrow(records)); U0 <- records$reference_initial[1]
  summary <- data.frame(case_id = case$id, noise = case$model$noise, attempts = nrow(records),
    reference = U0, mean_integral = I["mean"], integral_se = I["se"],
    integral_lower = I["lower"], integral_upper = I["upper"],
    stopped_healthy_loss_fraction = mean(records$healthy_loss),
    fixation_lower_bound = mean(records$healthy_loss),
    fixation_upper_bound = mean(records$healthy_loss) + mean(records$status == "reference_exit"),
    fixation = if (any(records$status == "reference_exit")) NA_real_ else w[1],
    fixation_lower = if (any(records$status == "reference_exit")) NA_real_ else w[2],
    fixation_upper = if (any(records$status == "reference_exit")) NA_real_ else w[3],
    mean_residual = D["mean"],
    residual_se = D["se"], residual_lower = D["lower"], residual_upper = D["upper"],
    outer_exits = sum(records$status == "reference_exit"),
    reference_escape_bound = reference$E[case$initial[1] + 1L, case$initial[case$model$stem[2]] + 1L],
    max_path_error_bound = max(records$integration_error_bound),
    max_enlargement_change = max(abs(records$integral - records$enlarged_integral)))
  bands <- c(0, .1, .25, .5, .75, .9, 1); frequency <- list()
  for (b in 1:6) {
    tau <- records[[paste0("band", b, "_time")]]
    ds <- records[[paste0("band", b, "_feedback")]]
    psi <- records[[paste0("band", b, "_response")]]
    value <- records[[paste0("band", b, "_integral")]]
    means <- c(mean(tau), mean(ds), mean(psi), mean(value)); T <- means[1]
    product <- if (T > 0) means[2] * means[3] / T else NA_real_
    gradient <- if (T > 0) c(-means[2] * means[3] / T^2, means[3] / T, means[2] / T, 0) else rep(NA_real_, 4)
    covariance <- if (length(tau) > 1) stats::cov(cbind(tau, ds, psi, value)) / length(tau) else matrix(NA_real_, 4, 4)
    g_cov <- c(0, 0, 0, 1) - gradient
    frequency[[b]] <- data.frame(case_id = case$id, noise = case$model$noise,
      band_lower = bands[b], band_upper = bands[b + 1L], mean_time = T,
      mean_feedback_departure = if (T > 0) means[2] / T else NA,
      integral = means[4], integral_se = if (length(value) > 1) sd(value) / sqrt(length(value)) else NA,
      product_of_means = product, covariance_contribution = means[4] - product,
      product_se = sqrt(max(0, drop(t(gradient) %*% covariance %*% gradient))),
      covariance_se = sqrt(max(0, drop(t(g_cov) %*% covariance %*% g_cov))))
  }
  time <- list()
  for (t in c(0, 2, 5, 10, 20, 40, 80, 160)) {
    tag <- gsub("[.]", "p", format(t, trim = TRUE, scientific = FALSE))
    i <- mean_interval(records[[paste0("t", tag, "_integral")]])
    u <- mean_interval(records[[paste0("t", tag, "_reference")]])
    residual <- records[[paste0("t", tag, "_reference")]] - U0 - records[[paste0("t", tag, "_integral")]]
    d <- mean_interval(residual)
    time[[length(time) + 1L]] <- data.frame(case_id = case$id, noise = case$model$noise, time = t,
      integral = i["mean"], integral_se = i["se"], reference_mean = u["mean"],
      residual = d["mean"], residual_se = d["se"],
      cumulative_healthy_loss = mean(records$healthy_loss == 1 & records$time <= t),
      cumulative_healthy_loss_lower = wilson(sum(records$healthy_loss == 1 & records$time <= t), nrow(records))[2],
      cumulative_healthy_loss_upper = wilson(sum(records$healthy_loss == 1 & records$time <= t), nrow(records))[3],
      cumulative_stem_loss = mean(records$status %in% c("healthy_loss", "leukemic_loss") & records$time <= t))
  }
  list(summary = summary, frequency = bind_rows(frequency), time = bind_rows(time),
       decomposition = occupation_decomposition(records, case))
}

occupation_decomposition <- function(records, case) {
  n <- nrow(records); product <- 0; influence <- numeric(n)
  for (b in 1:6) {
    tau <- records[[paste0("band", b, "_time")]]
    ds <- records[[paste0("band", b, "_feedback")]]
    psi <- records[[paste0("band", b, "_response")]]
    T <- mean(tau); D <- mean(ds); P <- mean(psi)
    if (T > 0) {
      product <- product + D * P / T
      influence <- influence - D * P / T^2 * (tau - T) +
        P / T * (ds - D) + D / T * (psi - P)
    }
  }
  integral <- records$integral; correction <- mean(integral)
  data.frame(case_id = case$id, noise = case$model$noise, attempts = n,
    product_of_means = product, covariance_contribution = correction - product,
    integral = correction,
    product_se = if (n > 1) sd(influence) / sqrt(n) else NA_real_,
    covariance_se = if (n > 1) sd(integral - correction - influence) / sqrt(n) else NA_real_,
    integral_se = if (n > 1) sd(integral) / sqrt(n) else NA_real_)
}
