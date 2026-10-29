reference_key <- function(model, limits) {
  eq <- equilibrium(model)
  paste(format(c(model$N, model$p, model$q, model$a, model$b, model$k,
                 eq$A, eq$B, model$asymmetric, limits), digits = 16), collapse = ":")
}

reference_for <- function(model, limits, cache) {
  key <- reference_key(model, limits)
  if (!exists(key, cache, inherits = FALSE)) assign(key, solve_reference(model, limits), cache)
  get(key, cache, inherits = FALSE)
}

reference_summary <- function(case, reference, enlarged = NULL) {
  m <- case$model; i <- case$initial[1] + 1L; j <- case$initial[m$stem[2]] + 1L
  data.frame(case_id = case$id, N = m$N, q = m$q, noise = m$noise,
    lower_bound = reference$U[i, j], escape_probability = reference$E[i, j],
    upper_bound = reference$U[i, j] + reference$E[i, j],
    enlarged_lower_bound = if (is.null(enlarged)) NA_real_ else enlarged$U[i, j],
    linear_solve_residual = reference$residual,
    outer_X1 = reference$limits[1], outer_Y1 = reference$limits[2])
}

selection_diagnostics <- function(cases) {
  bind_rows(lapply(cases, function(c) {
    m <- c$model; eq <- equilibrium(m)
    x <- c$initial[1] / m$N; y <- c$initial[m$stem[2]] / m$N
    projected <- if (x + y > 0) project_stems(x, y, m) else NA_real_
    data.frame(case_id = c$id, projected_healthy_stem_density = projected,
      projected_leukemic_stem_density = eq$r * (eq$K - projected),
      diffusion_fixation = if (c$stop_rule == "first_stem_loss")
        scale_function(m)$probability(projected) else NA_real_,
      founder_coefficient = if (c$suite %in% c("founders", "entry", "initialization"))
        scale_function(m)$founder_coefficient else NA_real_)
  }))
}

mechanism_contrast <- function(records_a, records_b, B = 2000L, seed = 1L) {
  summary_values <- function(x) {
    x <- x[is.finite(x$rms), , drop = FALSE]
    if (nrow(x) < 2) return(c(rms = NA, lag = NA, transverse = NA))
    tags <- vapply(seq(0, 2, .25), function(t)
      gsub("[.]", "p", format(t, trim = TRUE, scientific = FALSE)), "")
    columns <- paste0("t", tags, "_signal")
    values <- as.matrix(x[columns]); centred <- sweep(values, 2, colMeans(values), "-")
    left <- as.vector(centred[, -ncol(centred), drop = FALSE])
    right <- as.vector(centred[, -1, drop = FALSE])
    c(rms = mean(x$rms),
      lag = if (sd(left) > 0 && sd(right) > 0) cor(left, right) else NA,
      transverse = sd(x$transverse))
  }
  observed <- summary_values(records_b) - summary_values(records_a)
  set.seed(seed)
  draws <- replicate(B, {
    a <- records_a[sample.int(nrow(records_a), replace = TRUE), , drop = FALSE]
    b <- records_b[sample.int(nrow(records_b), replace = TRUE), , drop = FALSE]
    summary_values(b) - summary_values(a)
  })
  bind_rows(lapply(seq_along(observed), function(j) {
    v <- draws[j, ]; v <- v[is.finite(v)]
    data.frame(metric = names(observed)[j], difference = observed[j],
      lower = if (length(v)) unname(quantile(v, .025)) else NA_real_,
      upper = if (length(v)) unname(quantile(v, .975)) else NA_real_,
      bootstrap_valid = length(v), bootstrap_resamples = B)
  }))
}

run_simulations <- function(output, suite = "all", profile = "smoke", seed = 20261028L,
                            workers = 1L, replicates = NULL, force = FALSE,
                            max_events = 2e7) {
  all_cases <- paper_designs()
  selected <- select_designs(all_cases, if (suite == "staged") "staged" else suite, profile)
  if (!length(selected)) stop("No designs selected for suite: ", suite)
  write_csv(design_table(all_cases, "full"), file.path(output, "design", "paper_designs.csv"))
  write_csv(design_table(selected, profile, replicates), file.path(output, "design", "selected_designs.csv"))
  write_csv(selection_diagnostics(selected), file.path(output, "theory", "selection.csv"))
  cache <- new.env(parent = emptyenv()); records <- list(); summaries <- list()
  refrows <- list(); occup <- list(); entries <- list(); fluctuations <- list()
  for (case in selected) {
    attempts <- case_attempts(case, profile, replicates)
    message("Running ", case$id, ": ", attempts, " original attempts")
    ref <- enlarged <- NULL
    if (case$suite == "occupation") {
      ref <- reference_for(case$model, c(300L, 300L), cache)
      enlarged <- reference_for(case$model, c(400L, 400L), cache)
      refrows[[case$id]] <- reference_summary(case, ref, enlarged)
    }
    records[[case$id]] <- run_case(case, attempts, output, seed, workers, force,
                                    max_events, ref, enlarged)
    summaries[[case$id]] <- summarize_case(records[[case$id]], case, profile)
    if (case$suite == "entry" && !any(records[[case$id]]$censored))
      entries[[case$id]] <- summarize_entry(records[[case$id]], case)
    if (case$suite %in% c("mechanism", "noise"))
      fluctuations[[case$id]] <- summarize_fluctuations(records[[case$id]], case)
    if (case$suite == "occupation" && !any(records[[case$id]]$censored))
      occup[[case$id]] <- summarize_occupation(records[[case$id]], case, ref)
  }
  all_summary <- bind_rows(summaries)
  write_csv(all_summary, file.path(output, "summary", "all_cases.csv"))
  for (name in unique(all_summary$suite))
    write_csv(all_summary[all_summary$suite == name, , drop = FALSE],
              file.path(output, "summary", paste0(name, ".csv")))
  if (length(refrows)) write_csv(bind_rows(refrows), file.path(output, "theory", "reference_bounds.csv"))
  if (length(entries)) write_csv(bind_rows(entries), file.path(output, "summary", "entry_statistics.csv"))
  if (length(fluctuations))
    write_csv(bind_rows(fluctuations), file.path(output, "summary", "feedback_fluctuations.csv"))
  if (length(occup)) {
    write_csv(bind_rows(lapply(occup, `[[`, "summary")), file.path(output, "summary", "occupation.csv"))
    write_csv(bind_rows(lapply(occup, `[[`, "frequency")), file.path(output, "summary", "occupation_bands.csv"))
    write_csv(bind_rows(lapply(occup, `[[`, "time")), file.path(output, "summary", "occupation_time.csv"))
    write_csv(bind_rows(lapply(occup, `[[`, "decomposition")),
              file.path(output, "summary", "occupation_decomposition.csv"))
    for (case in selected) if (case$suite == "occupation") {
      x <- records[[case$id]]
      write_csv(data.frame(case_id = case$id, attempt = x$attempt,
        integral = x$integral, enlarged_integral = x$enlarged_integral,
        residual = x$residual, error_bound = x$integration_error_bound,
        departure_positive = x$positive_departure_integral,
        departure_negative = x$negative_departure_integral),
        file.path(output, "summary", paste0(case$id, "_path_integrals.csv")))
    }
  }
  if (suite %in% c("all", "staged", "entry")) {
    staged_cases <- selected[vapply(selected, function(c) c$suite == "entry" &&
        c$model$N %in% c(800, 1600) && c$threshold == c$model$N / 40, FALSE)]
    if (length(staged_cases) && !any(vapply(staged_cases, function(c)
        any(records[[c$id]]$censored), FALSE))) {
      stage <- bind_rows(lapply(staged_cases, function(c) staged_founder(c, records[[c$id]],
        output, seed, workers, max_events, bootstrap = 2000L, force = force)))
      write_csv(stage, file.path(output, "summary", "staged_founder.csv"))
    }
  }
  if (length(fluctuations)) {
    mechanism_cases <- selected[vapply(selected, function(c) c$suite == "mechanism", FALSE)]
    contrasts <- list()
    for (family in c("clearance", "transit")) for (N in c(50, 200)) {
      pair <- mechanism_cases[vapply(mechanism_cases, function(c) {
        speed <- if (c$model$n == 3) c$model$healthy$division[2] else c$model$healthy$clearance
        c$model$N == N && startsWith(c$id, paste0("mechanism_", family)) && speed %in% c(.25, 8)
      }, FALSE)]
      if (length(pair) == 2) {
        speeds <- vapply(pair, function(c) if (c$model$n == 3)
          c$model$healthy$division[2] else c$model$healthy$clearance, 0.)
        pair <- pair[order(speeds)]
        z <- mechanism_contrast(records[[pair[[1]]$id]], records[[pair[[2]]$id]],
          B = if (profile == "smoke") 200L else 2000L,
          seed = seed + 200000000L + pair[[1]]$index)
        z$family <- family; z$N <- N; z$slow_case <- pair[[1]]$id; z$fast_case <- pair[[2]]$id
        contrasts[[length(contrasts) + 1L]] <- z
      }
    }
    if (length(contrasts)) write_csv(bind_rows(contrasts),
      file.path(output, "summary", "mechanism_contrasts.csv"))
  }
  if (length(occup) == 2) {
    a <- records[["occupation_0"]]; b <- records[["occupation_16"]]
    if (!is.null(a) && !is.null(b)) {
      d <- independent_contrast(a$integral, b$integral)
      write_csv(data.frame(noise_low = 0, noise_high = 16, difference = d["difference"],
        se = d["se"], lower = d["lower"], upper = d["upper"],
        attempts_low = nrow(a), attempts_high = nrow(b)),
        file.path(output, "summary", "occupation_contrast.csv"))
    }
  }
  noise_cases <- selected[vapply(selected, function(c) c$suite == "noise", FALSE)]
  noise_contrasts <- list()
  for (N in unique(vapply(noise_cases, function(c) c$model$N, 0.))) {
    pair <- noise_cases[vapply(noise_cases, function(c) c$model$N == N &&
      c$model$noise %in% c(0, 4), FALSE)]
    if (length(pair) == 2) {
      pair <- pair[order(vapply(pair, function(c) c$model$noise, 0.))]
      a <- records[[pair[[1]]$id]]; b <- records[[pair[[2]]$id]]
      if (!any(a$censored) && !any(b$censored)) {
        difference <- mean(b$healthy_loss) - mean(a$healthy_loss)
        se <- sqrt(var(a$healthy_loss) / nrow(a) + var(b$healthy_loss) / nrow(b))
        noise_contrasts[[length(noise_contrasts) + 1L]] <- data.frame(N = N,
          low_noise = 0, high_noise = 4, fixation_difference = difference,
          se = se, lower = difference - 1.96 * se, upper = difference + 1.96 * se,
          low_attempts = nrow(a), high_attempts = nrow(b))
      }
    }
  }
  if (length(noise_contrasts)) write_csv(bind_rows(noise_contrasts),
    file.path(output, "summary", "terminal_noise_contrasts.csv"))
  latency_cases <- selected[vapply(selected, function(c) c$suite == "latency", FALSE)]
  if (length(latency_cases)) {
    shift <- tail(latency_shift(minimal_model(2000, q = 1, b = 1, dy = 1))$deterministic_shift, 1)
    by_path <- list(); by_case <- list(); quantiles <- list()
    for (case in latency_cases) {
      report <- summarize_latency(records[[case$id]], case, shift, bootstrap = 1500L,
                                  seed = seed + 450000000L + case$index)
      by_path[[case$id]] <- report$paths
      by_case[[case$id]] <- report$moments
      quantiles[[case$id]] <- report$quantiles
    }
    path_times <- bind_rows(by_path)
    if (!nrow(path_times)) path_times <- data.frame(case_id = character(),
      attempt = integer(), seed = integer(), crossing_time = numeric(), centered_time = numeric())
    write_csv(path_times, file.path(output, "summary", "latency_hit_times.csv"))
    write_csv(bind_rows(by_case), file.path(output, "summary", "latency_moments.csv"))
    write_csv(bind_rows(quantiles), file.path(output, "summary", "latency_quantiles.csv"))
  }
  invisible(list(cases = selected, records = records, summary = all_summary))
}

run_deterministic <- function(output, profile = "smoke", trajectory_designs = NULL) {
  require_ode()
  designs <- if (is.null(trajectory_designs)) manuscript_trajectory_designs() else trajectory_designs
  validate_trajectory_designs(designs)
  illustrations <- list(); illustrated <- list(); settings <- list(); segments <- list()
  for (id in names(designs)) {
    design <- designs[[id]]; zmodel <- design$model; family <- design$family
    z_eq <- equilibrium(zmodel); start <- design$initial
    expected_x <- project_stems(start[1], start[zmodel$stem[2]], zmodel)
    v <- ode_trajectory(start, zmodel, seq(0, 800, by = if (profile == "smoke") 4 else 1))
    v$family <- family; v$trajectory_id <- id; v$provenance <- design$provenance
    illustrated[[id]] <- v
    last <- v[nrow(v), ]
    illustrations[[id]] <- data.frame(family = family, trajectory_id = id,
      initial_healthy_stem_density = start[1],
      initial_leukemic_stem_density = start[zmodel$stem[2]],
      projected_healthy_stem_density = expected_x,
      recovered_healthy_stem_density = last$X1,
      projected_leukemic_stem_density = z_eq$r * (z_eq$K - expected_x),
      recovered_leukemic_stem_density = last$Y1)
    settings[[id]] <- data.frame(trajectory_id = id, family = family,
      compartment = zmodel$names, initial_density = start,
      feedback_weight = zmodel$weights, provenance = design$provenance)
    xi <- seq(0, z_eq$K, length.out = 301)
    segments[[id]] <- data.frame(trajectory_id = id, family = family,
      X1 = xi, Y1 = z_eq$r * (z_eq$K - xi))
  }
  dir.create(file.path(output, "deterministic"), recursive = TRUE, showWarnings = FALSE)
  saveRDS(designs, file.path(output, "deterministic", "trajectory_designs.rds"))
  write_csv(bind_rows(settings), file.path(output, "deterministic", "trajectory_settings.csv"))
  write_csv(bind_rows(segments), file.path(output, "deterministic", "coexistence_segments.csv"))
  write_csv(bind_rows(illustrations), file.path(output, "deterministic", "coexistence_recovery.csv"))
  write_csv(bind_rows(illustrated), file.path(output, "deterministic", "coexistence_trajectories.csv"))
  m <- minimal_model(100, q = 3, dy = 1.5)
  initial <- c(X1 = .25, X2 = .25, Y1 = .125, Y2 = .25)
  schedules <- data.frame(name = c("equal_independent", "common_division", "selective_1", "selective_2"),
                          mux = c(1, 0, .5, .5), muy = c(1, 0, 1.5, 2.5),
                          kappa = c(0, .6, 0, 0))
  therapy <- list(); trajectory <- list()
  for (j in seq_len(nrow(schedules))) {
    z <- schedules[j, ]; mu <- rep(0, m$dimension)
    mu[m$stem] <- c(z$mux, z$muy)
    out <- ode_trajectory(initial, m, seq(0, 100, by = if (profile == "smoke") .5 else .1),
      treatment = list(start = 0, end = 3, mu = mu, kappa = z$kappa))
    out$schedule <- z$name; trajectory[[j]] <- out
    at3 <- out[which.min(abs(out$time - 3)), ]; last <- out[nrow(out), ]
    prediction <- recovery_prediction(initial[m$stem], m,
      exposure = 3 * c(z$mux, z$muy), after = c(at3$X1, at3$Y1), counts = FALSE)
    therapy[[j]] <- data.frame(schedule = z$name, mu_x = z$mux, mu_y = z$muy,
      kappa = z$kappa, descendant_killing = 0,
      predicted_delta_h = prediction$predicted_delta_h,
      measured_delta_h = prediction$observed_delta_h,
      delta_h_error = prediction$discrepancy,
      predicted_recovered_leukemic_fraction = prediction$leukemic_fraction,
      observed_fraction_at_100 = last$Y1 / (last$X1 + last$Y1),
      recovered_healthy_stems = last$X1, recovered_leukemic_stems = last$Y1)
  }
  write_csv(bind_rows(therapy), file.path(output, "deterministic", "therapy_summary.csv"))
  write_csv(schedules, file.path(output, "deterministic", "therapy_settings.csv"))
  write_csv(bind_rows(trajectory), file.path(output, "deterministic", "therapy_trajectories.csv"))
  lmodel <- minimal_model(2000, q = 1, b = 1, dy = 1)
  latency <- latency_shift(lmodel, densities = 10^-(4:7))
  write_csv(latency, file.path(output, "deterministic", "latency_shift.csv"))
  prediction <- branching_prediction(lmodel)
  cshift <- tail(latency$deterministic_shift, 1)
  write_csv(data.frame(shift = cshift, survival = prediction["survival"],
    lambda = prediction["lambda"], conditional_mean = cshift +
      (0.5772156649015329 + log(prediction["survival"])) / prediction["lambda"],
    conditional_variance = pi^2 / (6 * prediction["lambda"]^2),
    conditional_skewness = 12 * sqrt(6) * 1.202056903159594 / pi^3),
    file.path(output, "deterministic", "latency_limit.csv"))
  stability <- list()
  for (N in c(50, 200)) for (speed in c(.25, .5, 1, 2, 4, 8)) for (family in c("clearance", "transit")) {
    dy0 <- if (N == 200) 1 else 4
    z <- if (family == "clearance") minimal_model(N, dx = speed, dy = speed * dy0, k = speed) else
      transit_model(N, speed, dy0)
    stability[[length(stability) + 1L]] <- data.frame(family = family, N = N, speed = speed,
      r = equilibrium(z)$r, minimum_return_rate = minimum_return_rate(z,
        intervals = if (profile == "smoke") 8L else 64L),
      routh_margin = if (family == "clearance") closed_routh_margin(z) else NA_real_)
  }
  write_csv(bind_rows(stability), file.path(output, "deterministic", "stability.csv"))
  invisible(list(therapy = bind_rows(therapy), latency = latency, stability = bind_rows(stability)))
}

run_audits <- function(output, profile = "smoke", seed = 20261028L) {
  set.seed(seed + 400000000L)
  n <- if (profile == "full") 100000L else 200L
  samples <- matrix(10^runif(4 * n, -2, 2), ncol = 4)
  margins <- numeric(n); band <- logical(n)
  a <- .75
  for (i in seq_len(n)) {
    p <- samples[i, 1]; q <- samples[i, 2]; dx <- samples[i, 3]; dy <- samples[i, 4]
    m <- minimal_model(100, p = p, q = q, dx = dx, dy = dy)
    margins[i] <- closed_routh_margin(m)
    band[i] <- dy / dx >= (2 * a - 1) / (2 * a) && dy / dx <= 2 * a / (2 * a - 1)
  }
  write_csv(data.frame(draw = seq_len(n), p = samples[, 1], q = samples[, 2],
    dx = samples[, 3], dy = samples[, 4], sufficient_band = band,
    routh_margin = margins, stable_on_closed_segment = margins > 0),
    file.path(output, "audits", "random_stability_draws.csv"))
  write_csv(data.frame(draws = n, stable_fraction = mean(margins > 0),
    unstable_percent = 100 * mean(margins <= 0),
    band_draws = sum(band), unstable_inside_sufficient_band = sum(band & margins <= 0)),
    file.path(output, "audits", "random_stability_summary.csv"))
  grid <- expand.grid(q = 10^seq(-1.5, 1.5, length.out = 25),
                      r = 10^seq(-1.5, 1.5, length.out = 25))
  grid$fixation_at_quarter <- grid$rare_clone_bias <- numeric(nrow(grid))
  for (i in seq_len(nrow(grid))) {
    m <- minimal_model(200, q = grid$q[i], dy = grid$q[i] * grid$r[i])
    sc <- scale_function(m)
    grid$fixation_at_quarter[i] <- sc$probability(sc$K / 2)
    grid$rare_clone_bias[i] <- sc$rare_clone_bias
  }
  write_csv(grid, file.path(output, "audits", "turnover_grid.csv"))
  unstable <- minimal_model(100, p = .05, q = 50, a = .9, dx = .5, dy = .02)
  eq <- equilibrium(unstable); psi_x <- .7 * (2 * unstable$a - 1)
  xi <- psi_x / (unstable$k * eq$A)
  vals <- eigen(model_jacobian(coexistence(xi, unstable), unstable), only.values = TRUE)$values
  write_csv(data.frame(real = Re(vals), imaginary = Im(vals), xi = xi, routh_margin =
    routh_margin(unstable, psi_x)), file.path(output, "audits", "unstable_example.csv"))
  invisible(TRUE)
}
