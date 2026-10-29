simulate_path <- function(initial, model, seed, stop_rule = "first_stem_loss",
                           threshold = NULL, burden = NULL, max_events = 2e7,
                           time_limit = Inf, observation_times = seq(0, 2, .25),
                           reference = NULL, enlarged_reference = NULL,
                           bands = c(0, .1, .25, .5, .75, .9, 1)) {
  validate_state(initial, model)
  if (sum(initial[model$stem]) == 0)
    stop("At least one stem pool must be present in the initial state.")
  stopifnot(stop_rule %in% c("first_stem_loss", "founder_entry", "establishment", "burden"),
            length(seed) == 1, is.finite(seed), seed > 0, seed < .Machine$integer.max,
            max_events >= 1, time_limit > 0, all(diff(bands) > 0),
            bands[1] == 0, tail(bands, 1) == 1)
  if (stop_rule %in% c("founder_entry", "establishment"))
    stopifnot(length(threshold) == 1, threshold >= 1, threshold == round(threshold))
  if (stop_rule == "burden") stopifnot(length(burden) == 1, burden >= 1)
  if (!is.null(enlarged_reference) && is.null(reference)) stop("A primary reference is required.")
  if (!is.null(reference) && stop_rule != "first_stem_loss")
    stop("The feedback integral stops at first stem loss or reference exit.")
  if (!is.null(reference)) validate_reference(reference, model)
  if (!is.null(enlarged_reference)) {
    validate_reference(enlarged_reference, model)
    if (any(enlarged_reference$limits < reference$limits))
      stop("The enlarged reference must contain the primary reference domain.")
  }
  if (!is.null(reference) && (initial[model$stem[1]] >= reference$limits[1] ||
                              initial[model$stem[2]] >= reference$limits[2]))
    stop("Initial stem counts must be strictly inside the reference rectangle.")
  RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(as.integer(seed))
  obs_times <- sort(unique(c(0, observation_times)))
  stopifnot(all(is.finite(obs_times)), all(obs_times >= 0))
  obs <- matrix(NA_real_, length(obs_times), 5 + model$dimension)
  colnames(obs) <- c("time", "signal", "integral", "reference", "healthy_loss", model$names)
  obs[, "time"] <- obs_times
  state <- as.numeric(initial); t <- 0; events <- 0L; integral <- 0; integral2 <- 0
  error_bound <- 0; positive <- 0; negative <- 0; squared_signal <- 0
  nb <- length(bands) - 1L; occupation <- matrix(0, nb, 4)
  colnames(occupation) <- c("time", "feedback", "response", "integral")
  oi <- 1L; status <- "running"; rms_survivor <- FALSE; transverse <- NA_real_
  cached_stems <- c(-1, -1); r <- r2 <- NULL
  repeat {
    x <- state[model$stem[1]]; y <- state[model$stem[2]]
    if (stop_rule == "burden" && state[model$dimension] >= burden) status <- "burden_hit"
    else if (stop_rule %in% c("founder_entry", "establishment") && y >= threshold) status <- "threshold_hit"
    else if (y == 0) status <- "leukemic_loss"
    else if (x == 0 && stop_rule %in% c("first_stem_loss", "founder_entry")) status <- "healthy_loss"
    else if (!is.null(reference) && (x >= reference$limits[1] || y >= reference$limits[2])) status <- "reference_exit"
    if (status != "running") break
    if (events >= max_events) { status <- "event_limit"; break }
    if (t >= time_limit) { status <- "time_limit"; break }
    s <- feedback(state, model)
    if (!is.null(reference) && any(c(x, y) != cached_stems)) {
      r <- reference_values_unchecked(x, y, reference)
      if (!is.null(enlarged_reference)) r2 <- reference_values_unchecked(x, y, enlarged_reference)
      cached_stems <- c(x, y)
    }
    rates <- event_rates(state, model); total <- sum(rates)
    if (!is.finite(total) || total <= 0 || any(rates < 0)) stop("Invalid event rates.")
    waiting <- rexp(1, total); next_time <- min(t + waiting, time_limit); dt <- next_time - t
    if (dt <= 0 || !is.finite(dt)) stop("Invalid or unresolvable event time.")
    g <- if (is.null(reference)) 0 else (s - r["signal"]) * r["Psi"]
    g2 <- if (is.null(enlarged_reference)) 0 else (s - r2["signal"]) * r2["Psi"]
    while (oi <= nrow(obs) && obs_times[oi] <= next_time) {
      obs[oi, -1] <- c(s, if (is.null(reference)) NA else integral + g * (obs_times[oi] - t),
                       if (is.null(reference)) NA else r["U"], 0, state)
      if (obs_times[oi] == 2) {
        rms_survivor <- x > 0 && y > 0
        J <- project_stems(x / model$N, y / model$N, model); eq <- equilibrium(model)
        expected <- eq$r * (eq$K - J)
        if (expected > 0) transverse <- log(y / model$N / expected)
      }
      oi <- oi + 1L
    }
    squared_signal <- squared_signal + max(0, min(next_time, 2) - min(t, 2)) *
      (s - 1 / (2 * model$a))^2
    if (!is.null(reference)) {
      integral <- integral + g * dt; integral2 <- integral2 + g2 * dt
      error_bound <- error_bound + abs(s - r["signal"]) * r["error_response"] * dt
      if (s >= r["signal"]) positive <- positive + g * dt else negative <- negative + g * dt
      band <- min(nb, findInterval(y / (x + y), bands, rightmost.closed = TRUE))
      occupation[band, ] <- occupation[band, ] +
        c(dt, (s - r["signal"]) * dt, r["Psi"] * dt, g * dt)
    }
    t <- next_time
    if (t >= time_limit) { status <- "time_limit"; break }
    chosen <- which(cumsum(rates) > runif(1) * total)[1]
    if (is.na(chosen)) stop("Could not select a cell event.")
    state <- state + model$jumps[chosen, ]; events <- events + 1L
    if (any(state < 0)) stop("A cell event produced a negative count.")
  }
  U0 <- Ufinal <- NA_real_
  if (!is.null(reference)) {
    U0 <- reference$U[initial[1] + 1L, initial[model$stem[2]] + 1L]
    Ufinal <- reference$U[state[1] + 1L, state[model$stem[2]] + 1L]
  }

  if (!status %in% c("time_limit", "event_limit")) {
    while (oi <= nrow(obs)) {
      obs[oi, -1] <- c(feedback(state, model), if (is.null(reference)) NA else integral,
                       Ufinal, as.numeric(status == "healthy_loss"), state)
      oi <- oi + 1L
    }
  }
  list(seed = seed, status = status, time = t, events = events,
       initial = setNames(initial, model$names), final = setNames(state, model$names),
       healthy_loss = as.integer(status == "healthy_loss"),
       hit = as.integer(status %in% c("threshold_hit", "burden_hit")),
       censored = status %in% c("event_limit", "time_limit"),
       rms = if (rms_survivor) sqrt(squared_signal / 2) else NA_real_,
       transverse = transverse, observations = obs,
       integral = if (is.null(reference)) NA_real_ else unname(integral),
       enlarged_integral = if (is.null(enlarged_reference)) NA_real_ else unname(integral2),
       integration_error_bound = if (is.null(reference)) NA_real_ else unname(error_bound),
       positive_departure_integral = if (is.null(reference)) NA_real_ else unname(positive),
       negative_departure_integral = if (is.null(reference)) NA_real_ else unname(negative),
       reference_initial = U0, reference_final = Ufinal,
       residual = if (is.null(reference)) NA_real_ else unname(Ufinal - U0 - integral),
       occupation = occupation, bands = bands)
}

path_record <- function(path, attempt, case_id) {
  row <- data.frame(case_id = case_id, attempt = attempt, seed = path$seed,
    status = path$status, time = path$time, events = path$events,
    healthy_loss = path$healthy_loss, hit = path$hit, censored = path$censored,
    rms = path$rms, transverse = path$transverse, integral = path$integral,
    enlarged_integral = path$enlarged_integral, integration_error_bound = path$integration_error_bound,
    positive_departure_integral = path$positive_departure_integral,
    negative_departure_integral = path$negative_departure_integral,
    reference_initial = path$reference_initial, reference_final = path$reference_final,
    residual = path$residual)
  for (name in names(path$final)) row[[paste0("final_", name)]] <- path$final[[name]]
  if (is.finite(path$reference_initial))
    for (b in seq_len(nrow(path$occupation))) for (j in colnames(path$occupation))
      row[[paste0("band", b, "_", j)]] <- path$occupation[b, j]
  if (nrow(path$observations) > 1) for (i in seq_len(nrow(path$observations))) {
    tag <- gsub("[.]", "p", format(path$observations[i, "time"], trim = TRUE, scientific = FALSE))
    for (j in c("signal", "integral", "reference", "healthy_loss"))
      row[[paste0("t", tag, "_", j)]] <- path$observations[i, j]
  }
  row
}
