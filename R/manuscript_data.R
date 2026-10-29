reported_numbers <- function(x) {
  x <- as.character(x)
  x <- gsub("$", "", x, fixed = TRUE)
  if (grepl("\\\\times", x)) {
    parts <- regmatches(x, gregexpr("-?[0-9]+(?:[.][0-9]+)?", x, perl = TRUE))[[1]]
    return(as.numeric(parts[1]) * 10^as.numeric(tail(parts, 1)))
  }
  if (grepl("^[0-9]+/[0-9]+$", x)) {
    parts <- as.numeric(strsplit(x, "/", fixed = TRUE)[[1]])
    return(parts[1] / parts[2])
  }
  as.numeric(regmatches(x, gregexpr("-?[0-9]+(?:[.][0-9]+)?", x, perl = TRUE))[[1]])
}

read_figure_csv <- function(directory, name) {
  file <- file.path(directory, name)
  if (!file.exists(file)) return(data.frame())
  read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
}

reported_figure_tables <- function() {
  directory <- system.file("extdata", "reported", package = "stemfeedback")
  if (!nzchar(directory)) stop("The installed reported-table transcriptions are missing.")
  table <- function(name) read_figure_csv(directory, paste0(name, ".csv"))
  number <- function(z, column) vapply(z[[column]], function(x) reported_numbers(x)[1], 0.)
  interval <- function(z, column) t(vapply(z[[column]], reported_numbers, numeric(3)))
  chain <- table("chain"); ci <- interval(chain, 6)
  baseline <- data.frame(case_id = c(paste0("balanced_", 1:8),
    paste0("equal_clearance_", 1:4), paste0("rare_", c(.25, 1, 4))),
    suite = "baseline", q = number(chain, 1), r = number(chain, 2), N = number(chain, 3),
    attempts = number(chain, 5), estimate = ci[, 1], lower = ci[, 2], upper = ci[, 3],
    diffusion_prediction = number(chain, 7))
  weak <- table("sel"); ci <- interval(weak, 5)
  selection <- data.frame(case_id = paste0("selection_", 1:5), suite = "weak_selection",
    q = number(weak, 1), theta = number(weak, 2), N = number(weak, 3),
    attempts = number(weak, 4), estimate = ci[, 1], lower = ci[, 2], upper = ci[, 3],
    diffusion_prediction = number(weak, 6))
  founder <- table("founder"); ci <- interval(founder, 4)
  founders <- data.frame(case_id = paste("founder", number(founder, 1), number(founder, 2), sep = "_"),
    suite = "founders", r = number(founder, 1), N = number(founder, 2),
    attempts = 100000L, healthy_losses = number(founder, 3),
    scaled_estimate = ci[, 1], scaled_lower = ci[, 2], scaled_upper = ci[, 3],
    founder_coefficient = number(founder, 7))
  stage <- table("foundersplit"); ci <- interval(stage, 5)
  staged <- data.frame(r = number(stage, 1), N = number(stage, 2),
    scaled_estimate = ci[, 1], scaled_lower = ci[, 2], scaled_upper = ci[, 3])
  families <- lapply(c("relax", "depth"), function(label) {
    z <- table(label); family <- if (label == "relax") "clearance" else "transit"
    column <- if (label == "relax") 8 else 7; ci <- interval(z, column)
    data.frame(case_id = paste(family, number(z, 1), number(z, 2), sep = "_"),
      suite = family, N = number(z, 1), speed = number(z, 2),
      attempts = number(z, column - 1L), estimate = ci[, 1], lower = ci[, 2], upper = ci[, 3])
  })
  local <- table("localrates")
  rates <- bind_rows(lapply(c(200, 50), function(N) bind_rows(lapply(1:2, function(j)
    data.frame(N = N, family = c("clearance", "transit")[j], speed = number(local, 1),
      local_return_rate = number(local, j + if (N == 200) 1 else 3))))))
  noise <- table("controlled"); ci <- interval(noise, 5)
  controlled <- data.frame(case_id = paste("noise", number(noise, 1), number(noise, 2), sep = "_"),
    suite = "noise", N = number(noise, 1), noise = number(noise, 2),
    attempts = number(noise, 3), healthy_losses = number(noise, 4),
    estimate = ci[, 1], lower = ci[, 2], upper = ci[, 3])
  fluctuations <- data.frame(case_id = controlled$case_id, N = controlled$N,
    noise = controlled$noise, signal_rms = number(noise, 6), rms_lower = NA_real_, rms_upper = NA_real_)
  band <- table("occupationfrequency")
  boundaries <- lapply(strsplit(band[[2]], "--", fixed = TRUE), as.numeric)
  bands <- data.frame(noise = number(band, 1),
    band_lower = vapply(boundaries, `[`, 0., 1), band_upper = vapply(boundaries, `[`, 0., 2),
    integral = number(band, 3), mean_feedback_departure = number(band, 4),
    product_of_means = number(band, 5), covariance_contribution = number(band, 6),
    integral_se = NA_real_)
  times <- table("occupationtimes")
  occurrence <- bind_rows(lapply(seq_len(nrow(times)), function(i) data.frame(
    noise = number(times, 1)[i], time = c(20, 80),
    cumulative_healthy_loss = c(number(times, 4)[i], number(times, 5)[i]),
    cumulative_healthy_loss_lower = NA_real_, cumulative_healthy_loss_upper = NA_real_)))
  est <- table("est"); ci <- interval(est, 6)
  establishment <- data.frame(suite = "establishment", N = number(est, 1), threshold = number(est, 2),
    b = number(est, 3), q = number(est, 4), attempts = number(est, 5),
    estimate = ci[, 1], lower = ci[, 2], upper = ci[, 3])
  lat <- table("lat")
  mean_se <- t(vapply(lat[[3]], reported_numbers, numeric(2)))

  latency <- data.frame(case_id = paste0("latency_", number(lat, 1)), N = number(lat, 1),
    hits = vapply(strsplit(lat[[2]], "/", fixed = TRUE), function(x) as.numeric(x[1]), 0.),
    empirical_mean = number(lat, 4), mean_se = mean_se[, 2],
    mean_lower = NA_real_, mean_upper = NA_real_, complete_run = TRUE)
  macro <- table("main_macros")
  value <- function(name) reported_numbers(macro$value[match(name, macro$name)])[1]
  decomposition <- data.frame(noise = c(0, 16),
    integral = c(value("OccIntegralZero"), value("OccIntegralAdded")),
    product_of_means = c(value("OccProductZero"), value("OccProductAdded")),
    covariance_contribution = c(value("OccAssociationZero"), value("OccAssociationAdded")),
    product_se = NA_real_, covariance_se = NA_real_, integral_se = NA_real_)
  refs <- data.frame(case_id = baseline$case_id, N = baseline$N,
    probability = number(chain, 8))
  refs <- refs[is.finite(refs$probability), ]
  frefs <- number(founder, 5) / founders$N
  refs <- rbind(refs, data.frame(case_id = founders$case_id[is.finite(frefs)],
    N = founders$N[is.finite(frefs)], probability = frefs[is.finite(frefs)]))
  list(all = bind_rows(c(list(baseline, selection, founders, controlled, establishment), families)),
    staged = staged, fluctuations = fluctuations, bands = bands,
    occupation_time = occurrence, decomposition = decomposition, rates = rates,
    references = refs, latency = latency, latency_paths = data.frame(),
    provenance = read_figure_csv(directory, "provenance.csv"))
}

manuscript_trajectory_designs <- function() {
  ans <- list()
  for (family in c("two_compartment", "deeper_hierarchy")) {
    m <- if (family == "two_compartment") minimal_model(100, q = 3, dy = 1.5) else
      cell_model(100, healthy = lineage(c(1, 2, 3), c(.75, .55, .45), 1),
        leukemic = lineage(c(3, 1.5), c(.75, .6), .7))
    eq <- equilibrium(m)
    pairs <- rbind(c(.04, .9), c(.4, .4), c(1.4, .05))
    for (i in 1:3) {
      initial <- c(eq$c * pairs[i, 1] * eq$K, eq$e * pairs[i, 2] * eq$r * eq$K)
      initial[m$dimension] <- initial[m$dimension] * .65
      ans[[paste0(family, "_", i)]] <- list(model = m, initial = initial, family = family,
        provenance = "Illustrative new starts and terminal-only weights; historical settings missing")
    }
  }
  ans
}

validate_trajectory_designs <- function(designs) {
  if (!length(designs) || is.null(names(designs)) || any(!nzchar(names(designs))) ||
      anyDuplicated(names(designs))) stop("Provide uniquely named trajectory designs.")
  for (design in designs) {
    if (!all(c("model", "initial", "family", "provenance") %in% names(design)))
      stop("Each trajectory needs model, initial densities, family and provenance.")
    validate_state(design$initial, design$model, FALSE)
    stopifnot(all(design$initial[design$model$stem] > 0), length(design$family) == 1L,
      nzchar(design$family), length(design$provenance) == 1L, nzchar(design$provenance))
  }
  invisible(TRUE)
}

write_manuscript_deterministic <- function(output, designs = NULL, profile = "quick") {
  profile <- match.arg(profile, c("smoke", "quick", "full"))
  run_deterministic(output, profile, designs)
  invisible(output)
}

prepare_manuscript_theory <- function(output, references = FALSE, reference_cases = NULL, multiplier = 4) {
  stopifnot(length(output) == 1L, nzchar(output), multiplier >= 2)
  cases <- paper_designs()
  rates <- bind_rows(lapply(cases, function(case) {
    if (!case$suite %in% c("clearance", "transit")) return(NULL)
    m <- case$model; J <- project_stems(case$initial[1] / m$N,
      case$initial[m$stem[2]] / m$N, m)
    data.frame(case_id = case$id, family = case$suite, N = m$N,
      speed = if (case$suite == "clearance") m$healthy$clearance else m$healthy$division[2],
      projected_start = J, local_return_rate = return_rate(J, m))
  }))
  write_csv(rates, file.path(output, "theory", "figure_local_rates.csv"))
  refrows <- data.frame()
  if (references) {
    ids <- if (is.null(reference_cases)) reported_figure_tables()$references$case_id else reference_cases
    if (!all(ids %in% names(cases))) stop("Unknown reference case ID.")
    refrows <- bind_rows(lapply(ids, function(id) {
      case <- cases[[id]]; m <- case$model; eq <- equilibrium(m)
      limits <- pmax(2L, ceiling(multiplier * m$N * c(eq$K, eq$r * eq$K)),
        case$initial[m$stem] + 1L)
      message("Reference for ", id, ": ", paste(limits, collapse = " x "))
      ref <- solve_reference(m, limits)
      v <- reference_values(case$initial[1], case$initial[m$stem[2]], ref, m)
      data.frame(case_id = id, N = m$N, probability = v["U"],
        escape = v["escape"], solve_residual = ref$residual, outer_X1 = limits[1], outer_Y1 = limits[2])
    }))
    write_csv(refrows, file.path(output, "theory", "figure_references.csv"))
  }
  invisible(list(rates = rates, references = refrows))
}

manuscript_figure_data <- function(output = NULL, source = c("run", "reported")) {
  source <- match.arg(source)
  if (source == "reported") ans <- reported_figure_tables() else {
    if (is.null(output)) stop("Run mode requires an output directory.")
    summary <- function(name) read_figure_csv(file.path(output, "summary"), paste0(name, ".csv"))
    ans <- list(all = summary("all_cases"), staged = summary("staged_founder"),
      fluctuations = summary("feedback_fluctuations"), bands = summary("occupation_bands"),
      occupation_time = summary("occupation_time"), decomposition = summary("occupation_decomposition"),
      latency = summary("latency_moments"), latency_paths = summary("latency_hit_times"),
      rates = read_figure_csv(file.path(output, "theory"), "figure_local_rates.csv"),
      references = read_figure_csv(file.path(output, "theory"), "figure_references.csv"),
      provenance = read_figure_csv(output, "run_manifest.csv"))
    if (nrow(ans$provenance) != 1L || !"profile" %in% names(ans$provenance))
      stop("Run mode requires a recognized experiment manifest.")
    if (nrow(ans$all)) {
      ans$all$theta <- ans$all$N * (ans$all$b - ans$all$a)
      ans$all$speed <- ifelse(ans$all$suite == "clearance", ans$all$dx, NA_real_)
      cases <- paper_designs()
      for (i in which(ans$all$suite == "transit"))
        ans$all$speed[i] <- cases[[ans$all$case_id[i]]]$model$healthy$division[2]
      est <- ans$all$suite == "establishment"
      ans$all$threshold <- NA_real_
      for (i in which(est)) ans$all$threshold[i] <- cases[[ans$all$case_id[i]]]$threshold
    }
  }
  ans$source <- source
  ans$trajectories <- ans$segments <- ans$trajectory_settings <- ans$therapy <- ans$latency_limit <- data.frame()
  ans$occupation_paths <- list()
  if (!is.null(output)) {
    directory <- file.path(output, "deterministic")
    ans$trajectories <- read_figure_csv(directory, "coexistence_trajectories.csv")
    ans$segments <- read_figure_csv(directory, "coexistence_segments.csv")
    ans$trajectory_settings <- read_figure_csv(directory, "trajectory_settings.csv")
    ans$therapy <- read_figure_csv(directory, "therapy_trajectories.csv")
    ans$latency_limit <- read_figure_csv(directory, "latency_limit.csv")
    if (source == "run") for (nu in c(0, 16))
      ans$occupation_paths[[as.character(nu)]] <- read_figure_csv(file.path(output, "paths"),
        paste0("occupation_", nu, ".csv"))
  }
  ans
}

run_manuscript <- function(output, profile = "smoke", seed = 20261028L, workers = 1L,
                           replicates = NULL, max_events = 2e7, references = FALSE, force = FALSE) {
  run_experiments(output, suite = "all", profile = profile, seed = seed, workers = workers,
    replicates = replicates, max_events = max_events, force = force, plots = FALSE)
  prepare_manuscript_theory(output, references = references)
  plot_manuscript_figures(output, source = "run")
  invisible(output)
}
