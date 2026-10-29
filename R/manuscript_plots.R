figure_rows <- function(z, field, values) {
  if (!nrow(z) || !field %in% names(z)) return(z[FALSE, , drop = FALSE])
  z[z[[field]] %in% values, , drop = FALSE]
}

figure_range <- function(..., fallback = c(0, 1)) {
  values <- unlist(list(...)); values <- values[is.finite(values)]
  if (!length(values)) return(fallback)
  limits <- range(values); padding <- max(diff(limits) * .08, max(abs(limits)) * .01, 1e-6)
  limits + c(-padding, padding)
}

figure_log_range <- function(...) {
  values <- unlist(list(...)); values <- values[is.finite(values) & values > 0]
  if (!length(values)) return(c(.1, 10))
  limits <- range(log(values)); padding <- max(diff(limits) * .05, .05)
  exp(limits + c(-padding, padding))
}

figure_intervals <- function(x, estimate, lower, upper, col = "#176B87", pch = 16) {
  graphics::points(x, estimate, col = col, pch = pch, cex = .8)
  ok <- is.finite(x) & is.finite(lower) & is.finite(upper) & upper > lower
  if (any(ok)) graphics::arrows(x[ok], lower[ok], x[ok], upper[ok],
    code = 3, angle = 90, length = .035, col = col)
}

figure_note <- function(text) graphics::mtext(text, side = 3, line = .25, adj = 0,
                                            cex = .62, col = "gray35")

figure_decomposition_range <- function(z) {
  figure_range(0, z$product_of_means, z$covariance_contribution, z$integral,
    z$product_of_means - 1.96 * z$product_se, z$product_of_means + 1.96 * z$product_se,
    z$covariance_contribution - 1.96 * z$covariance_se,
    z$covariance_contribution + 1.96 * z$covariance_se,
    z$integral - 1.96 * z$integral_se, z$integral + 1.96 * z$integral_se)
}

blank_figure_panel <- function(text) {
  graphics::plot.new(); graphics::plot.window(c(0, 1), c(0, 1))
  graphics::text(.5, .55, "Inputs unavailable", font = 2, cex = .95)
  graphics::text(.5, .4, paste(strwrap(text, 46), collapse = "\n"), cex = .75, col = "gray35")
}

figure_reference <- function(data, id) {
  z <- figure_rows(data$references, "case_id", id)
  if (nrow(z)) z$probability[1] else NA_real_
}

figure_resolved <- function(z, count) {
  nrow(z) == count && all(is.finite(z$estimate) & is.finite(z$lower) & is.finite(z$upper))
}

figure_limit <- function(data) {
  if (nrow(data$latency_limit)) return(data$latency_limit[1, ])
  data.frame(shift = -6.79531, survival = .5, lambda = 1/3,
    conditional_mean = -6.79531 + 3 * (.5772156649015329 + log(.5)))
}

latency_distribution_data <- function(values, shift, survival, lambda) {
  stopifnot(length(values) >= 1L, all(is.finite(values)), survival > 0, survival <= 1, lambda > 0)
  values <- sort(values); n <- length(values)
  F <- function(x) exp(-survival * exp(-lambda * (x - shift)))
  unique_values <- unique(values)
  counts <- tabulate(match(values, unique_values))
  right <- cumsum(counts) / n; left <- (cumsum(counts) - counts) / n
  ks <- max(abs(right - F(unique_values)), abs(left - F(unique_values)))
  grid <- sort(unique(c(values, seq(min(values) - 2, max(values) + 2, length.out = 401))))
  probability <- stats::ppoints(n)
  list(cdf = data.frame(time = grid, difference = stats::ecdf(values)(grid) - F(grid)),
    qq = data.frame(predicted = shift + (log(survival) - log(-log(probability))) / lambda,
                    empirical = values), ks = ks, hits = n)
}

plot_manuscript_figures <- function(output, source = c("run", "reported"), separate = FALSE,
                                    figures = c(paste0("Fig", 1:6), "FigS1"), data = NULL) {
  source <- match.arg(source)
  stopifnot(length(output) == 1L, nzchar(output), length(separate) == 1L,
    is.logical(separate), !is.na(separate), length(figures) > 0,
    !anyDuplicated(figures), all(figures %in% c(paste0("Fig", 1:6), "FigS1")))
  if (is.null(data)) data <- manuscript_figure_data(output, source)
  if (!identical(data$source, source)) stop("The figure-data source and requested source differ.")
  directory <- file.path(output, "manuscript_figures")
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  if (!nrow(data$rates) && source == "run")
    data$rates <- prepare_manuscript_theory(output)$rates
  coverage <- list(); files <- character(); owned_device <- NULL
  on.exit(if (!is.null(owned_device) && owned_device %in% grDevices::dev.list())
    grDevices::dev.off(owned_device), add = TRUE)
  open <- function(file) {
    grDevices::pdf(file, width = 9.5, height = 7.1, family = "Helvetica", useDingbats = FALSE)
    owned_device <<- grDevices::dev.cur(); files <<- c(files, file)
  }
  if (!separate) open(file.path(directory, paste0("manuscript_figures_", source, ".pdf")))
  record <- function(fig, panel, status = "available", detail = "") {
    graphics::mtext(panel, side = 3, line = 1.1, adj = 0, font = 2, cex = 1.1)
    coverage[[length(coverage) + 1L]] <<- data.frame(figure = fig, panel = panel,
      status = status, source = source, detail = detail)
  }
  missing <- function(fig, panel, detail) {
    blank_figure_panel(detail); record(fig, panel, "missing", detail)
  }
  probability_panel <- function(z, x, label, logx = FALSE, theory = NULL, prediction = NULL,
                                ylim = NULL, col = "#176B87", pch = 16, ticks = NULL) {
    if (is.null(ylim)) ylim <- figure_range(z$lower, z$upper, z$estimate, prediction)
    limits <- if (logx) figure_log_range(x, if (!is.null(theory)) theory$x) else
      figure_range(x, if (!is.null(theory)) theory$x)
    graphics::plot(NA, xlim = limits, ylim = ylim,
      xlab = label, ylab = "Leukemic fixation probability", log = if (logx) "x" else "",
      xaxt = if (is.null(ticks)) "s" else "n")
    if (!is.null(ticks)) graphics::axis(1, ticks, ticks)
    if (!is.null(theory)) graphics::lines(theory$x, theory$y, lwd = 1.4, col = "gray20")
    if (nrow(z)) figure_intervals(x, z$estimate, z$lower, z$upper, col, pch)
    if (!is.null(prediction)) graphics::points(x, prediction, pch = 4, col = "#C4522B", cex = .9)
  }
  for (fig in figures) {
    if (separate) open(file.path(directory, paste0(fig, "_", source, ".pdf")))
    graphics::par(mfrow = if (fig %in% c("Fig1", "Fig3", "Fig5")) c(1, 2) else c(2, 2),
      mar = c(4.1, 4.4, 2.5, .8), oma = c(2.4, 0, 2.4, 0), mgp = c(2.7, .8, 0),
      las = 1, bty = "l", cex.axis = .8, cex.lab = .9)
    if (fig == "Fig1") {
      z <- data$trajectories; seg <- data$segments
      if (!nrow(z)) missing(fig, "a", "Full trajectory starts, feedback weights and ODE trajectories are required.") else {
        graphics::plot(NA, xlim = range(z$X1[z$X1 > 0]), ylim = range(z$Y1[z$Y1 > 0]),
          log = "xy", xlab = "Healthy stem density", ylab = "Leukemic stem density")
        families <- unique(z$family); cols <- c("#176B87", "#C4522B")
        for (j in seq_along(families)) {
          family <- families[j]
          s <- figure_rows(seg, "family", family)
          if (nrow(s)) {
            s <- s[!duplicated(s$X1), ]; s <- s[s$X1 > 0 & s$Y1 > 0, ]
            graphics::lines(s$X1, s$Y1, col = grDevices::adjustcolor(cols[j], .3), lwd = 3)
          }
          for (id in unique(z$trajectory_id[z$family == family])) {
            v <- figure_rows(z, "trajectory_id", id)
            graphics::lines(v$X1, v$Y1, col = cols[j], lty = j, lwd = 1.1)
            graphics::points(v$X1[c(1, nrow(v))], v$Y1[c(1, nrow(v))], pch = c(1, 16), col = cols[j])
          }
        }
        graphics::legend("bottomleft", c("2+2 compartments", "4+3 compartments"),
          col = cols, lty = 1:2, bty = "n", cex = .7)
        illustrative <- any(grepl("[Ii]llustrative", z$provenance))
        if (illustrative) figure_note("Illustrative new starts and weights")
        record(fig, "a", if (illustrative) "illustrative" else "available",
          paste(unique(z$provenance), collapse = "; "))
      }
      z <- data$therapy
      if (!nrow(z)) missing(fig, "b", "Run write_manuscript_deterministic() to calculate the reported depletion schedules.") else {
        graphics::plot(NA, xlim = c(0, 30), ylim = c(0, 1), xlab = "Time (1/p)",
          ylab = "Leukemic stem frequency")
        graphics::rect(0, 0, 3, 1, col = "gray94", border = NA)
        ids <- c("equal_independent", "common_division", "selective_1", "selective_2")
        labels <- c("Equal hazards", "Per-division kill", "Selective ratio 3", "Selective ratio 5")
        cols <- c("#176B87", "#C4522B", "#6651A1", "#4B783E")
        for (j in seq_along(ids)) {
          v <- figure_rows(z, "schedule", ids[j])
          graphics::lines(v$time, v$Y1 / (v$X1 + v$Y1), col = cols[j], lty = j, lwd = 1.3)
        }
        graphics::abline(h = 1/3, col = "gray60", lty = 3)
        graphics::legend("right", labels, col = cols, lty = 1:4, bty = "n", cex = .65)
        record(fig, "b", detail = "Independent deterministic calculation with reported schedules")
      }
    }
    if (fig == "Fig2") {
      all <- data$all; base <- figure_rows(all, "suite", "baseline")
      z <- figure_rows(base, "case_id", c("balanced_1", "balanced_2", "balanced_5"))
      q <- exp(seq(log(.2), log(5), length.out = 81))
      curve <- vapply(q, function(q) {
        m <- minimal_model(200, q = q); diffusion_probability(start_from_stems(50, 50, m), m)
      }, 0.)
      probability_panel(z, z$q, "Division-rate ratio q/p", TRUE, list(x = q, y = curve),
        ylim = c(0, 1), ticks = c(.25, .5, 1, 2, 4))
      graphics::abline(h = .5, col = "gray65", lty = 3)
      record(fig, "a", if (figure_resolved(z, 3)) "available" else "partial", "Three balanced-start turnover designs")
      z <- base[grepl("^balanced_", base$case_id) & base$q == 4, , drop = FALSE]
      N <- c(50, 100, 200, 400, 800, 1600)
      prediction <- vapply(z$case_id, function(id) {
        case <- paper_designs()[[id]]; diffusion_probability(case$initial, case$model)
      }, 0.)
      probability_panel(z, z$N, "System size N", TRUE,
        list(x = N, y = rep(.35, 6)), prediction, ticks = N)
      record(fig, "b", if (figure_resolved(z, 6)) "available" else "partial", "Equal amplification; crosses use rounded starts")
      z <- base[grepl("^equal_clearance_", base$case_id) & base$q == 4, , drop = FALSE]
      prediction <- vapply(z$case_id, function(id) {
        case <- paper_designs()[[id]]; diffusion_probability(case$initial, case$model)
      }, 0.)
      ideal <- scale_function(minimal_model(200, q = 4, dy = 1))$probability(.25)
      probability_panel(z, z$N, "System size N", TRUE,
        list(x = c(200, 3200), y = rep(ideal, 2)), prediction, ticks = c(200, 800, 3200))
      record(fig, "c", if (figure_resolved(z, 3)) "available" else "partial", "Equal terminal clearance")
      z <- figure_rows(all, "suite", "weak_selection"); z <- figure_rows(z, "q", 1)
      zero <- figure_rows(base, "case_id", "balanced_2"); zero$theta <- rep(0, nrow(zero))
      z <- bind_rows(list(z, zero))
      theta <- seq(-2.5, 2.5, length.out = 81)
      curve <- vapply(theta, function(t) {
        m <- minimal_model(200, q = 1, b = .75 + t / 200)
        diffusion_probability(start_from_stems(50, 50, m), m)
      }, 0.)
      probability_panel(z, z$theta, "Weak self-renewal difference N(b-a)",
        theory = list(x = theta, y = curve), ylim = c(0, 1))
      record(fig, "d", if (figure_resolved(z, 5)) "available" else "partial", "Reuses the equal-self-renewal observation from panel a")
    }
    if (fig == "Fig3") for (j in 1:2) {
      r <- c(1, .25)[j]; z <- figure_rows(figure_rows(data$all, "suite", "founders"), "r", r)
      staged <- figure_rows(data$staged, "r", r)
      N <- exp(seq(log(100), log(1600), length.out = 81))
      prediction <- vapply(N, function(n) {
        m <- minimal_model(n, dy = 4 * r)
        n * diffusion_probability(c(n/2, n/2, 1, 0), m)
      }, 0.)
      limit <- scale_function(minimal_model(100, dy = 4 * r))$founder_coefficient
      graphics::plot(NA, xlim = c(90, 1800), log = "x", xaxt = "n",
        ylim = figure_range(0, z$scaled_upper, staged$scaled_upper, prediction),
        xlab = "System size N", ylab = "N x fixation probability")
      graphics::axis(1, c(100, 200, 400, 800, 1600), c(100, 200, 400, 800, 1600))
      graphics::lines(N, prediction, col = "gray20", lwd = 1.3)
      graphics::abline(h = limit, col = "gray50", lty = 2)
      if (nrow(z)) figure_intervals(z$N, z$scaled_estimate, z$scaled_lower, z$scaled_upper)
      if (nrow(staged)) figure_intervals(staged$N, staged$scaled_estimate, staged$scaled_lower,
        staged$scaled_upper, pch = 18)
      ref <- figure_rows(data$references, "case_id", z$case_id)
      if (nrow(ref)) graphics::points(ref$N, ref$N * ref$probability, pch = 0, col = "#C4522B")
      figure_note(paste("r =", r))
      graphics::legend("topleft", c("Direct", "Staged", "Reference", "Diffusion", "Founder limit"),
        pch = c(16, 18, 0, NA, NA), lty = c(NA, NA, NA, 1, 2),
        col = c("#176B87", "#176B87", "#C4522B", "gray20", "gray50"), bty = "n", cex = .6)
      record(fig, letters[j], if (nrow(z) == 3 && nrow(staged) == 2 && nrow(ref) == 2 &&
        all(is.finite(z$scaled_estimate)) && all(is.finite(staged$scaled_estimate)))
        "available" else "partial", "Bootstrap resample count for historical staged intervals remains unreported")
    }
    if (fig == "Fig4") for (j in 1:4) {
      N <- c(200, 50, 200, 50)[j]; local <- j > 2
      clear <- figure_rows(figure_rows(data$all, "suite", "clearance"), "N", N)
      transit <- figure_rows(figure_rows(data$all, "suite", "transit"), "N", N)
      if (!nrow(clear) && !nrow(transit)) {
        missing(fig, letters[j], paste("No clearance/transit outcomes for N =", N)); next
      }
      get_x <- function(z, family) if (!local) z$speed else {
        rates <- figure_rows(figure_rows(data$rates, "family", family), "N", N)
        rates$local_return_rate[match(z$speed, rates$speed)]
      }
      x <- get_x(clear, "clearance"); xt <- get_x(transit, "transit")
      id <- if (N == 200) "equal_clearance_2" else "balanced_3"
      case <- paper_designs()[[id]]; diff <- diffusion_probability(case$initial, case$model)
      ref <- figure_reference(data, id)
      graphics::plot(NA, xlim = range(c(x, if (local) xt else numeric()), na.rm = TRUE), log = "x",
        xaxt = if (local) "s" else "n",
        ylim = figure_range(clear$lower, clear$upper, if (local) c(transit$lower, transit$upper), diff, ref),
        xlab = if (local) "Local return rate at projected start" else "Clearance factor",
        ylab = "Leukemic fixation probability")
      if (!local) graphics::axis(1, c(.25, .5, 1, 2, 4, 8), c("1/4", "1/2", "1", "2", "4", "8"))
      graphics::abline(h = diff, lwd = 1.3, col = "gray20")
      if (is.finite(ref)) graphics::abline(h = ref, col = "#C4522B", lty = 2)
      if (nrow(clear)) figure_intervals(x, clear$estimate, clear$lower, clear$upper)
      if (local && nrow(transit)) figure_intervals(xt, transit$estimate, transit$lower, transit$upper,
        col = "#6651A1", pch = 17)
      figure_note(if (local) paste("N =", N, "| Clearance: circles | Transit: triangles") else paste("N =", N))
      record(fig, letters[j], if (figure_resolved(clear, 6) && (!local || figure_resolved(transit, 6)) &&
        is.finite(ref)) "available" else "partial", "Horizontal rates are local, not full-segment minima")
    }
    if (fig == "Fig5") {
      z <- figure_rows(data$fluctuations, "N", 100)
      if (!nrow(z)) missing(fig, "a", "No N=100 conditional RMS summaries are available.") else {
        z <- z[order(z$noise), ]; x <- log1p(z$noise)
        graphics::plot(x, z$signal_rms, type = "l", axes = FALSE,
          ylim = figure_range(z$rms_lower, z$rms_upper, z$signal_rms),
          xlab = "Terminal noise rate", ylab = "Conditional mean signal RMS", col = "#176B87")
        graphics::axis(1, at = x, labels = z$noise); graphics::axis(2); graphics::box(bty = "l")
        figure_intervals(x, z$signal_rms, z$rms_lower, z$rms_upper)
        complete <- all(is.finite(z$rms_lower) & is.finite(z$rms_upper))
        if (!complete) figure_note("RMS interval endpoints unavailable")
        record(fig, "a", if (complete && nrow(z) == 4) "available" else "partial",
          "Conditioned on both stem pools surviving to time 2")
      }
      z <- figure_rows(figure_rows(data$all, "suite", "noise"), "N", 100)
      if (!nrow(z)) missing(fig, "b", "No N=100 terminal-noise outcomes are available.") else {
        z <- z[order(z$noise), ]; x <- log1p(z$noise)
        graphics::plot(x, z$estimate, type = "l", axes = FALSE, col = "#176B87",
          ylim = figure_range(z$lower, z$upper), xlab = "Terminal noise rate", ylab = "Fixation probability")
        graphics::axis(1, at = x, labels = z$noise); graphics::axis(2); graphics::box(bty = "l")
        figure_intervals(x, z$estimate, z$lower, z$upper)
        record(fig, "b", if (nrow(z) == 4 && all(is.finite(z$estimate))) "available" else "partial",
          "All attempts; Wilson intervals")
      }
    }
    if (fig == "Fig6") {
      cols <- c("#176B87", "#C4522B"); z <- data$occupation_time
      if (!nrow(z) || !"integral" %in% names(z))
        missing(fig, "a", "Time-resolved path integrals and standard errors are not printed in the tables.") else {
        graphics::plot(NA, xlim = range(log1p(z$time)),
          ylim = figure_range(0, z$integral, z$integral - 1.96*z$integral_se, z$integral + 1.96*z$integral_se),
          xaxt = "n", xlab = "Time (1/p)", ylab = "Accumulated feedback correction")
        ticks <- c(0, 2, 10, 40, 160); graphics::axis(1, log1p(ticks), ticks)
        graphics::abline(h = 0, col = "gray70")
        for (j in 1:2) {
          v <- figure_rows(z, "noise", c(0, 16)[j]); v <- v[order(v$time), ]
          x <- log1p(v$time); low <- v$integral - 1.96*v$integral_se; high <- v$integral + 1.96*v$integral_se
          ok <- is.finite(low) & is.finite(high)
          if (any(ok)) graphics::polygon(c(x[ok], rev(x[ok])), c(low[ok], rev(high[ok])),
            col = grDevices::adjustcolor(cols[j], .15), border = NA)
          graphics::lines(x, v$integral, col = cols[j], lty = j, lwd = 1.3)
        }
        record(fig, "a", if (length(unique(z$noise)) == 2 && all(is.finite(z$integral_se)))
          "available" else "partial", "Exact holding-time integrals; pointwise normal intervals")
      }
      z <- data$bands
      if (!nrow(z)) missing(fig, "b", "Pathwise frequency-band contributions are required.") else {
        graphics::plot(NA, xlim = c(.5, 6.5), ylim = figure_range(z$integral,
          z$integral - 1.96*z$integral_se, z$integral + 1.96*z$integral_se), xaxt = "n",
          xlab = "Leukemic stem frequency (%)", ylab = "Band contribution")
        graphics::abline(h = 0, col = "gray70")
        graphics::axis(1, 1:6, c("0-10", "10-25", "25-50", "50-75", "75-90", "90-100"), cex.axis = .7)
        for (j in 1:2) {
          v <- figure_rows(z, "noise", c(0, 16)[j]); v <- v[order(v$band_lower), ]
          figure_intervals(seq_len(nrow(v)) + c(-.1, .1)[j], v$integral,
            v$integral - 1.96*v$integral_se, v$integral + 1.96*v$integral_se, cols[j], c(16, 17)[j])
        }
        complete <- nrow(z) == 12 && all(is.finite(z$integral_se))
        if (!complete) figure_note("Band standard errors unavailable")
        record(fig, "b", if (complete) "available" else "partial", "Six frequency bands per noise rate")
      }
      z <- data$decomposition
      if (!nrow(z)) missing(fig, "c", "Joint path records are required for the decomposition and total errors.") else {
        z <- z[order(z$noise), ]
        graphics::plot(NA, xlim = c(.4, 2.6), ylim = figure_decomposition_range(z), xaxt = "n",
          xlab = "Terminal noise rate", ylab = "Fixation correction")
        graphics::axis(1, seq_len(nrow(z)), z$noise); graphics::abline(h = 0, col = "gray65")
        for (i in seq_len(nrow(z))) {
          graphics::rect(i-.3, 0, i, z$product_of_means[i], col = "gray75", border = NA)
          graphics::rect(i, 0, i+.3, z$covariance_contribution[i], col = "gray40", border = NA)
          figure_intervals(c(i-.15, i+.15),
            c(z$product_of_means[i], z$covariance_contribution[i]),
            c(z$product_of_means[i]-1.96*z$product_se[i], z$covariance_contribution[i]-1.96*z$covariance_se[i]),
            c(z$product_of_means[i]+1.96*z$product_se[i], z$covariance_contribution[i]+1.96*z$covariance_se[i]),
            col = "gray20", pch = NA)
          figure_intervals(i, z$integral[i], z$integral[i]-1.96*z$integral_se[i],
            z$integral[i]+1.96*z$integral_se[i], col = "gray20", pch = 18)
        }
        graphics::legend("bottomleft", c("Product of means", "Covariance", "Total"),
          fill = c("gray75", "gray40", NA), pch = c(NA, NA, 18), bty = "n", cex = .65)
        complete <- nrow(z) == 2 && all(is.finite(z$product_se) & is.finite(z$covariance_se) & is.finite(z$integral_se))
        if (!complete) figure_note("Joint uncertainty unavailable")
        record(fig, "c", if (complete) "available" else "partial", "Joint delta-method errors retain covariance across bands")
      }
      paths <- data$occupation_paths
      raw <- length(paths) == 2L && all(vapply(paths, nrow, 0L) > 0)
      z <- data$occupation_time
      if (!raw && !nrow(z)) missing(fig, "d", "Original stopping times and outcomes, or tabulated probabilities, are required.") else {
        maximum <- if (raw) max(vapply(paths, function(v) mean(v$healthy_loss), 0.)) else
          max(z$cumulative_healthy_loss, na.rm = TRUE)
        graphics::plot(NA, xlim = log1p(c(0, 160)), ylim = c(0, min(1, max(.4, 1.05*maximum))), xaxt = "n",
          xlab = "Time (1/p)", ylab = "Healthy stem-loss probability")
        ticks <- c(0, 2, 10, 40, 160); graphics::axis(1, log1p(ticks), ticks)
        for (j in 1:2) {
          if (raw) {
            v <- paths[[as.character(c(0, 16)[j])]]
            if (any(v$censored)) next
            ts <- sort(unique(c(0, v$time[v$healthy_loss == 1 & v$time <= 160], 160)))
            prob <- vapply(ts, function(t) mean(v$healthy_loss == 1 & v$time <= t), 0.)
            graphics::lines(log1p(ts), prob, type = "s", col = cols[j], lty = j, lwd = 1.2)
          } else {
            v <- figure_rows(z, "noise", c(0, 16)[j])
            figure_intervals(log1p(v$time), v$cumulative_healthy_loss,
              v$cumulative_healthy_loss_lower, v$cumulative_healthy_loss_upper, cols[j], c(16, 17)[j])
          }
        }
        if (!raw) figure_note("Only tabulated observation times")
        complete <- raw && !any(vapply(paths, function(x)
          any(x$censored) || any(x$status == "reference_exit"), FALSE))
        record(fig, "d", if (complete) "available" else "partial", "Empirical cumulative healthy loss; reference exits remain unresolved")
      }
      graphics::legend("topright", c("Noise 0", "Noise 16"), col = cols, lty = 1:2, bty = "n", cex = .65)
    }
    if (fig == "FigS1") {
      z <- figure_rows(figure_rows(data$all, "suite", "establishment"), "N", 2000)
      z <- figure_rows(z, "threshold", 40)
      graphics::plot(NA, xlim = c(.77, 1.02), ylim = c(0, min(1, max(.55, z$upper, na.rm = TRUE))),
        xlab = "Leukemic self-renewal probability",
        ylab = "Establishment probability")
      b <- seq(.77, 1, length.out = 101)
      graphics::lines(b, 2 - 1.5/b, col = "gray20", lwd = 1.2)
      cols <- c("#176B87", "#C4522B", "#6651A1")
      for (j in 1:3) {
        v <- figure_rows(z, "q", c(.5, 1, 4)[j])
        if (nrow(v)) figure_intervals(v$b + c(-.004, 0, .004)[j], v$estimate,
          v$lower, v$upper, cols[j], c(16, 17, 15)[j])
      }
      graphics::legend("topleft", paste("q/p =", c(.5, 1, 4)), col = cols,
        pch = c(16, 17, 15), bty = "n", cex = .65)
      record(fig, "a", if (figure_resolved(z, 9)) "available" else "partial", "Branching survival curve and threshold-40 observations")
      paths <- data$latency_paths; moments <- data$latency; limit <- figure_limit(data)
      valid <- if (nrow(moments)) moments$case_id[moments$complete_run %in% TRUE] else character()
      paths <- figure_rows(paths, "case_id", valid)
      distributions <- if (nrow(paths)) lapply(split(paths$centered_time, paths$case_id),
        latency_distribution_data, shift = limit$shift, survival = limit$survival, lambda = limit$lambda) else list()
      if (!length(distributions)) missing(fig, "b", "Individual centered latency samples are needed for CDF differences.") else {
        ns <- vapply(distributions, `[[`, 0L, "hits")
        bound <- sqrt(-log(.05/2)/2) / sqrt(min(ns))
        graphics::plot(NA, xlim = figure_range(lapply(distributions, function(x) x$cdf$time)),
          ylim = figure_range(-bound, bound, lapply(distributions, function(x) x$cdf$difference)),
          xlab = "Centered time to burden", ylab = "Empirical CDF minus Gumbel CDF")
        graphics::abline(h = c(-bound, 0, bound), col = "gray60", lty = c(3, 2, 3))
        for (j in seq_along(distributions)) graphics::lines(distributions[[j]]$cdf,
          col = cols[j], lty = j, lwd = 1)
        graphics::legend("topright", names(distributions), col = cols[seq_along(distributions)],
          lty = seq_along(distributions), bty = "n", cex = .6)
        record(fig, "b", if (length(distributions) == 3) "available" else "partial", "Asymptotic Kolmogorov bound uses smallest hit count")
      }
      if (!length(distributions)) missing(fig, "c", "Raw latency samples are required; summary moments cannot determine quantiles.") else {
        ns <- vapply(names(distributions), function(id) moments$N[match(id, moments$case_id)], 0.)
        largest <- names(distributions)[which.max(ns)]; qq <- distributions[[largest]]$qq
        graphics::plot(qq$predicted, qq$empirical, pch = 16, cex = .35, col = "#6651A1",
          xlab = "Predicted Gumbel quantile", ylab = "Centered latency quantile")
        graphics::abline(0, 1, lty = 2, col = "gray20")
        figure_note(paste("Largest available N =", max(ns)))
        record(fig, "c", if (max(ns) == 8000) "available" else "partial", "All hit observations, using ppoints() plotting positions")
      }
      z <- moments
      if (!nrow(z)) missing(fig, "d", "Conditional mean latencies and bootstrap endpoints are required.") else {
        graphics::plot(NA, xlim = range(log(z$N)), ylim = figure_range(z$empirical_mean,
          z$mean_lower, z$mean_upper, limit$conditional_mean), xaxt = "n",
          xlab = "System size N", ylab = "Mean centered latency")
        graphics::axis(1, log(z$N), z$N); graphics::abline(h = limit$conditional_mean, lty = 2, col = "gray25")
        for (j in seq_len(nrow(z))) figure_intervals(log(z$N[j]), z$empirical_mean[j],
          z$mean_lower[j], z$mean_upper[j], cols[j], c(16, 17, 15)[j])
        complete <- all(is.finite(z$mean_lower) & is.finite(z$mean_upper))
        if (!complete) figure_note("Mean bootstrap endpoints unavailable")
        record(fig, "d", if (complete && nrow(z) == 3) "available" else "partial", "Conditional means; bootstrap intervals require original hit samples")
      }
    }
    graphics::mtext(paste(fig, if (source == "reported") "- reported summaries and new deterministic calculations" else "- new simulation results"),
      outer = TRUE, side = 3, line = .8, font = 2, cex = .95)
    profile <- if (source == "run" && nrow(data$provenance) && "profile" %in% names(data$provenance))
      paste("Profile:", data$provenance$profile[1], "/") else ""
    graphics::mtext(paste(profile, "Missing inputs and intervals are shown explicitly. Original manuscript artwork is unchanged."),
      outer = TRUE, side = 1, line = 1, cex = .66, col = "gray35")
    if (separate) { grDevices::dev.off(owned_device); owned_device <- NULL }
  }
  if (!separate) { grDevices::dev.off(owned_device); owned_device <- NULL }
  write_csv(bind_rows(coverage), file.path(directory, "coverage.csv"))
  saveRDS(list(source = source, provenance = data$provenance,
    trajectory_settings = data$trajectory_settings, package_version = as.character(utils::packageVersion("stemfeedback")),
    figure_code_fingerprint = code_fingerprint()),
    file.path(directory, "figure_provenance.rds"))
  invisible(files)
}
