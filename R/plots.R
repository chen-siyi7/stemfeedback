plot_results <- function(output, profile = "smoke") {
  opts <- list(out = output, profile = profile)

  figure_file <- function(name) {
    dir.create(file.path(opts$out, "figures"), recursive = TRUE, showWarnings = FALSE)
    file.path(opts$out, "figures", paste0(name, ".pdf"))
  }
  open_figure <- function(name, width = 6.2, height = 4.4) {
    grDevices::pdf(figure_file(name), width = width, height = height, family = "Helvetica",
                   useDingbats = FALSE)
    par(mar = c(4.2, 4.5, 1.8, .8), las = 1, bty = "l", cex.axis = .85, cex.lab = .95)
  }
  fresh_note <- function() mtext(paste("New", opts$profile, "run"), side = 3, adj = 1,
                                  line = .4, cex = .7, col = "gray40")
  summary_file <- file.path(opts$out, "summary", "all_cases.csv")
  if (file.exists(summary_file)) {
    z <- read.csv(summary_file, stringsAsFactors = FALSE)
    z <- z[is.finite(z$estimate) & is.finite(z$diffusion_prediction), , drop = FALSE]
    if (nrow(z)) {
      open_figure("fixation_vs_diffusion")
      plot(z$diffusion_prediction, z$estimate, xlim = c(0, 1), ylim = c(0, 1),
           pch = 16, col = "#264653", xlab = "Diffusion prediction",
           ylab = "Simulated leukemic fixation")
      segments(z$diffusion_prediction, z$lower, z$diffusion_prediction, z$upper,
               col = "#5E7783")
      abline(0, 1, col = "#B45F4D", lwd = 1.2)
      fresh_note(); dev.off()
    }
  }
  therapy_file <- file.path(opts$out, "deterministic", "therapy_trajectories.csv")
  if (file.exists(therapy_file)) {
    z <- read.csv(therapy_file, stringsAsFactors = FALSE)
    z$fraction <- z$Y1 / (z$X1 + z$Y1)
    labels <- c(equal_independent = "Equal depletion", common_division = "Division-linked",
                selective_1 = "Selective 1", selective_2 = "Selective 2")
    cols <- c("#1B4D3E", "#427EA1", "#A3692D", "#A3433F")
    open_figure("therapy_recovery")
    plot(NA, xlim = range(z$time), ylim = c(0, 1), xlab = "Time",
         ylab = "Leukemic share of stem cells")
    abline(v = 3, col = "gray65", lty = 3)
    for (j in seq_along(labels)) {
      v <- z[z$schedule == names(labels)[j], ]
      lines(v$time, v$fraction, col = cols[j], lwd = 1.5)
    }
    legend("right", legend = unname(labels), col = cols, lwd = 1.5,
           bty = "n", cex = .78)
    fresh_note(); dev.off()
  }
  occupation_file <- file.path(opts$out, "summary", "occupation_time.csv")
  if (file.exists(occupation_file)) {
    z <- read.csv(occupation_file, stringsAsFactors = FALSE)
    open_figure("feedback_integral_over_time")
    yr <- range(c(z$integral - 1.96 * z$integral_se,
                  z$integral + 1.96 * z$integral_se), na.rm = TRUE)
    if (!all(is.finite(yr)) || diff(yr) == 0) yr <- range(z$integral, 0)
    plot(NA, xlim = range(z$time), ylim = yr, xlab = "Time",
         ylab = "Mean feedback correction")
    cols <- c("#284B63", "#B45F4D")
    for (j in seq_along(sort(unique(z$noise)))) {
      nu <- sort(unique(z$noise))[j]; v <- z[z$noise == nu, ]
      lines(v$time, v$integral, col = cols[j], lwd = 1.5)
      points(v$time, v$integral, col = cols[j], pch = 16, cex = .65)
    }
    legend("bottomright", legend = paste("Terminal noise", sort(unique(z$noise))),
           col = cols[seq_along(unique(z$noise))], lwd = 1.5, bty = "n", cex = .8)
    fresh_note(); dev.off()
  }
  latency_file <- file.path(opts$out, "deterministic", "latency_shift.csv")
  if (file.exists(latency_file)) {
    z <- read.csv(latency_file, stringsAsFactors = FALSE)
    open_figure("deterministic_latency_shift")
    plot(-log10(z$initial_density), z$deterministic_shift,
         type = "b", pch = 16, col = "#264653", lwd = 1.5,
         xlab = expression(-log[10]("Initial leukemic density")),
         ylab = "Centered crossing time")
    fresh_note(); dev.off()
  }
  invisible(list.files(file.path(output, "figures"), pattern = "[.]pdf$", full.names = TRUE))
}
