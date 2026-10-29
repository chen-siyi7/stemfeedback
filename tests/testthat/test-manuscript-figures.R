test_that("reported transcriptions preserve values and missing inputs", {
  d <- manuscript_figure_data(source = "reported")
  expect_equal(nrow(figure_rows(d$all, "suite", "baseline")), 15)
  expect_equal(nrow(figure_rows(d$all, "suite", "weak_selection")), 5)
  expect_equal(nrow(figure_rows(d$all, "suite", "founders")), 6)
  expect_equal(nrow(figure_rows(d$all, "suite", "establishment")), 11)
  expect_equal(nrow(d$bands), 12)
  expect_equal(d$bands$band_lower[1:6], c(0, .1, .25, .5, .75, .9))
  expect_equal(d$bands$band_upper[1:6], c(.1, .25, .5, .75, .9, 1))
  expect_equal(d$bands$integral[6], -3.1e-6)
  expect_equal(d$decomposition$integral, c(-.0183, -.0692))
  expect_true(all(is.na(d$fluctuations$rms_lower)))
  expect_true(all(is.na(d$latency$mean_lower)))
  expect_equal(nrow(d$latency_paths), 0)
  expect_equal(nrow(d$occupation_time), 4)
  expect_false("integral" %in% names(d$occupation_time))
  expect_equal(nrow(d$provenance), 13)
  expect_equal(reported_numbers("1/4"), .25)
  expect_equal(reported_numbers("$3.5\\times10^{-5}$"), 3.5e-5)
  expect_true(is.na(reported_numbers("---")[1]))
})

test_that("CDF and QQ calculations use the conditional Gumbel distribution", {
  shift <- -6.79531; rho <- .5; lambda <- 1/3
  p <- stats::ppoints(20)
  values <- shift + (log(rho) - log(-log(p))) / lambda
  z <- latency_distribution_data(values, shift, rho, lambda)
  expect_equal(z$qq$empirical, z$qq$predicted, tolerance = 1e-12)
  expect_equal(z$ks, 1/40, tolerance = 1e-12)
  expect_equal(z$hits, 20)
  expect_true(all(abs(z$cdf$difference) <= z$ks + 1e-12))
  tied <- latency_distribution_data(c(0, 0, 1), 0, 1, 1)
  expected <- max(abs(c(0, 2/3, 2/3, 1) - exp(-exp(-c(0, 0, 1, 1)))))
  expect_equal(tied$ks, expected)
  expect_error(latency_distribution_data(c(0, NA), 0, 1, 1))
})

test_that("total decomposition uncertainty retains dependence between bands", {
  records <- data.frame(integral = c(-.1, -.2, -.3, -.5))
  tau <- c(1, 2, 4, 3); ds <- c(.1, .4, -.2, .3); psi <- c(.2, .5, .3, .8)
  for (b in 1:6) {
    records[[paste0("band", b, "_time")]] <- if (b <= 2) tau else 0
    records[[paste0("band", b, "_feedback")]] <- if (b <= 2) ds else 0
    records[[paste0("band", b, "_response")]] <- if (b <= 2) psi else 0
  }
  case <- paper_designs()[["occupation_0"]]
  z <- occupation_decomposition(records, case)
  T <- mean(tau); D <- mean(ds); P <- mean(psi)
  influence <- 2 * (-D*P/T^2*(tau-T) + P/T*(ds-D) + D/T*(psi-P))
  expect_equal(z$product_of_means, 2*D*P/T)
  expect_equal(z$product_se, sd(influence)/2)
  expect_equal(z$covariance_se, sd(records$integral - mean(records$integral) - influence)/2)
  expect_equal(z$product_of_means + z$covariance_contribution, z$integral)
  expect_equal(z$integral_se, sd(records$integral)/2)

  eps <- 1e-6; w <- rep(.25, 4); w[1] <- w[1]+eps; w[2] <- w[2]-eps
  perturbed <- 2*sum(w*ds)*sum(w*psi)/sum(w*tau)
  expect_equal((perturbed-z$product_of_means)/eps, influence[1]-influence[2], tolerance = 1e-6)
})

test_that("latency mean endpoints follow the manuscript bootstrap method", {
  case <- paper_designs()[["latency_500"]]
  records <- data.frame(status = rep("burden_hit", 5), time = c(1, 2, 4, 9, 30),
    censored = FALSE, attempt = 1:5, seed = 11:15)
  z <- summarize_latency(records, case, -6.79531, bootstrap = 211, seed = 44)
  expected <- bootstrap_interval(records$time - 3*log(500), mean, B = 211, seed = 44)
  expect_equal(z$moments$mean_lower, unname(expected["lower"]))
  expect_equal(z$moments$mean_upper, unname(expected["upper"]))
  expect_equal(z$moments$mean_bootstrap_resamples, 211)
  expect_equal(z$moments$mean_interval_method, "percentile_bootstrap")
})

test_that("local return rates match every plotted manuscript abscissa", {
  out <- tempfile("figure-theory-"); on.exit(unlink(out, recursive = TRUE), add = TRUE)
  z <- prepare_manuscript_theory(out)$rates
  reported <- manuscript_figure_data(source = "reported")$rates
  key <- function(x) paste(x$N, x$family, x$speed)
  expected <- reported$local_return_rate[match(key(z), key(reported))]
  expect_equal(nrow(z), 24)
  expect_true(all(is.finite(expected)))
  expect_true(all(abs(z$local_return_rate - expected) <= .00050001))
})

test_that("deterministic trajectories record full states, weights and models", {
  designs <- manuscript_trajectory_designs()
  expect_length(designs, 6)
  expect_true(all(vapply(designs, function(x) length(x$initial) == x$model$dimension, FALSE)))
  expect_true(all(vapply(designs, function(x) grepl("Illustrative", x$provenance), FALSE)))
  expect_error(validate_trajectory_designs(list(list(initial = 1))), "named")
  out <- tempfile("figure-ode-"); on.exit(unlink(out, recursive = TRUE), add = TRUE)
  write_manuscript_deterministic(out, designs, profile = "smoke")
  stored <- readRDS(file.path(out, "deterministic", "trajectory_designs.rds"))
  expect_identical(stored, designs)
  d <- manuscript_figure_data(out, source = "reported")
  expect_equal(length(unique(d$trajectories$trajectory_id)), 6)
  expect_equal(nrow(d$trajectory_settings), 33)
  expect_equal(unique(d$therapy$schedule),
    c("equal_independent", "common_division", "selective_1", "selective_2"))
  at3 <- d$therapy[d$therapy$time == 3, ]
  f <- at3$Y1 / (at3$X1 + at3$Y1)
  expect_equal(f[at3$schedule == "equal_independent"], .6056196, tolerance = 1e-6)
  expect_equal(f[at3$schedule == "selective_1"], .06770003, tolerance = 1e-6)
  expect_true(all(read.csv(file.path(out, "deterministic", "therapy_summary.csv"))$descendant_killing == 0))
  expect_equal(d$latency_limit$conditional_mean, -7.1431, tolerance = 1e-4)
})

test_that("all seven reported figures render without inventing missing data", {
  out <- tempfile("reported-figures-"); on.exit(unlink(out, recursive = TRUE), add = TRUE)
  before <- grDevices::dev.cur()
  files <- plot_manuscript_figures(out, source = "reported")
  expect_length(files, 1)
  expect_true(file.exists(files))
  expect_identical(readChar(files, 4, useBytes = TRUE), "%PDF")
  expect_identical(grDevices::dev.cur(), before)
  coverage <- read.csv(file.path(out, "manuscript_figures", "coverage.csv"))
  expect_equal(nrow(coverage), 22)
  expect_equal(length(unique(coverage$figure)), 7)
  missing <- paste(coverage$figure[coverage$status == "missing"], coverage$panel[coverage$status == "missing"])
  expect_setequal(missing, c("Fig1 a", "Fig1 b", "Fig6 a", "FigS1 b", "FigS1 c"))
  expect_true(all(coverage$status[coverage$figure == "Fig5" & coverage$panel == "a"] == "partial"))
  expect_error(plot_manuscript_figures(out, source = "reported", figures = "Fig8"))
})

test_that("new-run figures draw raw latency and outcome data", {

  d <- manuscript_figure_data(source = "reported"); d$source <- "run"
  d$provenance <- data.frame(profile = "synthetic unit test fixtures")
  d$latency_limit <- figure_limit(d)
  p <- stats::ppoints(20)
  centered <- d$latency_limit$shift +
    (log(.5) - log(-log(p))) / (1/3)
  d$latency_paths <- bind_rows(lapply(c(500, 2000, 8000), function(N)
    data.frame(case_id = paste0("latency_", N), centered_time = centered)))
  d$latency$mean_lower <- d$latency$empirical_mean-.2
  d$latency$mean_upper <- d$latency$empirical_mean+.2
  d$occupation_paths <- list("0" = data.frame(time = 1:4, healthy_loss = c(0, 1, 0, 1), censored = FALSE),
                            "16" = data.frame(time = 1:4, healthy_loss = c(1, 0, 0, 0), censored = FALSE))
  d$bands$integral_se <- 0
  out <- tempfile("new-figures-"); on.exit(unlink(out, recursive = TRUE), add = TRUE)
  expect_no_warning(files <- plot_manuscript_figures(out, source = "run",
    figures = c("Fig6", "FigS1"), data = d, separate = TRUE))
  expect_length(files, 2)
  expect_true(all(file.exists(files)))
  coverage <- read.csv(file.path(out, "manuscript_figures", "coverage.csv"))
  expect_equal(coverage$status[coverage$figure == "Fig6" & coverage$panel == "d"], "available")
  expect_true(all(coverage$status[coverage$figure == "FigS1"] == "available"))
})
test_that("decomposition axes include every component and total interval", {
  z <- data.frame(product_of_means = c(.01, .02),
    covariance_contribution = c(-.03, -.04), integral = c(-.02, -.02),
    product_se = c(.04, .03), covariance_se = c(.02, .08), integral_se = c(.1, .01))
  limits <- figure_decomposition_range(z)
  expect_lt(limits[1], -.216)
  expect_gt(limits[2], .176)
  z$product_se <- z$covariance_se <- z$integral_se <- NA_real_
  expect_true(all(is.finite(figure_decomposition_range(z))))
})
