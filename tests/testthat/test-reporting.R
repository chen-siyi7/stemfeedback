test_that("intervals, censoring and last-cell identities are correct", {
  check <- function(ok, why) expect_true(isTRUE(ok), info = why)
  near <- function(a, b, tol = 1e-10) check(max(abs(a - b)) < tol, "Numerical regression failed")

  near(unname(bootstrap_interval(7, B = 100)), c(7, 7, 7))

  zero <- cluster_score_interval(rep(0, 10))
  near(zero$upper, 1 - .025^(1 / 10))
  check(zero$lower == 0 && zero$upper > 0 && zero$bootstrap_resamples == 0,
        "Zero-success founder interval lost its uncertainty")
  one <- cluster_score_interval(rep(1, 10))
  near(one$lower, .025^(1 / 10))
  half <- cluster_score_interval(rep(.5, 10))
  check(half$lower < .5 && half$upper > .5, "Constant interior scores imply false certainty")
  mixed <- cluster_score_interval(c(rep(0, 90), rep(.5, 8), 1, 1), seed = 19)
  check(mixed$lower < mixed$upper && mixed$bootstrap_resamples == 2000,
        "Nondegenerate original-attempt bootstrap changed method")

  case <- paper_designs()[["latency_500"]]
  rec <- data.frame(status = rep("burden_hit", 4), time = c(13, 14, 16, 18),
                    censored = FALSE, attempt = 1:4, seed = 1:4)
  report <- summarize_latency(rec, case, -6.79531)
  near(report$moments$mean_crossing_time, mean(rec$time))
  near(report$moments$empirical_variance, var(rec$time))
  check(is.finite(report$moments$variance_lower) && is.finite(report$moments$variance_upper) &&
        report$moments$variance_bootstrap_resamples == 1500,
        "Latency variance bootstrap interval is missing")
  rec$status[4] <- "event_limit"; rec$censored[4] <- TRUE
  limited <- summarize_latency(rec, case, -6.79531)
  check(all(is.na(limited$quantiles$empirical)) &&
        is.na(limited$moments$empirical_mean) && is.na(limited$moments$empirical_variance) &&
        limited$moments$moment_status == "incomplete_censored",
        "Censored latency runs were reported as a complete conditional sample")

  m <- minimal_model(20, noise = 16)
  ref <- solve_reference(m, c(25, 25))
  for (stems in list(c(1, 1), c(1, 8), c(9, 1))) {
    at <- c(stems[1], 7, stems[2], 11)
    next_states <- sweep(m$jumps, 2, at, "+")
    delta <- ref$U[cbind(next_states[, 1] + 1, next_states[, 3] + 1)] -
      ref$U[at[1] + 1, at[3] + 1]
    near(sum(event_rates(at, m) * delta), occupation_kernel(at, m, ref), 1e-9)
  }

  m <- cell_model(30, healthy = lineage(c(2, 3), c(.75, .4), .7, c(.1, .7, 1)),
    leukemic = lineage(2, .75, 4, c(.5, 1)), noise = 8)
  for (stems in list(c(1, 1), c(1, 8), c(9, 1), c(12, 7))) {
    at <- c(stems[1], 3, 9, stems[2], 7)
    next_states <- sweep(m$jumps, 2, at, "+")
    fraction <- next_states[, 4] / (next_states[, 1] + next_states[, 4])
    near(sum(event_rates(at, m) * (fraction - at[4] / (at[1] + at[4]))), 0, 1e-10)
  }

})
