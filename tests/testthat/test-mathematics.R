test_that("mathematical and simulation invariants hold", {
  near <- function(x, y, tolerance = 1e-8) {
    expect_true(isTRUE(all.equal(unname(x), unname(y), tolerance = tolerance)),
      info = paste("Numerical check:", paste(x, collapse = ","), "versus", paste(y, collapse = ",")))
  }
  check <- function(condition, why) expect_true(isTRUE(condition), info = why)

  cases <- paper_designs()
  check(length(cases) == 98, "Paper design case count changed.")
  check(sum(vapply(cases, function(x) x$full_attempts, 0)) == 2253900,
        "Paper attempt ledger changed.")

  m <- minimal_model(100, p = 1, q = 4, a = .75, dx = 1, dy = 4)
  eq <- equilibrium(m)
  near(c(eq$K, eq$r, eq$signal), c(.5, 1, 2 / 3))
  z <- coexistence(.2, m)
  near(density_drift(z, m), rep(0, 4), 1e-12)
  near(sum(c(m$q / z[1], -m$p / z[3]) * density_drift(z, m)[c(1, 3)]), 0)

  state <- c(X1 = 23, X2 = 34, Y1 = 16, Y2 = 21)
  s <- feedback(state, m)
  d <- density_drift(state / m$N, m)
  near(d[1], m$p * (2 * m$a * s - 1) * state[1] / m$N)
  near(d[2], (2 * m$p * (1 - m$a * s) * state[1] - state[2]) / m$N)
  noisy <- minimal_model(100, q = 4, dy = 4, noise = 16)
  near(density_drift(state / 100, noisy), d)
  variance <- function(model) crossprod(model$jumps, model$jumps * event_rates(state, model))
  vdelta <- variance(noisy) - variance(m)
  near(diag(vdelta), c(0, 2 * 16 * state[2], 0, 2 * 16 * state[4]))
  near(vdelta - diag(diag(vdelta)), matrix(0, 4, 4))

  for (model in list(m, transit_model(100, rho = .5, dy0 = 4))) {
    point <- coexistence(equilibrium(model)$K * .37, model) * 1.07
    J <- model_jacobian(point, model); h <- 1e-6
    numericJ <- sapply(seq_along(point), function(j) {
      up <- down <- point; up[j] <- up[j] + h; down[j] <- down[j] - h
      (density_drift(up, model) - density_drift(down, model)) / (2 * h)
    })
    near(J, numericJ, 1e-7)
  }

  equal <- minimal_model(100, p = 1, q = 1, dy = 1)
  near(scale_function(equal)$probability(.15), .7, 1e-9)
  near(diffusion_probability(c(20, 18, 30, 21), equal), .6, 1e-9)
  f <- .75
  near(scale_function(m)$probability((1 - f) * eq$K),
       f * (2 * m$p + (m$q - m$p) * f) / (m$p + m$q), 1e-8)
  for (theta in c(0, 2)) {
    x <- minimal_model(200, q = 4, dy = 1, b = .75 + theta / 200)
    left <- minimal_model(200, q = 4, dy = 1 * (1 - 1e-8), b = x$b)
    right <- minimal_model(200, q = 4, dy = 1 * (1 + 1e-8), b = x$b)
    target <- scale_function(x)$probability(equilibrium(x)$K / 2)
    near(scale_function(left)$probability(equilibrium(left)$K / 2), target, 1e-6)
    near(scale_function(right)$probability(equilibrium(right)$K / 2), target, 1e-6)
  }

  small <- minimal_model(20, q = 4, dy = 4)
  ref <- solve_reference(small, c(25, 25), keep_system = TRUE)
  check(ref$residual < 1e-12, "Reference sparse solve residual is too large.")
  check(all(ref$U >= -1e-10 & ref$U <= 1 + 1e-10), "Invalid U boundary solution.")
  check(all(ref$E >= -1e-10 & ref$E <= 1 + 1e-10), "Invalid escape solution.")
  check(all(ref$U + ref$E <= 1 + 1e-9), "Disjoint boundary events exceed one.")
  near(as.matrix(ref$system$A %*% ref$system$solution), ref$system$rhs, 1e-11)
  check(ref$U[1, 6] == 1 && ref$U[26, 6] == 0 && ref$E[26, 6] == 1,
        "Reference boundary values changed.")

  for (model in list(small, minimal_model(20, q = 4, dy = 4, noise = 8))) {
    at <- c(5, 8, 4, 7); rates <- event_rates(at, model)
    x <- at[1]; y <- at[3]; old <- ref$U[x + 1, y + 1]
    new <- model$jumps + matrix(rep(at, each = nrow(model$jumps)), nrow = nrow(model$jumps))
    generator <- sum(rates * (ref$U[cbind(new[, 1] + 1, new[, 3] + 1)] - old))
    near(generator, occupation_kernel(at, model, ref), 1e-9)
  }

  initial <- c(5, 5, 5, 5)
  first <- simulate_path(initial, small, 13, reference = ref, observation_times = c(0, 2))
  second <- simulate_path(initial, small, 13, reference = ref, observation_times = c(0, 2))
  check(identical(first, second), "Path seeds do not reproduce the same trajectory.")
  censor <- simulate_path(initial, small, 13, max_events = 1, observation_times = 0)
  check(censor$censored && !censor$healthy_loss && !censor$hit,
        "Event cap was confused with a biological outcome.")

  residuals <- vapply(seq_len(200), function(i)
    simulate_path(initial, small, 1000L + i, reference = ref,
                  time_limit = 1, observation_times = 0)$residual, 0.)
  limit <- 4 * sd(residuals) / sqrt(length(residuals)) + 1e-10
  check(abs(mean(residuals)) < limit, "Finite-time reference identity failed Monte Carlo check.")
  complete <- lapply(seq_len(400), function(i)
    simulate_path(initial, small, 5000L + i, reference = ref, observation_times = 0))
  end_residuals <- vapply(complete, `[[`, 0., "residual")
  check(abs(mean(end_residuals)) < 4 * sd(end_residuals) / sqrt(length(end_residuals)),
        "Stopped reference identity failed Monte Carlo check.")
  check(all(vapply(complete, function(x) x$status %in%
    c("healthy_loss", "leukemic_loss", "reference_exit"), FALSE)),
    "An uncensored stopped path ended with an invalid status.")
  for (x in complete) {
    near(sum(x$occupation[, "time"]), x$time, 1e-9)
    near(sum(x$occupation[, "integral"]), x$integral, 1e-9)
    near(x$positive_departure_integral + x$negative_departure_integral,
         x$integral, 1e-9)
  }

  serial <- parallel_paths(1:4, function(i) simulate_path(initial, small, 7000L + i,
    observation_times = 0), workers = 1)
  parallel <- parallel_paths(1:4, function(i) simulate_path(initial, small, 7000L + i,
    observation_times = 0), workers = 2)
  check(identical(serial, parallel), "Worker count changed seeded trajectories.")

  extended <- minimal_model(100, q = 4, dy = 4, asymmetric = c(.2, .35),
    proliferation = function(s) 1 + .2 * s,
    proliferation_derivative = function(s) .2)
  e <- equilibrium(extended); zz <- coexistence(.4 * e$K, extended)
  near(density_drift(zz, extended), rep(0, 4), 1e-10)
  grad <- c(4 * (1 - .35) / zz[1], 0,
            -1 * (1 - .2) / zz[3], 0)
  near(sum(grad * density_drift(zz * 1.05, extended)), 0, 1e-10)

  ther <- minimal_model(100, q = 3, dy = 1.5)
  start <- c(X1 = .25, X2 = .25, Y1 = .125, Y2 = .25)
  for (policy in list(c(1, 1, 0, -6), c(0, 0, .6, 0), c(.5, 1.5, 0, 0),
                      c(.5, 2.5, 0, 3))) {
    out <- ode_trajectory(start, ther, c(0, .5, 1, 2, 3, 5, 20),
      treatment = list(start = 0, end = 3, mu = rep(policy[1:2], each = 2),
                       kappa = policy[3]))
    h <- conserved_coordinate(out$X1[nrow(out)], out$Y1[nrow(out)], ther) -
      conserved_coordinate(start[1], start[3], ther)
    near(h, policy[4], 1e-7)
  }
  shift <- latency_shift(minimal_model(2000, q = 1, b = 1, dy = 1),
                         densities = c(1e-4, 1e-5))
  near(shift$deterministic_shift[2], -6.795311, 1e-5)
  near(branching_prediction(minimal_model(2000, q = 1, b = .75), 40)["threshold_probability"],
       1 / 40, 1e-12)
  check(all(closed_routh_margin(minimal_model(100, dx = 1, dy = 1)) > 0),
        "Closed-segment stability test failed.")
})
