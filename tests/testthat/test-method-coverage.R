test_that("neutral mean time solves the backward equation with absorbing endpoints", {
  m <- minimal_model(N = 100, p = 1, q = 1, dx = 1, dy = 1)
  K <- equilibrium(m)$K
  expect_equal(neutral_mean_time(c(0, 1), m), c(0, 0))
  expect_equal(neutral_mean_time(.5, m), 2 * K * m$N * log(2) / m$p)
  f <- c(.1, .3, .5, .7, .9)
  expect_equal(neutral_mean_time(f, m), neutral_mean_time(1 - f, m))
  expect_equal(neutral_mean_time(f, m), m$N * neutral_mean_time(f, m, "slow"))
  h <- 1e-4
  second <- (neutral_mean_time(f + h, m, "slow") -
    2 * neutral_mean_time(f, m, "slow") + neutral_mean_time(f - h, m, "slow")) / h^2
  expect_equal(.5 * m$p * f * (1 - f) / K * second, rep(-1, length(f)),
               tolerance = 1e-6)
  larger <- minimal_model(N = 200, p = 1, q = 1, dx = 1, dy = 1)
  faster <- minimal_model(N = 100, p = 2, q = 2, dx = 2, dy = 2)
  expect_equal(neutral_mean_time(f, larger), 2 * neutral_mean_time(f, m))
  expect_equal(neutral_mean_time(f, faster), neutral_mean_time(f, m) / 2)
  expect_true(all(is.finite(neutral_mean_time(c(1e-15, 1 - 1e-15), m))))
})

test_that("neutral mean time rejects models outside its stated formula", {
  models <- list(minimal_model(q = 1, dy = 1, b = .76),
    minimal_model(q = 2, dy = 2), minimal_model(q = 1, dy = 2),
    minimal_model(q = 1, dy = 1, asymmetric = c(.1, .1)),
    minimal_model(q = 1, dy = 1, proliferation = function(s) 2))
  for (m in models) expect_error(neutral_mean_time(.5, m), "requires")
  m <- minimal_model(q = 1, dy = 1)
  for (f in list(-.1, 1.1, NA_real_, Inf, "half"))
    expect_error(neutral_mean_time(f, m), "fractions")
  expect_error(neutral_mean_time(.5, m, "days"), "arg")
})

test_that("recovery prediction matches the identity and independent ODE recovery", {
  m <- minimal_model(100, q = 3, dy = 1.5)
  initial <- c(.25, .25, .125, .25)
  policies <- list(c(1, 1, 0, -6), c(0, 0, .6, 0),
                   c(.5, 1.5, 0, 0), c(.5, 2.5, 0, 3))
  for (z in policies) {

    mu <- c(z[1], .2, z[2], .4)
    out <- ode_trajectory(initial, m, c(0, 3, 100),
      treatment = list(start = 0, end = 3, mu = mu, kappa = z[3]))
    measured <- c(out$X1[2], out$Y1[2])
    predicted <- recovery_prediction(initial[m$stem], m, 3 * z[1:2],
                                      after = measured, counts = FALSE)
    expect_equal(predicted$predicted_delta_h, z[4])
    expect_equal(predicted$observed_delta_h, z[4], tolerance = 2e-7)
    expect_equal(predicted$discrepancy, 0, tolerance = 2e-7)
    expect_equal(unname(predicted$recovered_stems), c(tail(out$X1, 1), tail(out$Y1, 1)),
                 tolerance = 2e-7)
    expect_equal(predicted$leukemic_fraction, tail(out$Y1 / (out$X1 + out$Y1), 1),
                 tolerance = 2e-7)
    expect_equal(conserved_coordinate(predicted$recovered_stems[1],
      predicted$recovered_stems[2], m), predicted$recovered_h, ignore_attr = TRUE)
    from_counts <- recovery_prediction(100 * initial[m$stem], m, 3 * z[1:2],
                                        after = 100 * measured)
    expect_equal(from_counts, predicted)
  }
  same <- recovery_prediction(initial[m$stem], m, counts = FALSE)
  expect_equal(same$leukemic_fraction, 1 / 3)
  expect_gt(recovery_prediction(initial[m$stem], m, c(3, 3), counts = FALSE)$leukemic_fraction,
             same$leukemic_fraction)
  expect_lt(recovery_prediction(initial[m$stem], m, c(1.5, 7.5), counts = FALSE)$leukemic_fraction,
             same$leukemic_fraction)
})

test_that("recovery prediction handles deeper hierarchies and validates measurements", {
  m <- transit_model(200)
  predicted <- recovery_prediction(c(20, 10), m, c(.3, .6))
  expect_equal(density_drift(predicted$recovered_state, m), rep(0, m$dimension),
               tolerance = 1e-12)
  expect_equal(project_stems(.1, .05, m),
    recovery_prediction(c(.1, .05), m, counts = FALSE)$recovered_stems[1],
    ignore_attr = TRUE)
  for (x in list(c(0, 1), c(-1, 1), c(1, NA), c(1, Inf), 1, c("1", "2")))
    expect_error(recovery_prediction(x, m), "measurements")
  for (x in list(c(-1, 0), c(1, NA), c(1, Inf), 1))
    expect_error(recovery_prediction(c(1, 1), m, exposure = x), "Exposure")
  expect_error(recovery_prediction(c(1, 1), m, after = c(0, 1)), "measurements")
  expect_error(recovery_prediction(c(1, 1), m, counts = NA), "counts")
  expect_error(recovery_prediction(c(1, 1), minimal_model(b = .8)), "requires")
  expect_error(recovery_prediction(c(1, 1), minimal_model(asymmetric = c(.1, 0))), "requires")
  equal <- minimal_model(q = 1, dy = 1)
  for (delta in c(-1000, -10, 0, 10, 1000)) {
    exposure <- c(max(0, -delta), max(0, delta))
    expect_equal(recovery_prediction(c(1, 1), equal, exposure)$leukemic_fraction,
                 stats::plogis(-delta), tolerance = 1e-12)
  }
})

test_that("the installed method map covers every exported function and help topic", {
  map <- utils::read.delim(system.file("METHOD_FUNCTION_MAP.tsv", package = "stemfeedback"),
                           stringsAsFactors = FALSE, na.strings = "")
  split_names <- function(x) unique(unlist(strsplit(x[!is.na(x)], ",", fixed = TRUE)))
  public <- split_names(map$public_functions)
  internal <- split_names(map$internal_functions)
  exports <- getNamespaceExports("stemfeedback")
  expect_setequal(public, exports)
  ns <- asNamespace("stemfeedback")
  for (name in public) {
    expect_true(is.function(getExportedValue("stemfeedback", name)), info = name)
    expect_true(length(utils::help(name, package = "stemfeedback")) > 0, info = name)
  }
  expect_true(all(vapply(internal, exists, FALSE, envir = ns, inherits = FALSE)))
  expect_false(anyDuplicated(map$method_id) > 0L)
})
