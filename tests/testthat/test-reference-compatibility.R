test_that("incompatible references fail before an integral is reported", {
  model <- minimal_model(20)
  ref <- solve_reference(model, c(25, 25))
  state <- c(5, 8, 4, 7)
  changed <- minimal_model(20, q = 1, dy = 1)
  expect_error(occupation_kernel(state, changed, ref), "incompatible")
  expect_error(reference_values(5, 4, ref, changed), "incompatible")
  expect_error(simulate_path(state, changed, 1, reference = ref), "incompatible")
  expect_error(simulate_competition(changed, state, n = 1, reference = ref), "incompatible")
  changed <- minimal_model(20, b = .76)
  expect_error(simulate_path(state, changed, 1, reference = ref), "incompatible")
  extended <- minimal_model(20, proliferation = function(s) 1 + .2 * s,
                             proliferation_derivative = function(s) .2)
  expect_error(simulate_path(state, extended, 1, reference = ref), "proliferation")
  ref$signature <- NULL
  expect_error(reference_values(5, 4, ref, model), "metadata")
})

test_that("compatible reference reuse preserves the exact generator identity", {
  model <- minimal_model(20)
  ref <- solve_reference(model, c(25, 25))

  models <- list(minimal_model(20, noise = 16),
                 minimal_model(20, dx = 2, dy = 8, k = 2))
  for (m in models) for (stems in list(c(1, 1), c(1, 8), c(9, 1), c(5, 4))) {
    state <- c(stems[1], 8, stems[2], 7)
    next_states <- sweep(m$jumps, 2, state, "+")
    difference <- ref$U[cbind(next_states[, 1] + 1, next_states[, 3] + 1)] -
      ref$U[state[1] + 1, state[3] + 1]
    generator <- sum(event_rates(state, m) * difference)
    expect_equal(occupation_kernel(state, m, ref), generator, tolerance = 1e-9)
    expect_equal(unname(reference_values(state[1], state[3], ref, m)["signal"]),
      1 / (1 + (state[1] + state[3]) / 20), tolerance = 1e-12)
  }
  larger <- solve_reference(model, c(30, 30))
  path <- simulate_path(c(5, 8, 4, 7), model, 1, reference = ref,
                         enlarged_reference = larger, time_limit = .1)
  expect_true(is.finite(path$integral))
  expect_error(simulate_path(c(5, 8, 4, 7), model, 1, reference = larger,
                              enlarged_reference = ref), "contain")
  expect_error(simulate_path(c(0, 2, 0, 2), model, 1), "stem pool")
})
test_that("reference probabilities and error bounds stay nonnegative after roundoff", {
  m <- minimal_model(200, q = 4, dy = 1)
  ref <- solve_reference(m, c(600L, 150L), keep_system = TRUE)
  expect_true(all(ref$U >= 0 & ref$U <= 1))
  expect_true(all(ref$E >= 0 & ref$E <= 1))
  expect_true(all(ref$U + ref$E <= 1))
  expect_true(all(ref$error_response >= 0))
  expect_gte(ref$probability_roundoff_adjustment, 0)
  expect_lte(ref$probability_roundoff_adjustment, 1e-10)
  expect_lt(ref$residual, 1e-10)
  expect_equal(ref$residual, max(abs(as.matrix(ref$system$A %*%
    ref$system$solution) - ref$system$rhs)))
})
