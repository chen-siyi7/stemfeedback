require_ode <- function() {
  if (!requireNamespace("deSolve", quietly = TRUE))
    stop('Install deSolve: install.packages("deSolve").')
}

therapy_drift <- function(z, model, mu = rep(0, model$dimension), kappa = 0) {
  stopifnot(length(mu) == model$dimension, all(mu >= 0), kappa >= 0, kappa <= 1)
  if (kappa > 0 && any(model$asymmetric > 0))
    stop("The cycle-specific treatment calculation assumes symmetric stem division.")
  dz <- density_drift(z, model) - mu * z
  if (kappa > 0) {
    s <- feedback(z, model, TRUE)
    g <- if (is.null(model$proliferation)) 1 else model$proliferation(s)
    for (ell in 1:2) {
      i <- model$stem[ell]; rate <- c(model$p, model$q)[ell] * g
      fraction <- c(model$a, model$b)[ell]
      dz[i] <- dz[i] - 2 * kappa * rate * fraction * s * z[i]
      dz[i + 1L] <- dz[i + 1L] - 2 * kappa * rate * (1 - fraction * s) * z[i]
    }
  }
  dz
}

ode_trajectory <- function(initial, model, times = seq(0, 100, .1),
                           treatment = NULL, rtol = 1e-10, atol = 1e-12) {
  require_ode(); validate_state(initial, model, FALSE)
  stopifnot(length(times) >= 2, all(diff(times) > 0))
  cuts <- range(times)
  if (!is.null(treatment)) cuts <- sort(unique(c(cuts, treatment$start, treatment$end)))
  cuts <- cuts[cuts >= min(times) & cuts <= max(times)]
  output <- list(); z <- setNames(initial, model$names)
  for (j in seq_len(length(cuts) - 1L)) {
    from <- cuts[j]; to <- cuts[j + 1L]; mu <- rep(0, model$dimension); kappa <- 0
    if (!is.null(treatment) && (from + to) / 2 >= treatment$start && (from + to) / 2 < treatment$end) {
      mu <- treatment$mu; kappa <- treatment$kappa
    }
    ts <- sort(unique(c(from, times[times > from & times < to], to)))
    f <- function(t, state, parms) list(therapy_drift(state, model, mu, kappa))
    out <- as.data.frame(deSolve::ode(z, ts, f, parms = NULL, rtol = rtol, atol = atol,
                                        method = "lsoda", maxsteps = 100000L))
    z <- unlist(out[nrow(out), -1], use.names = TRUE)
    if (j > 1) out <- out[-1, , drop = FALSE]
    output[[j]] <- out
  }
  ans <- do.call(rbind, output); rownames(ans) <- NULL
  if (min(as.matrix(ans[-1])) < -100 * atol) stop("ODE produced a negative cell density.")
  ans
}

invasion_profile <- function(model) {
  stopifnot(all(model$asymmetric == 0), is.null(model$proliferation))
  s <- 1 / (2 * model$a); lambda <- model$q * (model$b / model$a - 1)
  v <- rep(1, model$m); z <- model$leukemic
  if (model$m > 2) for (j in 2:(model$m - 1L))
    v[j] <- 2 * z$division[j - 1L] * (1 - z$self_renewal[j - 1L] * s) * v[j - 1L] /
      (lambda + z$division[j] * (1 - 2 * z$self_renewal[j] * s))
  j <- model$m - 1L
  v[model$m] <- 2 * z$division[j] * (1 - z$self_renewal[j] * s) * v[j] /
    (lambda + z$clearance)
  v
}

latency_shift <- function(model, burden = .05, densities = 10^-(4:7),
                          rtol = 1e-11, atol = 1e-14) {
  require_ode(); prediction <- branching_prediction(model)
  lambda <- unname(prediction["lambda"]); stopifnot(lambda > 0, burden > 0)
  eq <- equilibrium(model); v <- invasion_profile(model); rows <- list()
  for (eta in densities) {
    z <- setNames(c(eq$K * eq$c, eta * v), model$names)
    stopifnot(tail(z, 1) < burden)
    end <- max(100, 2 * abs(log(eta)) / lambda)
    f <- function(t, state, parms) list(density_drift(state, model))
    root <- function(t, state, parms) state[model$dimension] - burden
    event <- function(t, state, parms) state
    out <- deSolve::ode(z, c(0, end), f, parms = NULL, method = "lsodar", rootfunc = root,
                        events = list(func = event, root = TRUE, terminalroot = 1),
                        rtol = rtol, atol = atol, maxsteps = 100000L)
    crossing <- attr(out, "troot")
    if (!length(crossing)) stop("No deterministic burden crossing found.")
    t <- crossing[1]
    rows[[length(rows) + 1L]] <- data.frame(initial_density = eta, burden = burden,
      crossing_time = t, lambda = lambda, deterministic_shift = t + log(eta) / lambda)
  }
  ans <- do.call(rbind, rows)
  ans$change_from_previous <- c(NA, diff(ans$deterministic_shift)); ans
}

return_rate <- function(xi, model) {
  if (abs(model$b - model$a) > 1e-14) stop("Return rates are evaluated at equal stem self-renewal.")
  values <- eigen(model_jacobian(coexistence(xi, model), model), only.values = TRUE)$values
  tangent <- which.min(Mod(values))
  if (Mod(values[tangent]) > 1e-7 * max(1, Mod(values))) stop("No zero tangent eigenvalue found.")
  -max(Re(values[-tangent]))
}

minimum_return_rate <- function(model, intervals = 64L) {
  K <- equilibrium(model)$K; grid <- seq(0, K, length.out = intervals + 1L)
  values <- vapply(grid, return_rate, 0., model = model)
  minima <- vapply(seq_len(intervals), function(i)
    optimize(return_rate, c(grid[i], grid[i + 1L]), model = model, tol = 1e-10)$objective, 0.)
  min(values, minima)
}

routh_margin <- function(model, psi_x) {
  stopifnot(model$n == 2, model$m == 2)
  a <- model$a; dx <- model$healthy$clearance; dy <- model$leukemic$clearance
  psi_y <- 2 * a - 1 - psi_x
  (dx * (2 * a - psi_x) + dy * (1 + psi_x)) *
    (dx * dy + psi_x * dx * model$p + psi_y * dy * model$q) -
    2 * a * dx * dy * (psi_x * model$p + psi_y * model$q)
}

closed_routh_margin <- function(model) {
  a <- model$a; dx <- model$healthy$clearance; dy <- model$leukemic$clearance
  L0 <- 2 * a * dx + dy; L1 <- dy - dx
  R0 <- dx * dy + (2 * a - 1) * dy * model$q; R1 <- dx * model$p - dy * model$q
  candidates <- c(0, 2 * a - 1)
  if (L1 * R1 > 0) {
    vertex <- -(L0 * R1 + L1 * R0 - 2 * a * dx * dy * (model$p - model$q)) / (2 * L1 * R1)
    if (vertex > 0 && vertex < 2 * a - 1) candidates <- c(candidates, vertex)
  }
  min(vapply(candidates, function(x) routh_margin(model, x), 0.))
}
