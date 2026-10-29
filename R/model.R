lineage <- function(division = 1, self_renewal = 0.75, clearance = 1,
                    weights = c(rep(0, length(division)), 1)) {
  stopifnot(length(division) == length(self_renewal), length(division) >= 1,
            all(is.finite(division)), all(division > 0),
            all(is.finite(self_renewal)), all(self_renewal >= 0 & self_renewal <= 1),
            length(clearance) == 1, is.finite(clearance), clearance > 0,
            length(weights) == length(division) + 1, all(is.finite(weights)),
            all(weights >= 0), any(weights > 0))
  list(division = division, self_renewal = self_renewal,
       clearance = clearance, weights = weights, n = length(weights))
}

cell_model <- function(N = 100, healthy = lineage(), leukemic = lineage(4, .75, 4),
                       k = 1, noise = 0, asymmetric = c(0, 0),
                       proliferation = NULL, proliferation_derivative = NULL) {
  stopifnot(length(N) == 1, is.finite(N), N > 0, is.finite(k), k > 0,
            is.finite(noise), noise >= 0, length(asymmetric) == 2,
            all(is.finite(asymmetric)), all(asymmetric >= 0 & asymmetric < 1))
  a <- healthy$self_renewal[1]
  stopifnot(a > .5, leukemic$self_renewal[1] > 0,
            all(healthy$self_renewal[-1] < a),
            all(leukemic$self_renewal[-1] < a))
  n <- healthy$n; m <- leukemic$n; d <- n + m
  jumps <- list(); sources <- integer(); bases <- numeric(); fractions <- numeric()
  kinds <- character(); stem_events <- logical(); event_names <- character()
  add <- function(jump, source, base, fraction, kind, stem, label) {
    jumps[[length(jumps) + 1L]] <<- jump
    sources <<- c(sources, source); bases <<- c(bases, base)
    fractions <<- c(fractions, fraction); kinds <<- c(kinds, kind)
    stem_events <<- c(stem_events, stem); event_names <<- c(event_names, label)
  }
  for (ell in 1:2) {
    z <- if (ell == 1) healthy else leukemic
    offset <- if (ell == 1) 0L else n
    prefix <- if (ell == 1) "X" else "Y"
    for (i in seq_along(z$division)) {
      j <- offset + i; stem <- i == 1L
      alpha <- if (stem) asymmetric[ell] else 0
      jump <- numeric(d); jump[j] <- 1
      add(jump, j, z$division[i] * (1 - alpha), z$self_renewal[i],
          "renew", stem, paste0(prefix, i, "_renew"))
      jump[j] <- -1; jump[j + 1L] <- 2
      add(jump, j, z$division[i] * (1 - alpha), z$self_renewal[i],
          "advance", stem, paste0(prefix, i, "_differentiate"))
      if (stem && alpha > 0) {
        jump <- numeric(d); jump[j + 1L] <- 1
        add(jump, j, z$division[i] * alpha, 0, "constant", TRUE,
            paste0(prefix, "1_asymmetric"))
      }
    }
    j <- offset + z$n; jump <- numeric(d); jump[j] <- -1
    add(jump, j, z$clearance, 0, "constant", FALSE, paste0(prefix, "_clearance"))
    if (noise > 0) {
      jump[j] <- 1
      add(jump, j, noise, 0, "constant", FALSE, paste0(prefix, "_reservoir_gain"))
      jump[j] <- -1
      add(jump, j, noise, 0, "constant", FALSE, paste0(prefix, "_reservoir_loss"))
    }
  }
  structure(list(N = N, healthy = healthy, leukemic = leukemic, n = n, m = m,
                 dimension = d, stem = c(1L, n + 1L), terminal = c(n, d),
                 names = c(paste0("X", seq_len(n)), paste0("Y", seq_len(m))),
                 p = healthy$division[1], q = leukemic$division[1], a = a,
                 b = leukemic$self_renewal[1], k = k, noise = noise,
                 weights = c(healthy$weights, leukemic$weights),
                 asymmetric = asymmetric, proliferation = proliferation,
                 proliferation_derivative = proliferation_derivative,
                 jumps = do.call(rbind, jumps), sources = sources, bases = bases,
                 fractions = fractions, kinds = kinds, stem_events = stem_events,
                 event_names = event_names), class = "cell_model")
}

minimal_model <- function(N = 100, p = 1, q = 4, a = .75, b = a,
                          dx = 1, dy = q, k = 1, noise = 0, ...) {
  cell_model(N, lineage(p, a, dx), lineage(q, b, dy), k, noise, ...)
}

transit_model <- function(N = 200, rho = 1, dy0 = 1) {
  stopifnot(is.finite(rho), rho > 0)
  cell_model(N, lineage(c(1, rho), c(.75, 0), 2),
             lineage(c(4, 4 * rho), c(.75, 0), 2 * dy0))
}

feedback <- function(state, model, density = FALSE) {
  1 / (1 + model$k * sum(model$weights * state) / if (density) 1 else model$N)
}

event_rates <- function(state, model, density = FALSE) {
  s <- feedback(state, model, density)
  modifier <- rep(1, length(model$bases))
  i <- model$kinds == "renew"; modifier[i] <- model$fractions[i] * s
  i <- model$kinds == "advance"; modifier[i] <- 1 - model$fractions[i] * s
  rates <- model$bases * state[model$sources] * modifier
  if (!is.null(model$proliferation)) {
    g <- model$proliferation(s)
    if (length(g) != 1 || !is.finite(g) || g <= 0) stop("Proliferation must be positive.")
    rates[model$stem_events] <- rates[model$stem_events] * g
  }
  rates
}

profile <- function(z, a, stem_proliferation = 1) {
  s <- 1 / (2 * a); c <- rep(1, z$n)
  if (z$n > 2) for (i in 2:(z$n - 1L)) {
    upstream <- z$division[i - 1L] * if (i == 2L) stem_proliferation else 1
    upstream_fraction <- if (i == 2L) a else z$self_renewal[i - 1L]
    c[i] <- c[i - 1L] * 2 * upstream * (1 - upstream_fraction * s) /
      (z$division[i] * (1 - 2 * z$self_renewal[i] * s))
  }
  j <- z$n - 1L
  upstream <- z$division[j] * if (j == 1L) stem_proliferation else 1
  fraction <- if (j == 1L) a else z$self_renewal[j]
  c[z$n] <- c[j] * 2 * upstream * (1 - fraction * s) / z$clearance
  c
}

equilibrium <- function(model) {
  g <- if (is.null(model$proliferation)) 1 else model$proliferation(1 / (2 * model$a))
  c <- profile(model$healthy, model$a, g); e <- profile(model$leukemic, model$a, g)
  A <- sum(model$healthy$weights * c); B <- sum(model$leukemic$weights * e)
  M <- (2 * model$a - 1) / model$k
  list(c = c, e = e, A = A, B = B, M = M, K = M / A, r = A / B,
       signal = 1 / (2 * model$a))
}

coexistence <- function(xi, model, counts = FALSE) {
  eq <- equilibrium(model); stopifnot(xi >= 0, xi <= eq$K)
  z <- c(xi * eq$c, eq$r * (eq$K - xi) * eq$e)
  if (counts) z <- round(model$N * z)
  setNames(z, model$names)
}

start_from_stems <- function(x, y, model) {
  eq <- equilibrium(model)
  setNames(round(c(x * eq$c, y * eq$e)), model$names)
}

founder_start <- function(model, founders = 1, descendants = FALSE) {
  eq <- equilibrium(model); z <- round(c(model$N * eq$K * eq$c, rep(0, model$m)))
  z[model$stem[2]] <- founders
  if (descendants) z[model$n + seq_len(model$m)] <- round(founders * eq$e)
  setNames(z, model$names)
}

validate_state <- function(state, model, integer = TRUE) {
  if (length(state) != model$dimension || any(!is.finite(state)) || any(state < 0) ||
      (integer && any(state != round(state)))) stop("Invalid compartment state.")
  invisible(TRUE)
}

density_drift <- function(z, model) drop(crossprod(event_rates(z, model, TRUE), model$jumps))

model_jacobian <- function(z, model) {
  s <- feedback(z, model, TRUE); ds <- -model$k * s^2 * model$weights
  modifiers <- rep(1, length(model$bases)); slopes <- rep(0, length(modifiers))
  i <- model$kinds == "renew"; modifiers[i] <- model$fractions[i] * s; slopes[i] <- model$fractions[i]
  i <- model$kinds == "advance"; modifiers[i] <- 1 - model$fractions[i] * s; slopes[i] <- -model$fractions[i]
  if (!is.null(model$proliferation)) {
    if (is.null(model$proliferation_derivative)) stop("Supply proliferation_derivative for the Jacobian.")
    i <- model$stem_events; g <- model$proliferation(s); gp <- model$proliferation_derivative(s)
    slopes[i] <- slopes[i] * g + modifiers[i] * gp; modifiers[i] <- modifiers[i] * g
  }
  D <- outer(model$bases * z[model$sources] * slopes, ds)
  for (i in seq_along(model$bases)) D[i, model$sources[i]] <-
    D[i, model$sources[i]] + model$bases[i] * modifiers[i]
  crossprod(model$jumps, D)
}
