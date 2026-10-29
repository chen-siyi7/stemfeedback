softplus <- function(x) pmax(x, 0) + log1p(exp(-abs(x)))

conserved_coordinate <- function(x, y, model) {
  stopifnot(x > 0, y > 0)
  model$q * (1 - model$asymmetric[2]) * log(x) -
    model$p * (1 - model$asymmetric[1]) * log(y)
}

project_stems <- function(x, y, model) {
  eq <- equilibrium(model)
  if (x == 0 && y == 0) stop("The empty stem state has no selected coexistence mixture.")
  if (x == 0) return(0)
  if (y == 0) return(eq$K)
  project_coordinate(conserved_coordinate(x, y, model), model)
}

project_coordinate <- function(target, model) {
  eq <- equilibrium(model)
  p <- model$p * (1 - model$asymmetric[1]); q <- model$q * (1 - model$asymmetric[2])
  f <- function(t) q * (log(eq$K) - softplus(-t)) -
    p * (log(eq$r * eq$K) - softplus(t)) - target
  bound <- (abs(target - q * log(eq$K) + p * log(eq$r * eq$K)) +
              (p + q) * log(2)) / min(p, q) + 1
  eq$K * plogis(uniroot(f, c(-bound, bound), tol = 1e-12)$root)
}

diffusion_coefficients <- function(xi, model, theta = model$N * (model$b - model$a)) {
  eq <- equilibrium(model); eta <- eq$r * (eq$K - xi)
  stopifnot(xi > 0, eta > 0)
  p <- model$p * (1 - model$asymmetric[1]); q <- model$q * (1 - model$asymmetric[2])
  g <- if (is.null(model$proliferation)) 1 else model$proliferation(eq$signal)
  c(drift_h = g * p * q * (.5 * (1 / eta - 1 / xi) - theta / model$a),
    variance_h = g * p * q * (q / xi + p / eta))
}

scale_function <- function(model, theta = model$N * (model$b - model$a), tol = 1e-10) {
  eq <- equilibrium(model); K <- eq$K; r <- eq$r
  p <- model$p * (1 - model$asymmetric[1]); q <- model$q * (1 - model$asymmetric[2])
  d0 <- q * r * K; d1 <- p - q * r
  c0 <- (p + q) * (r - 1) + (2 * theta * r / model$a) * q * K
  c1 <- (2 * theta * r / model$a) * (p - q)
  log_derivative <- function(x) {
    t <- d1 * x / d0; f0 <- f1 <- numeric(length(x)); small <- abs(t) < 1e-4
    if (any(small)) {
      ts <- t[small]
      f0[small] <- x[small] / d0 * vapply(ts, function(z) sum((-z)^(0:7) / (1:8)), 0.)
      f1[small] <- x[small]^2 / d0 * vapply(ts, function(z) sum((-z)^(0:7) / (2:9)), 0.)
    }
    if (any(!small)) {
      f0[!small] <- log1p(t[!small]) / d1
      f1[!small] <- (t[!small] - log1p(t[!small])) * d0 / d1^2
    }
    log(q * K + (p - q) * x) + c0 * f0 + c1 * f1
  }
  offset <- max(log_derivative(c(0, K)), -optimize(function(x) -log_derivative(x), c(0, K))$objective)
  scaled <- function(x) exp(log_derivative(x) - offset)
  integral <- function(x) if (x <= 0) 0 else integrate(scaled, 0, x, rel.tol = tol,
                                                      subdivisions = 1000L)$value
  denominator <- integral(K)
  probability <- function(xi) vapply(xi, function(x) {
    stopifnot(x >= 0, x <= K)
    min(1, max(0, 1 - integral(x) / denominator))
  }, 0.)
  list(K = K, r = r, theta = theta, log_derivative = log_derivative,
       probability = probability, founder_coefficient = scaled(K) / (r * denominator),
       rare_clone_bias = K * scaled(K) / (r * denominator))
}

diffusion_probability <- function(state, model) {
  scale_function(model)$probability(project_stems(state[1] / model$N,
                                                state[model$stem[2]] / model$N, model))
}

branching_prediction <- function(model, threshold = Inf) {
  q <- model$q * (1 - model$asymmetric[2])
  g <- if (is.null(model$proliferation)) 1 else model$proliferation(1 / (2 * model$a))
  beta <- q * g * model$b / (2 * model$a)
  delta <- q * g * (1 - model$b / (2 * model$a)); lambda <- beta - delta
  survival <- max(0, 1 - delta / beta)
  if (is.infinite(threshold)) hitting <- survival else {
    stopifnot(threshold >= 1, threshold == round(threshold))
    ratio <- delta / beta
    hitting <- if (abs(ratio - 1) < 1e-12) 1 / threshold else
      if (ratio < 1) (1 - ratio) / (-expm1(threshold * log(ratio))) else {
        z <- threshold * log(ratio)
        exp(log(ratio - 1) - z - log1p(-exp(-z)))
      }
  }
  c(beta = beta, delta = delta, lambda = lambda, survival = survival, threshold_probability = hitting)
}

terminal_to_stem <- function(g, model) {
  eq <- equilibrium(model); output_ratio <- tail(eq$c, 1) / tail(eq$e, 1)
  output_ratio * g / (1 - g + output_ratio * g)
}

neutral_mean_time <- function(f, model, units = c("physical", "slow")) {
  units <- match.arg(units)
  if (!is.numeric(f) || any(!is.finite(f)) || any(f < 0 | f > 1))
    stop("Stem fractions must be finite numbers between zero and one.")
  eq <- equilibrium(model)
  if (abs(model$b - model$a) > 1e-14 ||
      abs(model$p - model$q) > 1e-12 * max(model$p, model$q) ||
      abs(eq$r - 1) > 1e-12 || any(model$asymmetric != 0) ||
      !is.null(model$proliferation))
    stop("Neutral mean time requires b=a, p=q, r=1 and unmodified symmetric stem division.")
  entropy <- numeric(length(f))
  inside <- f > 0 & f < 1
  entropy[inside] <- -f[inside] * log(f[inside]) -
    (1 - f[inside]) * log1p(-f[inside])
  2 * eq$K / model$p * entropy * if (units == "physical") model$N else 1
}
