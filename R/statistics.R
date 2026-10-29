wilson <- function(successes, attempts, level = .95) {
  if (attempts < 1) return(c(estimate = NA, lower = NA, upper = NA))
  stopifnot(successes >= 0, successes <= attempts)
  z <- qnorm((1 + level) / 2); p <- successes / attempts; d <- 1 + z^2 / attempts
  centre <- (p + z^2 / (2 * attempts)) / d
  half <- z * sqrt(p * (1 - p) / attempts + z^2 / (4 * attempts^2)) / d
  c(estimate = p, lower = max(0, centre - half), upper = min(1, centre + half))
}

mean_interval <- function(x, level = .95) {
  x <- x[is.finite(x)]; n <- length(x)
  if (!n) return(c(mean = NA, se = NA, lower = NA, upper = NA, n = 0))
  m <- mean(x); se <- if (n > 1) sd(x) / sqrt(n) else NA_real_
  d <- qnorm((1 + level) / 2) * se
  c(mean = m, se = se, lower = m - d, upper = m + d, n = n)
}

bootstrap_interval <- function(x, statistic = mean, B = 2000L, seed = 1L) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(estimate = NA, lower = NA, upper = NA))
  set.seed(seed)
  samples <- replicate(B, statistic(x[sample.int(length(x), length(x), replace = TRUE)]))
  c(estimate = statistic(x), lower = unname(quantile(samples, .025)),
    upper = unname(quantile(samples, .975)))
}

cluster_score_interval <- function(scores, B = 2000L, seed = 1L, level = .95) {
  stopifnot(length(scores) > 0, all(scores %in% c(0, .5, 1)),
            B >= 1, B == round(B), level > 0, level < 1)
  n <- length(scores); estimate <- mean(scores); alpha <- 1 - level

  if (all(scores == 0)) {
    lower <- 0; upper <- -expm1(log(alpha / 2) / n)
    method <- "exact_nonzero_probability_bound"; used <- 0L
  } else if (all(scores == 1)) {
    lower <- exp(log(alpha / 2) / n); upper <- 1
    method <- "exact_nonunit_probability_bound"; used <- 0L
  } else if (length(unique(scores)) == 1L) {
    half <- sqrt(log(2 / alpha) / (2 * n))
    lower <- max(0, estimate - half); upper <- min(1, estimate + half)
    method <- "bounded_mean_hoeffding"; used <- 0L
  } else {
    set.seed(seed)
    counts <- tabulate(2L * scores + 1L, nbins = 3L)
    resamples <- rmultinom(B, n, counts / n)
    draws <- (resamples[2, ] / 2 + resamples[3, ]) / n
    limits <- quantile(draws, c(alpha / 2, 1 - alpha / 2), names = FALSE)
    lower <- limits[1]; upper <- limits[2]
    method <- "original_attempt_percentile_bootstrap"; used <- B
  }
  list(estimate = estimate, lower = lower, upper = upper,
       method = method, bootstrap_resamples = used)
}

sample_skewness <- function(x) {
  if (length(x) < 3 || var(x) == 0) return(NA_real_)
  mean((x - mean(x))^3) / mean((x - mean(x))^2)^1.5
}

independent_contrast <- function(a, b) {
  stopifnot(length(a) > 1, length(b) > 1)
  difference <- mean(b) - mean(a)
  se <- sqrt(var(a) / length(a) + var(b) / length(b))
  c(difference = difference, se = se, lower = difference - 1.96 * se,
    upper = difference + 1.96 * se)
}

write_csv <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  write.csv(x, file, row.names = FALSE, na = "")
}

bind_rows <- function(xs) {
  xs <- Filter(function(x) !is.null(x) && nrow(x) > 0, xs)
  if (!length(xs)) return(data.frame())
  keys <- unique(unlist(lapply(xs, names)))
  do.call(rbind, lapply(xs, function(x) {
    for (key in setdiff(keys, names(x))) x[[key]] <- NA
    x[keys]
  }))
}
