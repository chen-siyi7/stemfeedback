recovery_prediction <- function(stems, model, exposure = c(0, 0),
                                after = NULL, counts = TRUE) {
  if (abs(model$b - model$a) > 1e-14 || any(model$asymmetric != 0))
    stop("Recovery prediction requires equal stem self-renewal and symmetric stem division.")
  if (!is.logical(counts) || length(counts) != 1L || is.na(counts))
    stop("counts must be TRUE or FALSE.")
  positive_pair <- function(x) is.numeric(x) && length(x) == 2L &&
    all(is.finite(x)) && all(x > 0)
  if (!positive_pair(stems) || (!is.null(after) && !positive_pair(after)))
    stop("Stem measurements must be two finite positive values, healthy then leukemic.")
  if (!is.numeric(exposure) || length(exposure) != 2L ||
      any(!is.finite(exposure)) || any(exposure < 0))
    stop("Exposure must contain two finite nonnegative integrated killing rates.")
  initial <- unname(stems) / if (counts) model$N else 1
  h0 <- conserved_coordinate(initial[1], initial[2], model)
  delta <- unname(model$p * exposure[2] - model$q * exposure[1])
  h1 <- h0 + delta
  xi <- project_coordinate(h1, model)
  state <- coexistence(xi, model)
  recovered <- state[model$stem]
  ans <- list(initial_h = h0, predicted_delta_h = delta, recovered_h = h1,
    recovered_stems = recovered, recovered_state = state,
    leukemic_fraction = unname(recovered[2] / sum(recovered)))
  if (!is.null(after)) {
    measured <- unname(after) / if (counts) model$N else 1
    observed <- conserved_coordinate(measured[1], measured[2], model) - h0
    ans$observed_delta_h <- observed
    ans$discrepancy <- observed - delta
  }
  ans
}
