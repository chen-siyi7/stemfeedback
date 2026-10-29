solve_reference <- function(model, limits = c(300L, 300L), keep_system = FALSE) {
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Install the Matrix package.")
  if (!is.null(model$proliferation)) stop("The occupation calculation uses unregulated stem division rates.")
  stopifnot(length(limits) == 2, all(limits >= 2), all(limits == round(limits)))
  Lx <- limits[1]; Ly <- limits[2]; nx <- Lx - 1L; ny <- Ly - 1L
  x <- rep(seq_len(nx), ny); y <- rep(seq_len(ny), each = nx)
  id <- seq_along(x); eq <- equilibrium(model)
  s <- 1 / (1 + model$k * (eq$A * x + eq$B * y) / model$N)
  p <- model$p * (1 - model$asymmetric[1]); q <- model$q * (1 - model$asymmetric[2])
  total <- p * x + q * y
  pxp <- p * x * model$a * s / total; pxm <- p * x * (1 - model$a * s) / total
  pyp <- q * y * model$b * s / total; pym <- q * y * (1 - model$b * s) / total
  rows <- id; cols <- id; values <- rep(1, length(id))
  add <- function(ok, target, probability) {
    rows <<- c(rows, id[ok]); cols <<- c(cols, target[ok]); values <<- c(values, -probability[ok])
  }
  add(x < nx, id + 1L, pxp); add(x > 1, id - 1L, pxm)
  add(y < ny, id + nx, pyp); add(y > 1, id - nx, pym)
  A <- Matrix::sparseMatrix(i = rows, j = cols, x = values,
                            dims = c(length(id), length(id)))
  rhsU <- ifelse(x == 1, pxm, 0)
  rhsE <- ifelse(x == nx, pxp, 0) + ifelse(y == ny, pyp, 0)
  solution <- as.matrix(Matrix::solve(A, cbind(rhsU, rhsE)))
  if (any(!is.finite(solution)) || any(solution < -1e-10) ||
      any(solution > 1 + 1e-10) || any(rowSums(solution) > 1 + 1e-10))
    stop("Reference solve left the probability simplex beyond numerical tolerance.")
  raw_solution <- solution
  solution[] <- pmin(1, pmax(0, solution))
  solution[, 2] <- pmin(solution[, 2], 1 - solution[, 1])
  roundoff_adjustment <- max(abs(solution - raw_solution))
  residual <- max(abs(as.matrix(A %*% solution) - cbind(rhsU, rhsE)))
  U <- E <- matrix(0, Lx + 1L, Ly + 1L)
  U[1, 2:(Ly + 1L)] <- 1
  U[2:Lx, 2:Ly] <- solution[, 1]
  E[Lx + 1L, 2:(Ly + 1L)] <- 1; E[2:(Lx + 1L), Ly + 1L] <- 1
  E[2:Lx, 2:Ly] <- solution[, 2]
  Psi <- B <- matrix(0, Lx + 1L, Ly + 1L)
  where <- cbind(x + 1L, y + 1L)
  Psi[where] <- p * model$a * x * (U[cbind(x + 2L, y + 1L)] - U[cbind(x, y + 1L)]) +
    q * model$b * y * (U[cbind(x + 1L, y + 2L)] - U[cbind(x + 1L, y)])
  B[where] <- p * model$a * x * (E[cbind(x + 2L, y + 1L)] + E[cbind(x, y + 1L)]) +
    q * model$b * y * (E[cbind(x + 1L, y + 2L)] + E[cbind(x + 1L, y)])
  out <- list(U = U, E = E, Psi = Psi, error_response = B, limits = limits,
              A_feedback = eq$A, B_feedback = eq$B, residual = residual,
              probability_roundoff_adjustment = roundoff_adjustment,
              signature = reference_signature(model))
  if (keep_system) out$system <- list(A = A, rhs = cbind(rhsU, rhsE), solution = solution)
  structure(out, class = "stemfeedback_reference")
}

reference_signature <- function(model) {
  if (!is.null(model$proliferation))
    stop("The feedback reference does not support proliferation multipliers.")
  eq <- equilibrium(model)
  c(p = model$p * (1 - model$asymmetric[1]),
    q = model$q * (1 - model$asymmetric[2]), a = model$a, b = model$b,
    feedback_x = model$k * eq$A / model$N,
    feedback_y = model$k * eq$B / model$N)
}

validate_reference <- function(reference, model) {
  expected <- reference_signature(model)
  if (!is.list(reference) || is.null(reference$signature))
    stop("Reference metadata are missing. Rebuild it with solve_reference().")
  if (!isTRUE(all.equal(reference$signature, expected, tolerance = 1e-12)))
    stop("Reference is incompatible with this model. Recompute it with solve_reference(model).")
  limits <- reference$limits
  if (length(limits) != 2L || any(!is.finite(limits)) ||
      any(limits < 2) || any(limits != round(limits)))
    stop("Invalid reference domain.")
  for (name in c("U", "E", "Psi", "error_response"))
    if (!is.matrix(reference[[name]]) ||
        !identical(dim(reference[[name]]), as.integer(limits + 1)))
      stop("Invalid reference matrix: ", name)
  invisible(TRUE)
}

reference_values <- function(x, y, reference, model) {
  validate_reference(reference, model)
  stopifnot(length(x) == 1L, length(y) == 1L, is.finite(x), is.finite(y),
            x == round(x), y == round(y))
  reference_values_unchecked(x, y, reference)
}

reference_values_unchecked <- function(x, y, reference) {
  if (x < 0 || y < 0 || x > reference$limits[1] || y > reference$limits[2])
    stop("Stem counts lie outside the reference domain.")
  i <- x + 1L; j <- y + 1L
  c(U = reference$U[i, j], escape = reference$E[i, j], Psi = reference$Psi[i, j],
    error_response = reference$error_response[i, j],
    signal = unname(1 / (1 + reference$signature["feedback_x"] * x +
                          reference$signature["feedback_y"] * y)))
}

occupation_kernel <- function(state, model, reference) {
  validate_state(state, model)
  validate_reference(reference, model)
  x <- state[model$stem[1]]; y <- state[model$stem[2]]
  if (x == 0 || y == 0) return(0)
  r <- reference_values_unchecked(x, y, reference)
  unname((feedback(state, model) - r["signal"]) * r["Psi"])
}
