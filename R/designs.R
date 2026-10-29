paper_designs <- function() {
  cases <- list(); index <- 0L
  add <- function(id, suite, model, initial, attempts, rule = "first_stem_loss",
                  threshold = NULL, burden = NULL, note = "Paper design") {
    index <<- index + 1L
    cases[[id]] <<- list(id = id, suite = suite, index = index, model = model,
      initial = initial, full_attempts = attempts, stop_rule = rule,
      threshold = threshold, burden = burden, note = note)
  }

  balanced <- data.frame(q = c(.25, 1, rep(4, 6)), N = c(200, 200, 50, 100, 200, 400, 800, 1600),
                         runs = c(2000, 2000, 2000, 2000, 2000, 2000, 1600, 2000))
  for (i in seq_len(nrow(balanced))) {
    z <- balanced[i, ]; m <- minimal_model(z$N, q = z$q)
    add(paste0("balanced_", i), "baseline", m, start_from_stems(round(z$N / 4), round(z$N / 4), m), z$runs)
  }
  clearance <- data.frame(q = c(.25, 4, 4, 4), N = c(200, 200, 800, 3200), runs = c(2000, 2000, 1600, 2000))
  for (i in seq_len(nrow(clearance))) {
    z <- clearance[i, ]; m <- minimal_model(z$N, q = z$q, dy = 1)
    add(paste0("equal_clearance_", i), "baseline", m,
        start_from_stems(round(z$N / 4), round(z$N / (4 * z$q)), m), z$runs)
  }
  for (q in c(.25, 1, 4)) {
    m <- minimal_model(200, q = q)
    add(paste0("rare_", q), "baseline", m, start_from_stems(95, 5, m), 6000)
  }
  weak <- data.frame(q = c(1, 1, 1, 1, 4), theta = c(-2, -1, 1, 2, 1))
  for (i in seq_len(nrow(weak))) {
    m <- minimal_model(200, q = weak$q[i], b = .75 + weak$theta[i] / 200)
    add(paste0("selection_", i), "weak_selection", m, start_from_stems(50, 50, m), 2000)
  }
  speeds <- c(.25, .5, 1, 2, 4, 8)
  for (N in c(200, 50)) for (speed in speeds) for (family in c("clearance", "transit")) {
    dy0 <- if (N == 200) 1 else 4
    m <- if (family == "clearance") minimal_model(N, dx = speed, dy = speed * dy0, k = speed) else
      transit_model(N, speed, dy0)
    initial <- start_from_stems(if (N == 200) 50 else 12, 12, m)
    id <- paste(family, N, speed, sep = "_")
    add(id, family, m, initial, 10000)
    if (speed %in% c(.25, 1, 8)) add(paste0("mechanism_", id), "mechanism", m, initial, 1500)
  }
  for (r in c(1, .25)) for (N in c(100, 200, 400)) {
    m <- minimal_model(N, dy = 4 * r)
    add(paste("founder", r, N, sep = "_"), "founders", m, founder_start(m), 100000)
  }
  for (r in c(1, .25)) for (N in c(100, 200)) {
    m <- minimal_model(N, dy = 4 * r)
    add(paste("initialization", r, N, sep = "_"), "initialization", m,
        founder_start(m, descendants = TRUE), 100000,
        note = "New explicit sensitivity grid; the original grid was not reported")
  }
  thresholds <- data.frame(N = c(400, 800, 1600, 1600), L = c(10, 20, 10, 40))
  for (r in c(1, .25)) for (i in seq_len(nrow(thresholds))) {
    N <- thresholds$N[i]; L <- thresholds$L[i]; m <- minimal_model(N, dy = 4 * r)
    add(paste("entry", r, N, L, sep = "_"), "entry", m, founder_start(m), 100000,
        "founder_entry", threshold = L)
  }
  noise <- data.frame(N = c(50, 50, 100, 100, 100, 100, 200, 200),
                      nu = c(0, 4, 0, 1, 4, 16, 0, 4),
                      runs = c(10000, 10000, 10000, 5000, 10000, 5000, 20000, 20000))
  for (i in seq_len(nrow(noise))) {
    m <- minimal_model(noise$N[i], noise = noise$nu[i]); initial <- rep(round(m$N / 4), 4)
    add(paste("noise", m$N, m$noise, sep = "_"), "noise", m, initial, noise$runs[i])
  }
  for (nu in c(0, 16)) {
    m <- minimal_model(100, noise = nu)
    add(paste0("occupation_", nu), "occupation", m, rep(25, 4), 10000)
  }
  establishment <- expand.grid(b = c(.8, .9, 1), q = c(.5, 1, 4))
  for (i in seq_len(nrow(establishment))) {
    m <- minimal_model(2000, q = establishment$q[i], b = establishment$b[i])
    add(paste0("establishment_", i), "establishment", m, founder_start(m), 2500,
        "establishment", threshold = 40)
  }
  for (NL in list(c(2000, 20), c(8000, 80))) {
    m <- minimal_model(NL[1], q = 1, b = .8)
    add(paste0("establishment_L", NL[2]), "establishment", m, founder_start(m), 2500,
        "establishment", threshold = NL[2])
  }
  for (N in c(500, 2000, 8000)) {
    m <- minimal_model(N, q = 1, b = 1)
    add(paste0("latency_", N), "latency", m, founder_start(m), 2400,
        "burden", burden = ceiling(.05 * N))
  }
  cases
}

select_designs <- function(cases, suite = "all", profile = "smoke") {
  if (suite == "staged") return(cases[vapply(cases, function(c)
    c$suite == "entry" && c$model$N %in% c(800, 1600) && c$threshold == c$model$N / 40, FALSE)])
  if (suite != "all") cases <- cases[vapply(cases, function(c) c$suite == suite, FALSE)]
  if (profile != "smoke") return(cases)

  keep <- vapply(cases, function(c) {
    m <- c$model
    switch(c$suite,
      baseline = c$id %in% c("balanced_1", "balanced_2", "balanced_3", "rare_1"),
      weak_selection = c$id %in% c("selection_1", "selection_4"),
      clearance = m$N == 50 && m$healthy$clearance %in% c(.25, 1, 8),
      transit = m$N == 50 && m$healthy$division[2] %in% c(.25, 1, 8),
      mechanism = m$N == 50 && (if (m$n == 2) m$healthy$clearance else m$healthy$division[2]) %in% c(.25, 8),
      founders = m$N == 100,
      initialization = m$N == 100,
      entry = m$N == 800,
      noise = c$id %in% c("noise_50_0", "noise_50_4", "noise_100_16"),
      occupation = TRUE,
      establishment = c$id %in% c("establishment_4", "establishment_5", "establishment_6"),
      latency = m$N == 500, FALSE)
  }, FALSE)
  cases[keep]
}

case_attempts <- function(case, profile, override = NULL) {
  if (!is.null(override)) return(as.integer(override))
  if (profile == "full") return(case$full_attempts)
  if (profile == "quick") return(min(100L, case$full_attempts))
  if (case$suite %in% c("founders", "initialization", "entry")) 10L else 3L
}

design_table <- function(cases, profile, override = NULL) {
  bind_rows(lapply(cases, function(c) {
    m <- c$model; e <- equilibrium(m)
    row <- data.frame(case_id = c$id, case_index = c$index, suite = c$suite,
      profile = profile, attempts = case_attempts(c, profile, override),
      paper_attempts = c$full_attempts, N = m$N, p = m$p, q = m$q, a = m$a, b = m$b,
      k = m$k, dx = m$healthy$clearance, dy = m$leukemic$clearance, noise = m$noise,
      healthy_compartments = m$n, leukemic_compartments = m$m, K = e$K, r = e$r,
      threshold = if (is.null(c$threshold)) NA else c$threshold,
      burden = if (is.null(c$burden)) NA else c$burden, note = c$note)
    for (i in seq_along(c$initial)) row[[paste0("initial_", m$names[i])]] <- c$initial[i]
    row
  }))
}
