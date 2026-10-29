test_that("installed API handles simulation, references and workers", {
  model <- minimal_model(20)
  initial <- c(5, 5, 5, 5)
  serial <- simulate_competition(model, initial, n = 4, seed = 7000,
                                 time_limit = 1, workers = 1)
  parallel <- simulate_competition(model, initial, n = 4, seed = 7000,
                                   time_limit = 1, workers = 2)
  expect_s3_class(serial, "stemfeedback_result")
  expect_identical(serial$records, parallel$records)
  expect_equal(serial$summary$attempts, 4)
  expect_null(serial$paths)
  expect_output(print(serial), "Stem-cell competition")
  ref <- solve_reference(model, c(25, 25))
  out <- simulate_competition(model, initial, n = 2, seed = 81,
                              reference = ref, time_limit = .2, keep_paths = TRUE)
  expect_length(out$paths, 2)
  expect_true(all(is.finite(out$records$integral)))
  expect_true(all(is.finite(out$records$residual)))
  expect_equal(out$feedback$reference_initial, ref$U[6, 6])
  expect_error(simulate_competition(model, initial, n = 0))
})

test_that("experiment resume preserves records and protects incompatible runs", {
  out <- tempfile("stemfeedback-api-")
  on.exit(unlink(out, recursive = TRUE), add = TRUE)
  expect_message(run_experiments(out, suite = "founders", replicates = 1,
                                 max_events = 1, plots = FALSE), "Running")
  files <- list.files(out, recursive = TRUE, full.names = TRUE)
  before <- tools::md5sum(files)
  expect_error(run_experiments(out, suite = "founders", replicates = 1,
                               max_events = 1, seed = 1234, plots = FALSE), "No files were changed")
  expect_identical(before, tools::md5sum(files))
  path <- file.path(out, "paths", "founder_1_100.csv")
  old <- read.csv(path)
  suppressMessages(run_experiments(out, suite = "founders", replicates = 2,
                                   max_events = 1, plots = FALSE))
  new <- read.csv(path)
  expect_equal(nrow(new), 2)
  expect_identical(old, new[1, , drop = FALSE])
  note <- file.path(out, "author-note.txt")
  writeLines("Keep this file", note)
  suppressMessages(run_experiments(out, suite = "founders", replicates = 1,
                                   max_events = 1, seed = 1234, force = TRUE, plots = FALSE))
  expect_true(file.exists(note))
  manifest <- read.csv(file.path(out, "run_manifest.csv"))
  expect_equal(manifest$base_seed, 1234)
  expect_equal(manifest$package_version, as.character(utils::packageVersion("stemfeedback")))
  expect_true(nzchar(manifest$code_fingerprint))
})

test_that("plotting uses package arguments and closes its device", {
  out <- tempfile("stemfeedback-plot-")
  dir.create(file.path(out, "summary"), recursive = TRUE)
  on.exit(unlink(out, recursive = TRUE), add = TRUE)
  write.csv(data.frame(estimate = c(.3, .6), lower = c(.2, .5),
    upper = c(.4, .7), diffusion_prediction = c(.35, .65)),
    file.path(out, "summary", "all_cases.csv"), row.names = FALSE)
  before <- grDevices::dev.cur()
  files <- plot_results(out)
  expect_length(files, 1)
  expect_true(file.exists(files))
  expect_identical(readChar(files, 4, useBytes = TRUE), "%PDF")
  expect_identical(grDevices::dev.cur(), before)
})

test_that("socket workers can load the installed package", {
  sockets <- parallel_paths
  platform <- .Platform; platform$OS.type <- "windows"
  environment(sockets) <- list2env(list(.Platform = platform),
                                  parent = environment(parallel_paths))
  result <- sockets(1:2, function(i) stemfeedback::minimal_model(N = 10 + i)$N,
                     workers = 2)
  expect_equal(unlist(result), c(11, 12))
})
