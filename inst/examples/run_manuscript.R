library(stemfeedback)
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("Usage: Rscript run_manuscript.R OUTPUT [PROFILE] [WORKERS] [ATTEMPTS]")
output <- args[1]
profile <- if (length(args) >= 2) args[2] else "smoke"
workers <- if (length(args) >= 3) as.integer(args[3]) else 1L
attempts <- if (length(args) >= 4) as.integer(args[4]) else NULL
run_manuscript(output, profile = profile, workers = workers, replicates = attempts)
