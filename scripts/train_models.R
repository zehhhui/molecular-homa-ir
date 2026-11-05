#!/usr/bin/env Rscript
# Fit baseline models, cross-fit integration weights, and apply frozen models at follow-up.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[1] != '--config') {
  stop('Usage: Rscript scripts/train_models.R --config config/default.json')
}
file_arg <- grep('^--file=', commandArgs(FALSE), value = TRUE)
repo <- dirname(dirname(normalizePath(sub('^--file=', '', file_arg[1]))))
source(file.path(repo, 'R', 'model_core.R'))
cfg <- read_model_config(args[2], repo)
run_model_pipeline(cfg)
