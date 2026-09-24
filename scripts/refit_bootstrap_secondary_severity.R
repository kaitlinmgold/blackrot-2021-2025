#!/usr/bin/env Rscript
# Calibrate secondary severity group and interaction tests.
# Inputs: model file, output directory, filename prefix, simulations, cores, seed.

suppressPackageStartupMessages({
  library(glmmTMB)
  library(jsonlite)
  library(parallel)
})

args <- commandArgs(trailingOnly = TRUE)
model_file <- if (length(args) >= 1) args[[1]] else stop("model RDS required")
out_dir <- if (length(args) >= 2) args[[2]] else stop("output directory required")
prefix <- if (length(args) >= 3) args[[3]] else "secondary"
n_sim <- if (length(args) >= 4) as.integer(args[[4]]) else 1999L
n_cores <- if (length(args) >= 5) as.integer(args[[5]]) else 4L
seed <- if (length(args) >= 6) as.integer(args[[6]]) else 20260841L
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

bundle <- readRDS(model_file)
data <- bundle$data
models <- bundle$severity

ctrl_default <- glmmTMBControl(
  optCtrl = list(iter.max = 20000L, eval.max = 20000L)
)
ctrl_bfgs <- glmmTMBControl(
  optimizer = optim,
  optArgs = list(method = "BFGS"),
  optCtrl = list(maxit = 20000L)
)

usable <- function(fit) {
  !inherits(fit, "try-error") &&
    isTRUE(fit$sdr$pdHess) &&
    identical(unname(fit$fit$convergence), 0L) &&
    is.finite(as.numeric(logLik(fit)))
}

fit_retry <- function(formula, data_i) {
  fit <- try(suppressWarnings(glmmTMB(
    formula,
    dispformula = ~Year,
    data = data_i,
    family = gaussian(),
    REML = FALSE,
    control = ctrl_default
  )), silent = TRUE)
  if (usable(fit)) return(list(fit = fit, optimizer = "nlminb"))
  fit <- try(suppressWarnings(glmmTMB(
    formula,
    dispformula = ~Year,
    data = data_i,
    family = gaussian(),
    REML = FALSE,
    control = ctrl_bfgs
  )), silent = TRUE)
  if (usable(fit)) return(list(fit = fit, optimizer = "BFGS"))
  list(fit = NULL, optimizer = NA_character_)
}

run_test <- function(null_fit, alternative_fit, null_formula, alternative_formula,
                     label, stem, seed_offset) {
  observed <- max(0, 2 * (
    as.numeric(logLik(alternative_fit)) - as.numeric(logLik(null_fit))
  ))
  set.seed(seed + seed_offset)
  simulations <- simulate(null_fit, nsim = n_sim, seed = seed + seed_offset)
  one <- function(index) {
    data_i <- data
    data_i$HA_Logit_Severity <- simulations[[index]]
    null_result <- fit_retry(null_formula, data_i)
    alt_result <- fit_retry(alternative_formula, data_i)
    if (is.null(null_result$fit) || is.null(alt_result$fit)) {
      return(data.frame(
        Simulation = index,
        Statistic = NA_real_,
        Null_Optimizer = null_result$optimizer,
        Alternative_Optimizer = alt_result$optimizer,
        Status = "fit_failure"
      ))
    }
    statistic <- max(0, 2 * (
      as.numeric(logLik(alt_result$fit)) - as.numeric(logLik(null_result$fit))
    ))
    data.frame(
      Simulation = index,
      Statistic = statistic,
      Null_Optimizer = null_result$optimizer,
      Alternative_Optimizer = alt_result$optimizer,
      Status = "ok"
    )
  }
  started <- proc.time()[[3]]
  rows <- if (.Platform$OS.type == "windows" || n_cores <= 1L) {
    lapply(seq_len(n_sim), one)
  } else {
    mclapply(seq_len(n_sim), one, mc.cores = n_cores, mc.preschedule = TRUE)
  }
  table_out <- do.call(rbind, rows)
  elapsed <- proc.time()[[3]] - started
  valid <- table_out$Statistic[is.finite(table_out$Statistic)]
  exceedances <- sum(valid >= observed)
  p_bootstrap <- (exceedances + 1) / (length(valid) + 1)
  summary <- list(
    test = label,
    observed_likelihood_ratio = observed,
    calibration = paste(
      "Parametric-bootstrap likelihood-ratio test; both nested models were",
      "refitted by ML in every simulation with year-specific residual dispersion"
    ),
    requested_simulations = n_sim,
    valid_simulations = length(valid),
    failed_simulations = n_sim - length(valid),
    exceedances = exceedances,
    bootstrap_p = p_bootstrap,
    monte_carlo_se = sqrt(p_bootstrap * (1 - p_bootstrap) / (length(valid) + 1)),
    seed = seed + seed_offset,
    elapsed_seconds = elapsed,
    status_counts = as.list(table(table_out$Status)),
    quantiles = as.list(quantile(valid, c(0.5, 0.9, 0.95, 0.99)))
  )
  write_json(
    summary,
    file.path(out_dir, paste0(prefix, "_", stem, ".json")),
    pretty = TRUE,
    auto_unbox = TRUE,
    digits = 15
  )
  write.csv(
    table_out,
    file.path(out_dir, paste0(prefix, "_", stem, "_statistics.csv")),
    row.names = FALSE
  )
  summary
}

year_formula <- formula(models$year_ml)
additive_formula <- formula(models$additive_ml)
full_formula <- formula(models$full_ml)

message("Bootstrapping secondary severity group effect...")
group_summary <- run_test(
  models$year_ml,
  models$additive_ml,
  year_formula,
  additive_formula,
  "Secondary analysis-group effect on transformed cluster severity",
  "severity_group_lrt_refit_bootstrap",
  1L
)
message("Bootstrapping secondary severity interaction...")
interaction_summary <- run_test(
  models$additive_ml,
  models$full_ml,
  additive_formula,
  full_formula,
  "Secondary year-by-analysis-group interaction for transformed cluster severity",
  "severity_interaction_lrt_refit_bootstrap",
  2L
)

combined <- list(group = group_summary, interaction = interaction_summary)
write_json(
  combined,
  file.path(out_dir, paste0(prefix, "_severity_refit_bootstrap_summary.json")),
  pretty = TRUE,
  auto_unbox = TRUE,
  digits = 15
)
cat(toJSON(combined, pretty = TRUE, auto_unbox = TRUE, digits = 15), "\n")
