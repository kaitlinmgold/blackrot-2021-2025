#!/usr/bin/env Rscript
# Calibrate a secondary incidence interaction using the efficient score.
# Inputs: model file, output directory, filename prefix, simulations, cores, seed.

suppressPackageStartupMessages({
  library(glmmTMB)
  library(jsonlite)
  library(parallel)
})

args <- commandArgs(trailingOnly = TRUE)
model_file <- if (length(args) >= 1) args[[1]] else stop("model RDS required")
out_dir <- if (length(args) >= 2) args[[2]] else stop("output directory required")
prefix <- if (length(args) >= 3) args[[3]] else "incidence"
n_sim <- if (length(args) >= 4) as.integer(args[[4]]) else 9999L
n_cores <- if (length(args) >= 5) as.integer(args[[5]]) else 4L
seed <- if (length(args) >= 6) as.integer(args[[6]]) else 20260831L
data_dir <- if (length(args) >= 7) args[[7]] else NULL
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

bundle <- readRDS(model_file)
models <- bundle$incidence
inc_add <- models$additive
inc_full <- models$full

if (!is.null(bundle$data)) {
  data <- bundle$data
  analysis_data <- data
} else {
  if (is.null(data_dir)) stop("data directory is required for a bundle without embedded data")
  d <- read.csv(
    file.path(data_dir, "clean_plot_data_public.csv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  d$Year <- factor(d$Year, levels = 2021:2025)
  d$Program_Class <- factor(
    d$Program_Class,
    levels = c("Untreated", "Biological only", "Mixed", "Conventional only")
  )
  d$BlockYear <- factor(d$BlockYear)
  d$Public_Program_ID <- factor(d$Public_Program_ID)
  d$Cluster_Negative_Count <- d$N_Clusters_Rated - d$Cluster_Positive_Count
  analysis_data <- d
}

n_add_beta <- length(fixef(inc_add)$cond)
n_full_beta <- length(fixef(inc_full)$cond)
n_interaction <- n_full_beta - n_add_beta
interaction_index <- (n_add_beta + 1):n_full_beta

insert_interaction_zeros <- function(additive_parameters) {
  c(
    additive_parameters[seq_len(n_add_beta)],
    rep(0, n_interaction),
    additive_parameters[(n_add_beta + 1):length(additive_parameters)]
  )
}

score_from_objective <- function(objective, null_parameters) {
  gradient <- objective$gr(null_parameters)
  hessian <- optimHess(null_parameters, objective$fn, objective$gr)
  nuisance <- setdiff(seq_along(null_parameters), interaction_index)
  info_nn <- hessian[nuisance, nuisance, drop = FALSE]
  info_in <- hessian[interaction_index, nuisance, drop = FALSE]
  info_ii <- hessian[interaction_index, interaction_index, drop = FALSE]
  efficient_information <- info_ii - info_in %*% solve(info_nn, t(info_in))
  efficient_score <- gradient[interaction_index] -
    info_in %*% solve(info_nn, gradient[nuisance])
  statistic <- as.numeric(
    t(efficient_score) %*% solve(efficient_information, efficient_score)
  )
  list(
    statistic = statistic,
    max_abs_nuisance_score = max(abs(gradient[nuisance])),
    min_eigenvalue = min(
      eigen(efficient_information, symmetric = TRUE, only.values = TRUE)$values
    ),
    condition_number = kappa(efficient_information)
  )
}

ctrl <- glmmTMBControl(optCtrl = list(iter.max = 20000L, eval.max = 20000L))
ctrl_bfgs <- glmmTMBControl(
  optimizer = optim,
  optArgs = list(method = "BFGS"),
  optCtrl = list(maxit = 20000L)
)
full_structure <- glmmTMB(
  formula = formula(inc_full),
  data = analysis_data,
  family = betabinomial(link = "logit"),
  control = ctrl,
  doFit = FALSE
)

observed_null <- insert_interaction_zeros(inc_add$fit$par)
stopifnot(length(observed_null) == length(inc_full$fit$par))
observed <- score_from_objective(inc_full$obj, observed_null)

set.seed(seed)
simulations <- simulate(inc_add, nsim = n_sim, seed = seed)

score_one <- function(index) {
  simulated <- simulations[[index]]
  result <- tryCatch({
    refitted_null <- refit(inc_add, simulated)
    if (!isTRUE(refitted_null$sdr$pdHess) || refitted_null$fit$convergence != 0) {
      retry_data <- analysis_data
      retry_data$Cluster_Positive_Count <- simulated[, 1]
      retry_data$Cluster_Negative_Count <- simulated[, 2]
      refitted_null <- suppressWarnings(glmmTMB(
        formula = formula(inc_add),
        data = retry_data,
        family = betabinomial(link = "logit"),
        control = ctrl_bfgs
      ))
      if (!isTRUE(refitted_null$sdr$pdHess) || refitted_null$fit$convergence != 0) {
        return(c(statistic = NA_real_, nuisance = NA_real_, eigen = NA_real_, condition = NA_real_))
      }
    }
    null_parameters <- insert_interaction_zeros(refitted_null$fit$par)
    structure_i <- full_structure
    structure_i$data.tmb$yobs <- as.numeric(simulated[, 1])
    structure_i$data.tmb$size <- as.numeric(rowSums(simulated))
    objective_i <- glmmTMB:::fitTMB(structure_i, doOptim = FALSE)
    score <- score_from_objective(objective_i, null_parameters)
    c(
      statistic = score$statistic,
      nuisance = score$max_abs_nuisance_score,
      eigen = score$min_eigenvalue,
      condition = score$condition_number
    )
  }, error = function(error) {
    c(statistic = NA_real_, nuisance = NA_real_, eigen = NA_real_, condition = NA_real_)
  })
  result
}

message(sprintf(
  "Running %d refitted-null parametric score simulations on %d cores...",
  n_sim,
  n_cores
))
started <- proc.time()[[3]]
bootstrap_matrix <- do.call(
  rbind,
  mclapply(seq_len(n_sim), score_one, mc.cores = n_cores, mc.preschedule = TRUE)
)
elapsed <- proc.time()[[3]] - started

valid <- is.finite(bootstrap_matrix[, "statistic"])
valid_statistics <- bootstrap_matrix[valid, "statistic"]
exceedances <- sum(valid_statistics >= observed$statistic)
p_bootstrap <- (exceedances + 1) / (length(valid_statistics) + 1)
mcse <- sqrt(p_bootstrap * (1 - p_bootstrap) / (length(valid_statistics) + 1))

summary_payload <- list(
  test = "Efficient score test of the year-by-group interaction",
  null_model = "Additive beta-binomial GLMM with block-within-year and program random intercepts",
  calibration = "Parametric bootstrap with the additive null refitted for every simulated dataset",
  observed_score_statistic = observed$statistic,
  df = n_interaction,
  asymptotic_p = pchisq(observed$statistic, n_interaction, lower.tail = FALSE),
  requested_simulations = n_sim,
  valid_simulations = length(valid_statistics),
  failed_simulations = n_sim - length(valid_statistics),
  exceedances = exceedances,
  bootstrap_p = p_bootstrap,
  monte_carlo_se = mcse,
  seed = seed,
  cores = n_cores,
  elapsed_seconds = elapsed,
  observed_information_min_eigenvalue = observed$min_eigenvalue,
  observed_information_condition_number = observed$condition_number,
  max_abs_nuisance_score_observed = observed$max_abs_nuisance_score,
  bootstrap_quantiles = as.list(quantile(
    valid_statistics,
    c(0.5, 0.9, 0.95, 0.99)
  ))
)

write_json(
  summary_payload,
  file.path(out_dir, paste0(prefix, "_interaction_refit_score_bootstrap.json")),
  pretty = TRUE,
  auto_unbox = TRUE,
  digits = 15
)
write.csv(
  data.frame(
    Simulation = seq_len(n_sim),
    Score_Statistic = bootstrap_matrix[, "statistic"],
    Max_Abs_Nuisance_Score = bootstrap_matrix[, "nuisance"],
    Information_Min_Eigenvalue = bootstrap_matrix[, "eigen"],
    Information_Condition_Number = bootstrap_matrix[, "condition"]
  ),
  file.path(out_dir, paste0(prefix, "_interaction_refit_score_bootstrap_statistics.csv")),
  row.names = FALSE
)
cat(toJSON(summary_payload, pretty = TRUE, auto_unbox = TRUE, digits = 15), "\n")
