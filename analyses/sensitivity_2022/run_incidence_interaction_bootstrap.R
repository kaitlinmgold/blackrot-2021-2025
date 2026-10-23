#!/usr/bin/env Rscript

## Calibrates the interaction test after excluding the eight unresolved
## 2022 records, leaving 256 plots. Refits the null model in each sample.

suppressPackageStartupMessages({
  library(glmmTMB)
  library(emmeans)
  library(jsonlite)
  library(parallel)
})

args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[[1]] else stop("data directory required")
model_file <- if (length(args) >= 2) args[[2]] else stop("primary model RDS required")
out_dir <- if (length(args) >= 3) args[[3]] else stop("output directory required")
n_sim <- if (length(args) >= 4) as.integer(args[[4]]) else 9999L
n_cores <- if (length(args) >= 5) as.integer(args[[5]]) else 4L
seed <- if (length(args) >= 6) as.integer(args[[6]]) else 20260831L
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(normalizePath(out_dir) != normalizePath(data_dir))
if(file.exists(file.path(out_dir,"interaction_sensitivity_summary.json"))) stop("Use a fresh output directory.")

options(contrasts = c("contr.treatment", "contr.poly"))

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
d$Cluster_Negative_Count <-
  d$N_Clusters_Rated - d$Cluster_Positive_Count

stopifnot(nrow(d) == 256L, nlevels(d$Public_Program_ID) == 64L)
stopifnot(all(d$Cluster_Negative_Count >= 0L))
stopifnot(n_sim >= 99L)

saved <- readRDS(model_file)
inc_saved <- saved$incidence
sev_saved <- saved$severity

ctrl_default <- glmmTMBControl(
  optCtrl = list(iter.max = 20000L, eval.max = 20000L)
)
ctrl_bfgs <- glmmTMBControl(
  optimizer = optim,
  optArgs = list(method = "BFGS"),
  optCtrl = list(maxit = 20000L)
)

is_usable_fit <- function(fit) {
  !inherits(fit, "try-error") &&
    isTRUE(fit$sdr$pdHess) &&
    identical(unname(fit$fit$convergence), 0L) &&
    is.finite(as.numeric(logLik(fit)))
}

fit_with_retry <- function(formula, data, family, dispformula = ~1,
                           REML = FALSE) {
  fit <- try(
    suppressWarnings(glmmTMB(
      formula,
      data = data,
      family = family,
      dispformula = dispformula,
      REML = REML,
      control = ctrl_default
    )),
    silent = TRUE
  )
  if (is_usable_fit(fit)) {
    return(list(fit = fit, optimizer = "nlminb", status = "ok"))
  }
  fit_bfgs <- try(
    suppressWarnings(glmmTMB(
      formula,
      data = data,
      family = family,
      dispformula = dispformula,
      REML = REML,
      control = ctrl_bfgs
    )),
    silent = TRUE
  )
  if (is_usable_fit(fit_bfgs)) {
    return(list(fit = fit_bfgs, optimizer = "BFGS", status = "ok"))
  }
  list(fit = NULL, optimizer = NA_character_, status = "fit_failure")
}

run_parallel <- function(indices, fun) {
  if (.Platform$OS.type == "windows" || n_cores <= 1L) {
    lapply(indices, fun)
  } else {
    mclapply(indices, fun, mc.cores = n_cores, mc.preschedule = TRUE)
  }
}

monte_carlo_summary <- function(statistics, observed, requested, label,
                                calibration, seed_used) {
  valid <- statistics[is.finite(statistics)]
  exceedances <- sum(valid >= observed)
  failed <- requested - length(valid)
  p_value <- (exceedances + 1) / (length(valid) + 1)
  list(
    test = label,
    calibration = calibration,
    observed_statistic = unname(observed),
    requested_simulations = requested,
    valid_simulations = length(valid),
    failed_simulations = failed,
    exceedances = exceedances,
    bootstrap_p = p_value,
    monte_carlo_se = sqrt(p_value * (1 - p_value) / (length(valid) + 1)),
    missing_result_p_lower = (exceedances + 1) / (requested + 1),
    missing_result_p_upper =
      (exceedances + failed + 1) / (requested + 1),
    seed = seed_used,
    quantiles = as.list(quantile(valid, c(0.5, 0.9, 0.95, 0.99)))
  )
}

write_result <- function(summary, table, stem) {
  write_json(
    summary,
    file.path(out_dir, paste0(stem, ".json")),
    pretty = TRUE,
    auto_unbox = TRUE,
    digits = 15
  )
  write.csv(
    table,
    file.path(out_dir, paste0(stem, "_statistics.csv")),
    row.names = FALSE
  )
}

## -------------------------------------------------------------------------
## Incidence: conventional refit parametric bootstrap of efficient score.
## -------------------------------------------------------------------------

message("Preparing refitted incidence interaction score bootstrap...")

inc_formula_add <- cbind(Cluster_Positive_Count, Cluster_Negative_Count) ~
  Year + Program_Class + (1 | BlockYear) + (1 | Public_Program_ID)
inc_formula_full <- update(inc_formula_add, . ~ . + Year:Program_Class)

inc_add <- inc_saved$additive
inc_full <- inc_saved$full

full_structure <- glmmTMB(
  inc_formula_full,
  family = betabinomial(link = "logit"),
  data = d,
  doFit = FALSE
)

n_additive_beta <- length(fixef(inc_add)$cond)
n_full_beta <- length(fixef(inc_full)$cond)
n_interaction_beta <- n_full_beta - n_additive_beta
interaction_index <- (n_additive_beta + 1L):n_full_beta

insert_interaction_zeros <- function(additive_parameters) {
  c(
    additive_parameters[seq_len(n_additive_beta)],
    rep(0, n_interaction_beta),
    additive_parameters[(n_additive_beta + 1L):length(additive_parameters)]
  )
}

score_components <- function(objective, null_parameters) {
  nuisance_index <- setdiff(seq_along(null_parameters), interaction_index)
  gradient <- objective$gr(null_parameters)
  hessian <- optimHess(null_parameters, objective$fn, objective$gr)
  hessian <- 0.5 * (hessian + t(hessian))

  information_nn <- hessian[nuisance_index, nuisance_index, drop = FALSE]
  information_in <- hessian[interaction_index, nuisance_index, drop = FALSE]
  information_ii <- hessian[interaction_index, interaction_index, drop = FALSE]

  efficient_information <- information_ii -
    information_in %*% solve(information_nn, t(information_in))
  efficient_information <-
    0.5 * (efficient_information + t(efficient_information))
  eigenvalues <- eigen(
    efficient_information,
    symmetric = TRUE,
    only.values = TRUE
  )$values
  if (any(!is.finite(eigenvalues)) || min(eigenvalues) <= 0) {
    stop("non-positive efficient information")
  }

  efficient_score <- gradient[interaction_index] -
    information_in %*% solve(
      information_nn,
      gradient[nuisance_index]
    )
  statistic <- as.numeric(
    t(efficient_score) %*%
      solve(efficient_information, efficient_score)
  )
  c(
    statistic = statistic,
    max_abs_nuisance_score = max(abs(gradient[nuisance_index])),
    information_min_eigenvalue = min(eigenvalues),
    information_condition_number = kappa(efficient_information)
  )
}

observed_null_parameters <- insert_interaction_zeros(inc_add$fit$par)
observed_score <- score_components(inc_full$obj, observed_null_parameters)

set.seed(seed)
inc_simulations <- simulate(inc_add, nsim = n_sim, seed = seed)

incidence_one <- function(index) {
  simulated <- inc_simulations[[index]]
  data_i <- d
  data_i$Cluster_Positive_Count <- simulated[, 1]
  data_i$Cluster_Negative_Count <- simulated[, 2]

  additive_result <- fit_with_retry(
    inc_formula_add,
    data_i,
    family = betabinomial(link = "logit")
  )
  if (!identical(additive_result$status, "ok")) {
    return(data.frame(
      Simulation = index,
      Statistic = NA_real_,
      Max_Abs_Nuisance_Score = NA_real_,
      Information_Min_Eigenvalue = NA_real_,
      Information_Condition_Number = NA_real_,
      Optimizer = NA_character_,
      Status = additive_result$status
    ))
  }

  additive_fit <- additive_result$fit
  null_parameters <- insert_interaction_zeros(additive_fit$fit$par)
  structure_i <- full_structure
  structure_i$data.tmb$yobs <- as.numeric(simulated[, 1])
  structure_i$data.tmb$size <- as.numeric(rowSums(simulated))
  objective <- glmmTMB:::fitTMB(structure_i, doOptim = FALSE)
  values <- try(score_components(objective, null_parameters), silent = TRUE)
  if (inherits(values, "try-error") || any(!is.finite(values))) {
    return(data.frame(
      Simulation = index,
      Statistic = NA_real_,
      Max_Abs_Nuisance_Score = NA_real_,
      Information_Min_Eigenvalue = NA_real_,
      Information_Condition_Number = NA_real_,
      Optimizer = additive_result$optimizer,
      Status = "score_failure"
    ))
  }
  data.frame(
    Simulation = index,
    Statistic = unname(values[["statistic"]]),
    Max_Abs_Nuisance_Score =
      unname(values[["max_abs_nuisance_score"]]),
    Information_Min_Eigenvalue =
      unname(values[["information_min_eigenvalue"]]),
    Information_Condition_Number =
      unname(values[["information_condition_number"]]),
    Optimizer = additive_result$optimizer,
    Status = "ok"
  )
}

inc_started <- proc.time()[[3]]
inc_rows <- run_parallel(seq_len(n_sim), incidence_one)
inc_table <- do.call(rbind, inc_rows)
inc_elapsed <- proc.time()[[3]] - inc_started
inc_summary <- monte_carlo_summary(
  inc_table$Statistic,
  observed_score[["statistic"]],
  n_sim,
  "Efficient score test of the year-by-program-class interaction",
  paste(
    "Conventional parametric bootstrap; the additive beta-binomial null",
    "was refitted and the nuisance projection and information matrix were",
    "recomputed in every simulation"
  ),
  seed
)
inc_summary$asymptotic_df <- n_interaction_beta
inc_summary$asymptotic_p <- pchisq(
  observed_score[["statistic"]],
  n_interaction_beta,
  lower.tail = FALSE
)
inc_summary$observed_max_abs_nuisance_score <-
  observed_score[["max_abs_nuisance_score"]]
inc_summary$observed_information_min_eigenvalue <-
  observed_score[["information_min_eigenvalue"]]
inc_summary$observed_information_condition_number <-
  observed_score[["information_condition_number"]]
inc_summary$elapsed_seconds <- inc_elapsed
inc_summary$status_counts <- as.list(table(inc_table$Status))
inc_summary$optimizer_counts <- list(
  nlminb = sum(inc_table$Optimizer == "nlminb", na.rm = TRUE),
  BFGS = sum(inc_table$Optimizer == "BFGS", na.rm = TRUE),
  missing = sum(is.na(inc_table$Optimizer))
)
write_result(
  inc_summary,
  inc_table,
  "incidence_interaction_score_refit_bootstrap"
)

overall_summary <- list(
  software=list(R=R.version.string,glmmTMB=as.character(packageVersion("glmmTMB")),emmeans=as.character(packageVersion("emmeans"))),
  simulations_per_test=n_sim, cores=n_cores, retained_plots=nrow(d),
  excluded_program_ids=c("Y2022_P02","Y2022_P03"),
  incidence_interaction=inc_summary,
  scope="Exact original incidence interaction calibration repeated after excluding eight ambiguous 2022 plots; severity calibration is archived separately."
)
write_json(overall_summary,file.path(out_dir,"interaction_sensitivity_summary.json"),pretty=TRUE,auto_unbox=TRUE,digits=15)
writeLines(capture.output(sessionInfo()),file.path(out_dir,"R_session_info.txt"))
cat(toJSON(overall_summary,pretty=TRUE,auto_unbox=TRUE,digits=15),"\n")
