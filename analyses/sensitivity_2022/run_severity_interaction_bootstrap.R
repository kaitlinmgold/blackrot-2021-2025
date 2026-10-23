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
## Severity: refitted parametric-bootstrap LRTs under the selected variance
## model. The class and interaction tests use their respective null models.
## -------------------------------------------------------------------------

message("Preparing refitted severity likelihood-ratio bootstraps...")

sev_formula_year <- HA_Logit_Severity ~
  Year + (1 | BlockYear) + (1 | Public_Program_ID)
sev_formula_add <- update(sev_formula_year, . ~ . + Program_Class)
sev_formula_full <- update(sev_formula_add, . ~ . + Year:Program_Class)
sev_dispersion <- ~Year + Program_Class

severity_lrt_bootstrap <- function(null_fit, null_formula, alternative_formula,
                                   label, seed_used, stem) {
  observed_alternative <- switch(
    stem,
    severity_class_lrt_refit_bootstrap = sev_saved$additive_ml,
    severity_interaction_lrt_refit_bootstrap = sev_saved$full_ml,
    stop("unrecognized severity test stem")
  )
  observed_statistic <- max(
    0,
    2 * (
      as.numeric(logLik(observed_alternative)) -
        as.numeric(logLik(null_fit))
    )
  )

  set.seed(seed_used)
  simulations <- simulate(null_fit, nsim = n_sim, seed = seed_used)
  one <- function(index) {
    data_i <- d
    data_i$HA_Logit_Severity <- simulations[[index]]
    null_result <- fit_with_retry(
      null_formula,
      data_i,
      family = gaussian(),
      dispformula = sev_dispersion,
      REML = FALSE
    )
    alternative_result <- fit_with_retry(
      alternative_formula,
      data_i,
      family = gaussian(),
      dispformula = sev_dispersion,
      REML = FALSE
    )
    if (!identical(null_result$status, "ok") ||
        !identical(alternative_result$status, "ok")) {
      return(data.frame(
        Simulation = index,
        Statistic = NA_real_,
        Null_Optimizer = null_result$optimizer,
        Alternative_Optimizer = alternative_result$optimizer,
        Status = "fit_failure"
      ))
    }
    statistic <- max(
      0,
      2 * (
        as.numeric(logLik(alternative_result$fit)) -
          as.numeric(logLik(null_result$fit))
      )
    )
    data.frame(
      Simulation = index,
      Statistic = statistic,
      Null_Optimizer = null_result$optimizer,
      Alternative_Optimizer = alternative_result$optimizer,
      Status = "ok"
    )
  }

  started <- proc.time()[[3]]
  rows <- run_parallel(seq_len(n_sim), one)
  table_out <- do.call(rbind, rows)
  elapsed <- proc.time()[[3]] - started
  summary <- monte_carlo_summary(
    table_out$Statistic,
    observed_statistic,
    n_sim,
    label,
    paste(
      "Parametric-bootstrap likelihood-ratio test; both nested models were",
      "refitted by ML in every simulation while retaining the selected",
      "year-plus-program-class residual dispersion model"
    ),
    seed_used
  )
  summary$elapsed_seconds <- elapsed
  summary$status_counts <- as.list(table(table_out$Status))
  write_result(summary, table_out, stem)
  list(summary = summary, table = table_out)
}

sev_interaction_result <- severity_lrt_bootstrap(
  sev_saved$additive_ml,
  sev_formula_add,
  sev_formula_full,
  "Year-by-program-class interaction for transformed cluster severity",
  seed + 2L,
  "severity_interaction_lrt_refit_bootstrap"
)

overall_summary <- list(
  software=list(R=R.version.string,glmmTMB=as.character(packageVersion("glmmTMB")),emmeans=as.character(packageVersion("emmeans"))),
  simulations_per_test=n_sim, cores=n_cores, retained_plots=nrow(d),
  excluded_program_ids=c("Y2022_P02","Y2022_P03"),
  severity_interaction=sev_interaction_result$summary,
  scope="Exact original severity interaction calibration repeated after excluding eight ambiguous 2022 plots; incidence calibration is archived separately."
)
write_json(overall_summary,file.path(out_dir,"interaction_sensitivity_summary.json"),pretty=TRUE,auto_unbox=TRUE,digits=15)
writeLines(capture.output(sessionInfo()),file.path(out_dir,"R_session_info.txt"))
cat(toJSON(overall_summary,pretty=TRUE,auto_unbox=TRUE,digits=15),"\n")
