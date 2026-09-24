#!/usr/bin/env Rscript

# Focused sensitivity to excluding every record assigned to the two unresolved
# 2022 program identities. Uses the retained primary families, hierarchy,
# transformation, residual-dispersion structure, EMMs, and contrast procedures.
# It does NOT test an interaction or calibrate a severity omnibus test.
# Original files are read-only inputs. All outputs go to the explicit third
# argument, which must differ from the original data/model directories.

suppressPackageStartupMessages({
  library(glmmTMB)
  library(emmeans)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: Rscript run_additive_sensitivity.R DATA_DIR ORIGINAL_PRIMARY_DIR NEW_OUTPUT_DIR")
}
data_dir <- normalizePath(args[[1]], mustWork = TRUE)
primary_dir <- normalizePath(args[[2]], mustWork = TRUE)
dir.create(args[[3]], recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(args[[3]], mustWork = TRUE)
stopifnot(out_dir != data_dir, out_dir != primary_dir)
if (file.exists(file.path(out_dir, "sensitivity_primary_additive_models.rds"))) {
  stop("An existing sensitivity model file is present; choose a new output directory.")
}

options(contrasts = c("contr.treatment", "contr.poly"))
original_data <- read.csv(file.path(data_dir, "clean_plot_data_public.csv"),
                          stringsAsFactors = FALSE, check.names = FALSE)
excluded_ids <- c("Y2022_P02", "Y2022_P03")
excluded <- original_data$Public_Program_ID %in% excluded_ids
stopifnot(nrow(original_data) == 264L, sum(excluded) == 8L)
stopifnot(all(original_data$Year[excluded] == 2022L))
stopifnot(all(original_data$Program_Class[excluded] == "Conventional only"))
d <- original_data[!excluded, , drop = FALSE]
d$Year <- factor(d$Year, levels = 2021:2025)
d$Program_Class <- factor(d$Program_Class,
  levels = c("Untreated", "Biological only", "Mixed", "Conventional only"))
d$BlockYear <- factor(d$BlockYear)
d$Public_Program_ID <- factor(d$Public_Program_ID)
d$Cluster_Negative_Count <- d$N_Clusters_Rated - d$Cluster_Positive_Count
stopifnot(nrow(d) == 256L, nlevels(d$Public_Program_ID) == 64L)
stopifnot(all(d$Cluster_Negative_Count >= 0L))
stopifnot(max(abs(d$HA_Logit_Severity -
  log((d$Cluster_Severity_Pct + 0.5) / (100.5 - d$Cluster_Severity_Pct)))) < 1e-10)

write.csv(original_data[excluded, ], file.path(out_dir, "excluded_2022_rows.csv"), row.names = FALSE)
write.csv(d, file.path(out_dir, "analysis_data_256_plots.csv"), row.names = FALSE)
write.csv(as.data.frame(with(d, table(Year, Program_Class))),
          file.path(out_dir, "plot_support.csv"), row.names = FALSE)
program_data <- unique(d[, c("Year", "Public_Program_ID", "Program_Class")])
write.csv(as.data.frame(with(program_data, table(Year, Program_Class))),
          file.path(out_dir, "program_support.csv"), row.names = FALSE)

ctrl <- glmmTMBControl(optCtrl = list(iter.max = 20000L, eval.max = 20000L))
ctrl_bfgs <- glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"),
                           optCtrl = list(maxit = 20000L))
usable <- function(x) {
  !inherits(x, "try-error") && isTRUE(x$sdr$pdHess) &&
    identical(unname(x$fit$convergence), 0L) && is.finite(as.numeric(logLik(x)))
}
fit_model <- function(formula, family, REML = FALSE, dispformula = ~1) {
  x <- try(glmmTMB(formula, data = d, family = family, REML = REML,
                   dispformula = dispformula, control = ctrl), silent = TRUE)
  optimizer <- "nlminb"
  if (!usable(x)) {
    x <- try(glmmTMB(formula, data = d, family = family, REML = REML,
                     dispformula = dispformula, control = ctrl_bfgs), silent = TRUE)
    optimizer <- "BFGS"
  }
  if (!usable(x)) stop("Both optimizers failed for ", deparse(formula))
  list(model = x, optimizer = optimizer)
}

started <- proc.time()[[3]]
inc_result <- fit_model(
  cbind(Cluster_Positive_Count, Cluster_Negative_Count) ~
    Year + Program_Class + (1 | BlockYear) + (1 | Public_Program_ID),
  betabinomial(link = "logit"))
sev_result <- fit_model(
  HA_Logit_Severity ~ Year + Program_Class +
    (1 | BlockYear) + (1 | Public_Program_ID),
  gaussian(), REML = TRUE, dispformula = ~Year + Program_Class)
inc <- inc_result$model
sev <- sev_result$model

inv_ha <- function(z) pmin(100, pmax(0, 101 * plogis(z) - 0.5))
inc_emm <- emmeans(inc, ~Program_Class, weights = "equal", type = "response")
inc_est <- as.data.frame(summary(inc_emm, infer = c(TRUE, TRUE)))
inc_pairs <- as.data.frame(summary(pairs(inc_emm), infer = c(TRUE, TRUE), adjust = "holm"))
sev_emm <- emmeans(sev, ~Program_Class, weights = "equal", ddf = "kenward-roger", data = d)
sev_est <- as.data.frame(summary(sev_emm, infer = c(TRUE, TRUE)))
lower_name <- intersect(c("lower.CL", "asymp.LCL"), names(sev_est))[[1]]
upper_name <- intersect(c("upper.CL", "asymp.UCL"), names(sev_est))[[1]]
sev_est$Estimate_Pct <- inv_ha(sev_est$emmean)
sev_est$Lower_Pct <- inv_ha(sev_est[[lower_name]])
sev_est$Upper_Pct <- inv_ha(sev_est[[upper_name]])
sev_pairs <- as.data.frame(summary(pairs(sev_emm), infer = c(TRUE, TRUE), adjust = "holm"))
write.csv(inc_est, file.path(out_dir, "incidence_pooled_class_estimates.csv"), row.names = FALSE)
write.csv(inc_pairs, file.path(out_dir, "incidence_pooled_class_contrasts.csv"), row.names = FALSE)
write.csv(sev_est, file.path(out_dir, "severity_pooled_class_estimates_KR_approximate.csv"), row.names = FALSE)
write.csv(sev_pairs, file.path(out_dir, "severity_pooled_class_contrasts_KR_approximate.csv"), row.names = FALSE)

original <- readRDS(file.path(primary_dir, "primary_models.rds"))
parameters <- function(model, endpoint, dataset) {
  parts <- fixef(model)
  rows <- lapply(c("cond", "disp"), function(component) {
    v <- parts[[component]]
    if (!length(v)) return(NULL)
    data.frame(Endpoint = endpoint, Dataset = dataset, Component = component,
               Parameter = names(v), Value = as.numeric(v))
  })
  rv <- VarCorr(model)$cond
  for (group in names(rv)) {
    rows[[length(rows) + 1L]] <- data.frame(
      Endpoint = endpoint, Dataset = dataset, Component = "random_SD",
      Parameter = group, Value = unname(attr(rv[[group]], "stddev")[[1]]))
  }
  if (endpoint == "Incidence") {
    rows[[length(rows) + 1L]] <- data.frame(Endpoint = endpoint, Dataset = dataset,
      Component = "beta_binomial", Parameter = "rho", Value = 1 / (1 + sigma(model)))
  }
  do.call(rbind, rows)
}
parameter_table <- rbind(
  parameters(original$incidence$additive, "Incidence", "Original264"),
  parameters(inc, "Incidence", "Exclude8"),
  parameters(original$severity$additive_reml, "Severity", "Original264"),
  parameters(sev, "Severity", "Exclude8"))
write.csv(parameter_table, file.path(out_dir, "parameter_comparison_long.csv"), row.names = FALSE)
parameter_wide <- reshape(parameter_table, idvar = c("Endpoint", "Component", "Parameter"),
                          timevar = "Dataset", direction = "wide")
parameter_wide$Change <- parameter_wide$Value.Exclude8 - parameter_wide$Value.Original264
write.csv(parameter_wide, file.path(out_dir, "parameter_comparison.csv"), row.names = FALSE)

compare_table <- function(old_name, new_table, key, output_name) {
  old <- read.csv(file.path(primary_dir, old_name), stringsAsFactors = FALSE, check.names = FALSE)
  joined <- merge(old, new_table, by = key, suffixes = c("_Original264", "_Exclude8"), all = TRUE)
  write.csv(joined, file.path(out_dir, output_name), row.names = FALSE)
  joined
}
inc_comp <- compare_table("incidence_pooled_class_estimates.csv", inc_est, "Program_Class", "incidence_EMM_comparison.csv")
sev_comp <- compare_table("severity_pooled_class_estimates_kenward_roger.csv", sev_est, "Program_Class", "severity_EMM_comparison.csv")
inc_pair_comp <- compare_table("incidence_pooled_class_contrasts.csv", inc_pairs, "contrast", "incidence_pairwise_comparison.csv")
sev_pair_comp <- compare_table("severity_pooled_class_contrasts_kenward_roger.csv", sev_pairs, "contrast", "severity_pairwise_comparison.csv")

fit_info <- function(model, optimizer) list(
  optimizer = optimizer, convergence = unname(model$fit$convergence),
  pd_hessian = isTRUE(model$sdr$pdHess), logLik = as.numeric(logLik(model)),
  fixed_effects = as.list(fixef(model)$cond))
summary <- list(
  status = "Additive sensitivity completed; no interaction tests or new omnibus bootstrap inference performed",
  excluded_program_ids = excluded_ids, excluded_plots = 8L, retained_plots = nrow(d),
  retained_programs = nlevels(d$Public_Program_ID),
  incidence = fit_info(inc, inc_result$optimizer), severity = fit_info(sev, sev_result$optimizer),
  max_absolute_incidence_EMM_change_percentage_points =
    max(abs(100 * (inc_comp$prob_Exclude8 - inc_comp$prob_Original264))),
  max_absolute_severity_EMM_change_percentage_points =
    max(abs(sev_comp$Estimate_Pct_Exclude8 - sev_comp$Estimate_Pct_Original264)),
  incidence_pairwise_significance_changes =
    sum((inc_pair_comp$p.value_Original264 < 0.05) != (inc_pair_comp$p.value_Exclude8 < 0.05)),
  severity_pairwise_significance_changes =
    sum((sev_pair_comp$p.value_Original264 < 0.05) != (sev_pair_comp$p.value_Exclude8 < 0.05)),
  warnings = c("Severity KR confidence intervals and pairwise P values are approximate for modeled heterogeneous dispersion.",
               "Stable additive estimates would not establish stability of the interaction test or equivalence among classes.",
               "Exclusion leaves only one conventional program in 2022, as well as in 2024 and 2025."),
  elapsed_seconds = proc.time()[[3]] - started,
  versions = list(R = R.version.string, glmmTMB = as.character(packageVersion("glmmTMB")),
                  emmeans = as.character(packageVersion("emmeans"))))
write_json(summary, file.path(out_dir, "additive_sensitivity_summary.json"), pretty = TRUE,
           auto_unbox = TRUE, digits = 15)
saveRDS(list(data = d, incidence = list(additive = inc),
             severity = list(additive_reml = sev)),
        file.path(out_dir, "sensitivity_primary_additive_models.rds"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "R_session_info.txt"))
message("Focused additive sensitivity completed in ", round(summary$elapsed_seconds, 1), " seconds.")
