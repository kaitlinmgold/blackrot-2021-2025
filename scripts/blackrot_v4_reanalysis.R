#!/usr/bin/env Rscript
# Fit primary class models from data/clean_plot_data_public.csv.
# Writes estimates, contrasts, diagnostics, and model objects to the second argument.
# Use the bootstrap script for calibrated omnibus tests.

suppressPackageStartupMessages({
  library(glmmTMB)
  library(emmeans)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[[1]] else stop("data directory required")
out_dir <- if (length(args) >= 2) args[[2]] else stop("output directory required")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

options(contrasts = c("contr.treatment", "contr.poly"))
set.seed(20260831)

plot_file <- file.path(data_dir, "clean_plot_data_public.csv")
d <- read.csv(plot_file, stringsAsFactors = FALSE, check.names = FALSE)

d$Year <- factor(d$Year, levels = 2021:2025)
d$Program_Class <- factor(
  d$Program_Class,
  levels = c("Untreated", "Biological only", "Mixed", "Conventional only")
)
d$BlockYear <- factor(d$BlockYear)
d$Public_Program_ID <- factor(d$Public_Program_ID)
d$Cluster_Negative_Count <- d$N_Clusters_Rated - d$Cluster_Positive_Count

stopifnot(nrow(d) == 264L)
stopifnot(all(d$Cluster_Negative_Count >= 0L))
stopifnot(all(!is.na(d$HA_Logit_Severity)))

ctrl <- glmmTMBControl(
  optCtrl = list(iter.max = 20000L, eval.max = 20000L)
)

fit_ok <- function(model) {
  list(
    pd_hessian = isTRUE(model$sdr$pdHess),
    convergence = unname(model$fit$convergence),
    message = model$fit$message %||% "",
    logLik = as.numeric(logLik(model)),
    AIC = as.numeric(AIC(model)),
    BIC = as.numeric(BIC(model)),
    npar = attr(logLik(model), "df")
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x

lrt_row <- function(label, small, large) {
  tab <- anova(small, large)
  data.frame(
    Test = label,
    Chisq = unname(tab$Chisq[2]),
    df = unname(tab$`Chi Df`[2]),
    P = unname(tab$`Pr(>Chisq)`[2]),
    stringsAsFactors = FALSE
  )
}

as_plain_df <- function(x) {
  out <- as.data.frame(x)
  rownames(out) <- NULL
  out
}

inverse_ha <- function(z) {
  ez <- exp(z)
  pmax(0, pmin(100, (100.5 * ez - 0.5) / (1 + ez)))
}

message("Fitting beta-binomial incidence models...")
inc_formula_year <- cbind(Cluster_Positive_Count, Cluster_Negative_Count) ~
  Year + (1 | BlockYear) + (1 | Public_Program_ID)
inc_formula_add <- update(inc_formula_year, . ~ . + Program_Class)
inc_formula_full <- update(inc_formula_add, . ~ . + Year:Program_Class)

inc_year <- glmmTMB(
  inc_formula_year, data = d, family = betabinomial(link = "logit"), control = ctrl
)
inc_add <- glmmTMB(
  inc_formula_add, data = d, family = betabinomial(link = "logit"), control = ctrl
)
inc_full <- suppressWarnings(glmmTMB(
  inc_formula_full, data = d, family = betabinomial(link = "logit"), control = ctrl
))

inc_tests <- rbind(
  lrt_row("Program class", inc_year, inc_add),
  lrt_row("Year x class (nonregular diagnostic LRT)", inc_add, inc_full)
)
write.csv(inc_tests, file.path(out_dir, "incidence_model_tests_glmmTMB.csv"), row.names = FALSE)

inc_emm <- emmeans(inc_add, ~ Program_Class, weights = "equal", type = "response")
inc_est <- as_plain_df(summary(inc_emm, infer = c(TRUE, TRUE)))
inc_pairs <- as_plain_df(summary(pairs(inc_emm), infer = c(TRUE, TRUE), adjust = "holm"))
write.csv(inc_est, file.path(out_dir, "incidence_pooled_class_estimates.csv"), row.names = FALSE)
write.csv(inc_pairs, file.path(out_dir, "incidence_pooled_class_contrasts.csv"), row.names = FALSE)

inc_year_emm <- emmeans(inc_add, ~ Program_Class | Year, type = "response")
inc_year_est <- as_plain_df(summary(inc_year_emm, infer = c(TRUE, TRUE)))
write.csv(inc_year_est, file.path(out_dir, "incidence_descriptive_year_estimates.csv"), row.names = FALSE)

inc_phi <- unname(sigma(inc_add))
inc_rho <- 1 / (1 + inc_phi)
inc_diagnostics <- list(
  software = list(
    R = R.version.string,
    glmmTMB = as.character(packageVersion("glmmTMB")),
    emmeans = as.character(packageVersion("emmeans"))
  ),
  year_only = fit_ok(inc_year),
  additive = fit_ok(inc_add),
  full_interaction = fit_ok(inc_full),
  beta_binomial_phi = inc_phi,
  intraplot_rho = inc_rho,
  additive_random_effects = capture.output(print(VarCorr(inc_add)))
)
write_json(
  inc_diagnostics,
  file.path(out_dir, "incidence_fit_diagnostics.json"),
  pretty = TRUE,
  auto_unbox = TRUE,
  digits = 15
)

message("Fitting heterogeneous-variance Gaussian severity models...")
sev_formula_full <- HA_Logit_Severity ~ Year * Program_Class +
  (1 | BlockYear) + (1 | Public_Program_ID)

sev_hom_reml <- glmmTMB(
  sev_formula_full, dispformula = ~1, data = d, family = gaussian(),
  REML = TRUE, control = ctrl
)
sev_year_reml <- glmmTMB(
  sev_formula_full, dispformula = ~Year, data = d, family = gaussian(),
  REML = TRUE, control = ctrl
)
sev_class_reml <- glmmTMB(
  sev_formula_full, dispformula = ~Program_Class, data = d, family = gaussian(),
  REML = TRUE, control = ctrl
)
sev_year_class_reml <- glmmTMB(
  sev_formula_full, dispformula = ~Year + Program_Class, data = d, family = gaussian(),
  REML = TRUE, control = ctrl
)

variance_models <- list(
  homogeneous = sev_hom_reml,
  year = sev_year_reml,
  class = sev_class_reml,
  year_plus_class = sev_year_class_reml
)
variance_comparison <- do.call(
  rbind,
  lapply(names(variance_models), function(name) {
    model <- variance_models[[name]]
    data.frame(
      Structure = name,
      logLik = unname(logLik(model)),
      npar = attr(logLik(model), "df"),
      AIC = AIC(model),
      BIC = BIC(model),
      pd_hessian = isTRUE(model$sdr$pdHess),
      convergence = unname(model$fit$convergence),
      stringsAsFactors = FALSE
    )
  })
)
write.csv(
  variance_comparison,
  file.path(out_dir, "severity_reml_variance_structure_comparison.csv"),
  row.names = FALSE
)

dispersion_lrt <- rbind(
  lrt_row("Homogeneous vs year-specific dispersion", sev_hom_reml, sev_year_reml),
  lrt_row("Year-specific vs year-plus-class dispersion", sev_year_reml, sev_year_class_reml),
  lrt_row("Class-specific vs year-plus-class dispersion", sev_class_reml, sev_year_class_reml)
)
write.csv(
  dispersion_lrt,
  file.path(out_dir, "severity_reml_dispersion_likelihood_ratio_tests.csv"),
  row.names = FALSE
)

sev_formula_year <- HA_Logit_Severity ~ Year +
  (1 | BlockYear) + (1 | Public_Program_ID)
sev_formula_add <- update(sev_formula_year, . ~ . + Program_Class)

sev_year_ml <- glmmTMB(
  sev_formula_year, dispformula = ~Year + Program_Class, data = d,
  family = gaussian(), REML = FALSE, control = ctrl
)
sev_add_ml <- glmmTMB(
  sev_formula_add, dispformula = ~Year + Program_Class, data = d,
  family = gaussian(), REML = FALSE, control = ctrl
)
sev_full_ml <- glmmTMB(
  sev_formula_full, dispformula = ~Year + Program_Class, data = d,
  family = gaussian(), REML = FALSE, control = ctrl
)

sev_naive_lrt <- rbind(
  lrt_row("Program class", sev_year_ml, sev_add_ml),
  lrt_row("Year x class", sev_add_ml, sev_full_ml)
)
write.csv(
  sev_naive_lrt,
  file.path(out_dir, "severity_naive_ml_likelihood_ratio_diagnostic.csv"),
  row.names = FALSE
)

sev_joint_tests <- as_plain_df(joint_tests(
  sev_year_class_reml,
  ddf = "kenward-roger"
))
write.csv(
  sev_joint_tests,
  file.path(out_dir, "severity_kenward_roger_joint_tests.csv"),
  row.names = FALSE
)

sev_add_reml <- glmmTMB(
  sev_formula_add,
  dispformula = ~Year + Program_Class,
  data = d,
  family = gaussian(),
  REML = TRUE,
  control = ctrl
)
sev_pooled_emm <- emmeans(
  sev_add_reml,
  ~ Program_Class,
  weights = "equal",
  ddf = "kenward-roger",
  data = d
)
sev_pooled_est <- as_plain_df(summary(sev_pooled_emm, infer = c(TRUE, TRUE)))
sev_pooled_lower <- intersect(c("lower.CL", "asymp.LCL"), names(sev_pooled_est))[[1]]
sev_pooled_upper <- intersect(c("upper.CL", "asymp.UCL"), names(sev_pooled_est))[[1]]
sev_pooled_est$Estimate_Pct <- inverse_ha(sev_pooled_est$emmean)
sev_pooled_est$Lower_Pct <- inverse_ha(sev_pooled_est[[sev_pooled_lower]])
sev_pooled_est$Upper_Pct <- inverse_ha(sev_pooled_est[[sev_pooled_upper]])
sev_pooled_pairs <- as_plain_df(summary(
  pairs(sev_pooled_emm),
  infer = c(TRUE, TRUE),
  adjust = "holm"
))
write.csv(
  sev_pooled_est,
  file.path(out_dir, "severity_pooled_class_estimates_kenward_roger.csv"),
  row.names = FALSE
)
write.csv(
  sev_pooled_pairs,
  file.path(out_dir, "severity_pooled_class_contrasts_kenward_roger.csv"),
  row.names = FALSE
)

sev_emm <- emmeans(
  sev_year_class_reml,
  ~ Program_Class | Year,
  ddf = "kenward-roger",
  data = d
)
sev_est_link <- as_plain_df(summary(sev_emm, infer = c(TRUE, TRUE)))
sev_lower_col <- intersect(c("lower.CL", "asymp.LCL"), names(sev_est_link))[[1]]
sev_upper_col <- intersect(c("upper.CL", "asymp.UCL"), names(sev_est_link))[[1]]
sev_est <- sev_est_link
sev_est$Estimate_Pct <- inverse_ha(sev_est$emmean)
sev_est$Lower_Pct <- inverse_ha(sev_est[[sev_lower_col]])
sev_est$Upper_Pct <- inverse_ha(sev_est[[sev_upper_col]])
sev_pairs <- as_plain_df(summary(pairs(sev_emm), infer = c(TRUE, TRUE), adjust = "holm"))
write.csv(sev_est, file.path(out_dir, "severity_year_estimates_kenward_roger.csv"), row.names = FALSE)
write.csv(sev_pairs, file.path(out_dir, "severity_year_contrasts_kenward_roger.csv"), row.names = FALSE)

pearson_resid <- residuals(sev_year_class_reml, type = "pearson")
conditional_resid <- residuals(sev_year_class_reml, type = "response")
resid_frame <- data.frame(
  Public_Plot_ID = d$Public_Plot_ID,
  Year = d$Year,
  Program_Class = d$Program_Class,
  Fitted_Link = fitted(sev_year_class_reml),
  Residual_Response = conditional_resid,
  Residual_Pearson = pearson_resid
)
write.csv(resid_frame, file.path(out_dir, "severity_conditional_residuals_glmmTMB.csv"), row.names = FALSE)

brown_forsythe <- function(residual, group, label) {
  group <- factor(group)
  centered <- abs(residual - ave(residual, group, FUN = median))
  fit <- lm(centered ~ group)
  tab <- anova(fit)
  data.frame(
    Residual_Model = label,
    Grouping = deparse(substitute(group)),
    F = unname(tab$`F value`[1]),
    df1 = unname(tab$Df[1]),
    df2 = unname(tab$Df[2]),
    P = unname(tab$`Pr(>F)`[1]),
    stringsAsFactors = FALSE
  )
}

bf_test <- function(residual, group, residual_model, grouping_label) {
  group <- factor(group)
  centered <- abs(residual - ave(residual, group, FUN = median))
  tab <- anova(lm(centered ~ group))
  data.frame(
    Residual_Model = residual_model,
    Grouping = grouping_label,
    F = unname(tab$`F value`[1]),
    df1 = unname(tab$Df[1]),
    df2 = unname(tab$Df[2]),
    P = unname(tab$`Pr(>F)`[1]),
    stringsAsFactors = FALSE
  )
}

homogeneous_pearson <- residuals(sev_hom_reml, type = "pearson")
residual_tests <- rbind(
  bf_test(homogeneous_pearson, d$Year, "Homogeneous", "Year"),
  bf_test(homogeneous_pearson, d$Program_Class, "Homogeneous", "Program class"),
  bf_test(
    homogeneous_pearson,
    interaction(d$Year, d$Program_Class, drop = TRUE),
    "Homogeneous",
    "Year x program class"
  ),
  bf_test(pearson_resid, d$Year, "Year plus class", "Year"),
  bf_test(pearson_resid, d$Program_Class, "Year plus class", "Program class"),
  bf_test(
    pearson_resid,
    interaction(d$Year, d$Program_Class, drop = TRUE),
    "Year plus class",
    "Year x program class"
  )
)
write.csv(
  residual_tests,
  file.path(out_dir, "severity_residual_variance_diagnostics.csv"),
  row.names = FALSE
)

sev_diagnostics <- list(
  homogeneous_reml = fit_ok(sev_hom_reml),
  year_reml = fit_ok(sev_year_reml),
  class_reml = fit_ok(sev_class_reml),
  year_plus_class_reml = fit_ok(sev_year_class_reml),
  year_ml = fit_ok(sev_year_ml),
  additive_ml = fit_ok(sev_add_ml),
  full_ml = fit_ok(sev_full_ml),
  additive_reml = fit_ok(sev_add_reml),
  selected_random_effects = capture.output(print(VarCorr(sev_year_class_reml))),
  selected_dispersion_coefficients = unname(fixef(sev_year_class_reml)$disp),
  shapiro_pearson = unclass(shapiro.test(pearson_resid))
)
write_json(
  sev_diagnostics,
  file.path(out_dir, "severity_fit_diagnostics.json"),
  pretty = TRUE,
  auto_unbox = TRUE,
  digits = 15
)

message("Checking alternate optimizer agreement...")
ctrl_bfgs <- glmmTMBControl(
  optimizer = optim,
  optArgs = list(method = "BFGS"),
  optCtrl = list(maxit = 20000L)
)
inc_add_bfgs <- glmmTMB(
  inc_formula_add,
  data = d,
  family = betabinomial(link = "logit"),
  control = ctrl_bfgs
)
sev_year_class_reml_bfgs <- glmmTMB(
  sev_formula_full,
  dispformula = ~Year + Program_Class,
  data = d,
  family = gaussian(),
  REML = TRUE,
  control = ctrl_bfgs
)
optimizer_agreement <- data.frame(
  Model = c("Incidence additive", "Severity full mean, year-plus-class dispersion"),
  Default_LogLik = c(as.numeric(logLik(inc_add)), as.numeric(logLik(sev_year_class_reml))),
  BFGS_LogLik = c(as.numeric(logLik(inc_add_bfgs)), as.numeric(logLik(sev_year_class_reml_bfgs))),
  Absolute_LogLik_Difference = c(
    abs(as.numeric(logLik(inc_add)) - as.numeric(logLik(inc_add_bfgs))),
    abs(as.numeric(logLik(sev_year_class_reml)) - as.numeric(logLik(sev_year_class_reml_bfgs)))
  ),
  Max_Absolute_Conditional_Beta_Difference = c(
    max(abs(fixef(inc_add)$cond - fixef(inc_add_bfgs)$cond)),
    max(abs(fixef(sev_year_class_reml)$cond - fixef(sev_year_class_reml_bfgs)$cond))
  ),
  Default_PD_Hessian = c(inc_add$sdr$pdHess, sev_year_class_reml$sdr$pdHess),
  BFGS_PD_Hessian = c(inc_add_bfgs$sdr$pdHess, sev_year_class_reml_bfgs$sdr$pdHess),
  stringsAsFactors = FALSE
)
write.csv(
  optimizer_agreement,
  file.path(out_dir, "primary_optimizer_agreement.csv"),
  row.names = FALSE
)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "R_session_info.txt"))
saveRDS(
  list(
    incidence = list(year = inc_year, additive = inc_add, full = inc_full),
    severity = list(
      homogeneous_reml = sev_hom_reml,
      year_reml = sev_year_reml,
      class_reml = sev_class_reml,
      year_plus_class_reml = sev_year_class_reml,
      additive_reml = sev_add_reml,
      year_ml = sev_year_ml,
      additive_ml = sev_add_ml,
      full_ml = sev_full_ml
    )
  ),
  file.path(out_dir, "primary_models.rds")
)

message("Completed standard-R primary-model fits.")
