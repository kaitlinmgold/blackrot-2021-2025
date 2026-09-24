#!/usr/bin/env Rscript
# Fit secondary metal-status and component models from the study CSVs.
# Writes summaries and model objects to the second argument.
# Program attributes are not randomized factorial effects.

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

ctrl <- glmmTMBControl(optCtrl = list(iter.max = 20000L, eval.max = 20000L))

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

add_pct_columns <- function(frame, estimate_col, lower_candidates, upper_candidates) {
  lower_col <- intersect(lower_candidates, names(frame))[[1]]
  upper_col <- intersect(upper_candidates, names(frame))[[1]]
  frame$Estimate_Pct <- inverse_ha(frame[[estimate_col]])
  frame$Lower_Pct <- inverse_ha(frame[[lower_col]])
  frame$Upper_Pct <- inverse_ha(frame[[upper_col]])
  frame
}

score_interaction <- function(additive, full) {
  n_add_beta <- length(fixef(additive)$cond)
  n_full_beta <- length(fixef(full)$cond)
  n_interaction <- n_full_beta - n_add_beta
  add_par <- additive$fit$par
  null_par <- c(
    add_par[seq_len(n_add_beta)],
    rep(0, n_interaction),
    add_par[(n_add_beta + 1):length(add_par)]
  )
  gradient <- full$obj$gr(null_par)
  hessian <- optimHess(null_par, full$obj$fn, full$obj$gr)
  interaction_index <- (n_add_beta + 1):n_full_beta
  nuisance <- c(seq_len(n_add_beta), (n_full_beta + 1):length(null_par))
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
    df = n_interaction,
    p = pchisq(statistic, n_interaction, lower.tail = FALSE),
    max_abs_nuisance_score = max(abs(gradient[nuisance])),
    min_eigenvalue = min(eigen(efficient_information, symmetric = TRUE, only.values = TRUE)$values),
    condition_number = kappa(efficient_information)
  )
}

run_metal_analysis <- function(data, prefix, group_levels) {
  data$Analysis_Group <- factor(data$Analysis_Group, levels = group_levels)
  data$Year <- droplevels(data$Year)
  data$BlockYear <- droplevels(data$BlockYear)
  data$Public_Program_ID <- droplevels(data$Public_Program_ID)

  inc_year <- glmmTMB(
    cbind(Cluster_Positive_Count, Cluster_Negative_Count) ~
      Year + (1 | BlockYear) + (1 | Public_Program_ID),
    data = data, family = betabinomial(link = "logit"), control = ctrl
  )
  inc_add <- glmmTMB(
    cbind(Cluster_Positive_Count, Cluster_Negative_Count) ~
      Year + Analysis_Group + (1 | BlockYear) + (1 | Public_Program_ID),
    data = data, family = betabinomial(link = "logit"), control = ctrl
  )
  inc_full <- suppressWarnings(glmmTMB(
    cbind(Cluster_Positive_Count, Cluster_Negative_Count) ~
      Year * Analysis_Group + (1 | BlockYear) + (1 | Public_Program_ID),
    data = data, family = betabinomial(link = "logit"), control = ctrl
  ))
  inc_score <- score_interaction(inc_add, inc_full)
  inc_tests <- rbind(
    lrt_row("Analysis group", inc_year, inc_add),
    data.frame(
      Test = "Year x analysis group (efficient score)",
      Chisq = inc_score$statistic,
      df = inc_score$df,
      P = inc_score$p
    )
  )
  write.csv(inc_tests, file.path(out_dir, paste0(prefix, "_incidence_model_tests.csv")), row.names = FALSE)

  inc_emm <- emmeans(inc_add, ~Analysis_Group, weights = "equal", type = "response")
  inc_est <- as_plain_df(summary(inc_emm, infer = c(TRUE, TRUE)))
  inc_pairs <- as_plain_df(summary(pairs(inc_emm), infer = c(TRUE, TRUE), adjust = "holm"))
  inc_year_est <- as_plain_df(summary(
    emmeans(inc_add, ~Analysis_Group | Year, type = "response"),
    infer = c(TRUE, TRUE)
  ))
  write.csv(inc_est, file.path(out_dir, paste0(prefix, "_incidence_pooled_estimates.csv")), row.names = FALSE)
  write.csv(inc_pairs, file.path(out_dir, paste0(prefix, "_incidence_pooled_contrasts.csv")), row.names = FALSE)
  write.csv(inc_year_est, file.path(out_dir, paste0(prefix, "_incidence_year_estimates.csv")), row.names = FALSE)

  sev_formula_full <- HA_Logit_Severity ~ Year * Analysis_Group +
    (1 | BlockYear) + (1 | Public_Program_ID)
  sev_hom_reml <- glmmTMB(
    sev_formula_full, dispformula = ~1, data = data, family = gaussian(),
    REML = TRUE, control = ctrl
  )
  sev_year_reml <- glmmTMB(
    sev_formula_full, dispformula = ~Year, data = data, family = gaussian(),
    REML = TRUE, control = ctrl
  )
  sev_year_group_reml <- glmmTMB(
    sev_formula_full, dispformula = ~Year + Analysis_Group, data = data,
    family = gaussian(), REML = TRUE, control = ctrl
  )
  variance_comparison <- data.frame(
    Structure = c("Homogeneous", "Year", "Year plus analysis group"),
    logLik = c(logLik(sev_hom_reml), logLik(sev_year_reml), logLik(sev_year_group_reml)),
    npar = c(
      attr(logLik(sev_hom_reml), "df"),
      attr(logLik(sev_year_reml), "df"),
      attr(logLik(sev_year_group_reml), "df")
    ),
    AIC = c(AIC(sev_hom_reml), AIC(sev_year_reml), AIC(sev_year_group_reml)),
    BIC = c(BIC(sev_hom_reml), BIC(sev_year_reml), BIC(sev_year_group_reml))
  )
  write.csv(
    variance_comparison,
    file.path(out_dir, paste0(prefix, "_severity_variance_comparison.csv")),
    row.names = FALSE
  )
  variance_tests <- rbind(
    lrt_row("Homogeneous vs year-specific dispersion", sev_hom_reml, sev_year_reml),
    lrt_row(
      "Year-specific vs year-plus-analysis-group dispersion",
      sev_year_reml,
      sev_year_group_reml
    )
  )
  write.csv(
    variance_tests,
    file.path(out_dir, paste0(prefix, "_severity_variance_tests.csv")),
    row.names = FALSE
  )

  sev_year_ml <- glmmTMB(
    HA_Logit_Severity ~ Year + (1 | BlockYear) + (1 | Public_Program_ID),
    dispformula = ~Year, data = data, family = gaussian(),
    REML = FALSE, control = ctrl
  )
  sev_add_ml <- glmmTMB(
    HA_Logit_Severity ~ Year + Analysis_Group +
      (1 | BlockYear) + (1 | Public_Program_ID),
    dispformula = ~Year, data = data, family = gaussian(),
    REML = FALSE, control = ctrl
  )
  sev_full_ml <- glmmTMB(
    sev_formula_full, dispformula = ~Year, data = data,
    family = gaussian(), REML = FALSE, control = ctrl
  )
  sev_naive_lrt <- rbind(
    lrt_row("Analysis group", sev_year_ml, sev_add_ml),
    lrt_row("Year x analysis group", sev_add_ml, sev_full_ml)
  )
  write.csv(
    sev_naive_lrt,
    file.path(out_dir, paste0(prefix, "_severity_naive_ml_lrt_diagnostic.csv")),
    row.names = FALSE
  )

  sev_add_reml <- glmmTMB(
    HA_Logit_Severity ~ Year + Analysis_Group +
      (1 | BlockYear) + (1 | Public_Program_ID),
    dispformula = ~Year, data = data, family = gaussian(),
    REML = TRUE, control = ctrl
  )
  sev_pool_emm <- emmeans(
    sev_add_reml, ~Analysis_Group, weights = "equal",
    ddf = "kenward-roger", data = data
  )
  sev_pool_est <- add_pct_columns(
    as_plain_df(summary(sev_pool_emm, infer = c(TRUE, TRUE))),
    "emmean", c("lower.CL", "asymp.LCL"), c("upper.CL", "asymp.UCL")
  )
  sev_pool_pairs <- as_plain_df(summary(
    pairs(sev_pool_emm), infer = c(TRUE, TRUE), adjust = "holm"
  ))
  sev_year_emm <- emmeans(
    sev_add_reml, ~Analysis_Group | Year,
    ddf = "kenward-roger", data = data
  )
  sev_year_est <- add_pct_columns(
    as_plain_df(summary(sev_year_emm, infer = c(TRUE, TRUE))),
    "emmean", c("lower.CL", "asymp.LCL"), c("upper.CL", "asymp.UCL")
  )
  write.csv(sev_pool_est, file.path(out_dir, paste0(prefix, "_severity_pooled_estimates.csv")), row.names = FALSE)
  write.csv(sev_pool_pairs, file.path(out_dir, paste0(prefix, "_severity_pooled_contrasts.csv")), row.names = FALSE)
  write.csv(sev_year_est, file.path(out_dir, paste0(prefix, "_severity_year_estimates.csv")), row.names = FALSE)

  sev_full_reml <- glmmTMB(
    sev_formula_full, dispformula = ~Year, data = data, family = gaussian(),
    REML = TRUE, control = ctrl
  )
  sev_full_emm <- emmeans(
    sev_full_reml, ~Analysis_Group | Year,
    ddf = "kenward-roger", data = data
  )
  sev_joint_tests <- as_plain_df(joint_tests(
    sev_full_reml,
    ddf = "kenward-roger"
  ))
  write.csv(
    sev_joint_tests,
    file.path(out_dir, paste0(prefix, "_severity_kenward_roger_joint_tests.csv")),
    row.names = FALSE
  )
  sev_full_est <- add_pct_columns(
    as_plain_df(summary(sev_full_emm, infer = c(TRUE, TRUE))),
    "emmean", c("lower.CL", "asymp.LCL"), c("upper.CL", "asymp.UCL")
  )
  sev_full_pairs <- as_plain_df(summary(
    pairs(sev_full_emm), infer = c(TRUE, TRUE), adjust = "holm"
  ))
  write.csv(
    sev_full_est,
    file.path(out_dir, paste0(prefix, "_severity_interaction_year_estimates.csv")),
    row.names = FALSE
  )
  write.csv(
    sev_full_pairs,
    file.path(out_dir, paste0(prefix, "_severity_interaction_year_contrasts.csv")),
    row.names = FALSE
  )

  diagnostics <- list(
    n_plots = nrow(data),
    n_programs = length(unique(data$Public_Program_ID)),
    support = as.data.frame(with(data, table(Year, Analysis_Group))),
    incidence = list(
      year_pd_hessian = inc_year$sdr$pdHess,
      additive_pd_hessian = inc_add$sdr$pdHess,
      full_pd_hessian = inc_full$sdr$pdHess,
      full_convergence = inc_full$fit$convergence,
      full_message = inc_full$fit$message,
      phi = unname(sigma(inc_add)),
      rho = 1 / (1 + unname(sigma(inc_add))),
      score = inc_score,
      random_effects = capture.output(print(VarCorr(inc_add)))
    ),
    severity = list(
      selected_structure = "year-specific residual dispersion",
      selected_pd_hessian = sev_year_reml$sdr$pdHess,
      random_effects = capture.output(print(VarCorr(sev_year_reml)))
    )
  )
  write_json(
    diagnostics,
    file.path(out_dir, paste0(prefix, "_diagnostics.json")),
    pretty = TRUE, auto_unbox = TRUE, digits = 15
  )

  saveRDS(
    list(
      data = data,
      incidence = list(year = inc_year, additive = inc_add, full = inc_full),
      severity = list(
        homogeneous_reml = sev_hom_reml,
        year_reml = sev_year_reml,
        year_group_reml = sev_year_group_reml,
        additive_reml = sev_add_reml,
        full_reml = sev_full_reml,
        year_ml = sev_year_ml,
        additive_ml = sev_add_ml,
        full_ml = sev_full_ml
      )
    ),
    file.path(out_dir, paste0(prefix, "_models.rds"))
  )
}

message("Running any-metal analyses...")
metal <- subset(
  d,
  as.integer(as.character(Year)) %in% 2023:2025 &
    (Program_Class == "Untreated" | Program_Class == "Biological only")
)
metal$Analysis_Group <- ifelse(
  metal$Program_Class == "Untreated",
  "Untreated",
  metal$Biological_Metal_Group
)
run_metal_analysis(
  metal,
  "metal_any",
  c("Untreated", "Non-metal biological", "Metal-containing biological")
)

message("Running strict-metal analyses...")
strict <- subset(
  metal,
  Program_Class == "Untreated" |
    Biological_Metal_Group == "Non-metal biological" |
    Strict_Metal_Only == 1
)
strict$Analysis_Group <- ifelse(
  strict$Program_Class == "Untreated",
  "Untreated",
  ifelse(strict$Strict_Metal_Only == 1, "Metal-only biological", "Non-metal biological")
)
run_metal_analysis(
  strict,
  "metal_strict",
  c("Untreated", "Non-metal biological", "Metal-only biological")
)

message("Running exploratory mechanism-indicator analyses...")
bio <- droplevels(subset(d, Program_Class == "Biological only"))
indicator_names <- c(
  "Microbe_Present", "Biochemical_Present", "Metal_Present", "Unknown_Biological_Present"
)

inc_mech_base <- glmmTMB(
  cbind(Cluster_Positive_Count, Cluster_Negative_Count) ~
    Year + (1 | BlockYear) + (1 | Public_Program_ID),
  data = bio, family = betabinomial(link = "logit"), control = ctrl
)
inc_mech_full <- glmmTMB(
  cbind(Cluster_Positive_Count, Cluster_Negative_Count) ~
    Year + Microbe_Present + Biochemical_Present + Metal_Present +
    Unknown_Biological_Present + (1 | BlockYear) + (1 | Public_Program_ID),
  data = bio, family = betabinomial(link = "logit"), control = ctrl
)

sev_mech_base_ml <- glmmTMB(
  HA_Logit_Severity ~ Year + (1 | BlockYear) + (1 | Public_Program_ID),
  dispformula = ~Year, data = bio, family = gaussian(), REML = FALSE, control = ctrl
)
sev_mech_full_ml <- glmmTMB(
  HA_Logit_Severity ~ Year + Microbe_Present + Biochemical_Present + Metal_Present +
    Unknown_Biological_Present + (1 | BlockYear) + (1 | Public_Program_ID),
  dispformula = ~Year, data = bio, family = gaussian(), REML = FALSE, control = ctrl
)
sev_mech_full_reml <- update(sev_mech_full_ml, REML = TRUE)

mechanism_tests <- rbind(
  transform(lrt_row("Four nonexclusive indicators", inc_mech_base, inc_mech_full), Endpoint = "Incidence"),
  transform(lrt_row("Four nonexclusive indicators", sev_mech_base_ml, sev_mech_full_ml), Endpoint = "Severity")
)
mechanism_tests <- mechanism_tests[, c("Endpoint", "Test", "Chisq", "df", "P")]
write.csv(mechanism_tests, file.path(out_dir, "mechanism_global_tests_v4.csv"), row.names = FALSE)

inc_coef <- as.data.frame(summary(inc_mech_full)$coefficients$cond)
inc_coef$Term <- rownames(inc_coef)
inc_coef <- inc_coef[inc_coef$Term %in% indicator_names, ]
inc_coef$Endpoint <- "Incidence"
inc_coef$Holm_P <- p.adjust(inc_coef$`Pr(>|z|)`, method = "holm")

sev_coef <- as.data.frame(summary(sev_mech_full_reml, ddf = "kenward-roger")$coefficients$cond)
sev_coef$Term <- rownames(sev_coef)
sev_coef <- sev_coef[sev_coef$Term %in% indicator_names, ]
sev_coef$Endpoint <- "Severity"
sev_coef$Holm_P <- p.adjust(sev_coef$`Pr(>|t|)`, method = "holm")
names(sev_coef)[names(sev_coef) == "Pr(>|t|)"] <- "Raw_P"
names(inc_coef)[names(inc_coef) == "Pr(>|z|)"] <- "Raw_P"
mechanism_coefficients <- rbind(
  inc_coef[, c("Endpoint", "Term", "Estimate", "Std. Error", "z value", "Raw_P", "Holm_P")],
  transform(
    sev_coef[, c("Endpoint", "Term", "Estimate", "Std. Error", "t value", "Raw_P", "Holm_P")],
    `z value` = `t value`
  )[, c("Endpoint", "Term", "Estimate", "Std. Error", "z value", "Raw_P", "Holm_P")]
)
names(mechanism_coefficients)[names(mechanism_coefficients) == "z value"] <- "Test_Statistic"
write.csv(
  mechanism_coefficients,
  file.path(out_dir, "mechanism_coefficients_v4.csv"),
  row.names = FALSE
)

write_json(
  list(
    incidence_pd_hessian = inc_mech_full$sdr$pdHess,
    severity_pd_hessian = sev_mech_full_reml$sdr$pdHess,
    incidence_rho = 1 / (1 + unname(sigma(inc_mech_full))),
    incidence_random_effects = capture.output(print(VarCorr(inc_mech_full))),
    severity_random_effects = capture.output(print(VarCorr(sev_mech_full_reml)))
  ),
  file.path(out_dir, "mechanism_diagnostics_v4.json"),
  pretty = TRUE, auto_unbox = TRUE, digits = 15
)

message("Completed secondary R reanalyses.")
