#!/usr/bin/env Rscript

# Regenerate secondary, non-primary descriptive outputs for the black rot v4
# manuscript from the de-identified public CSV files. This script does not
# create figures and does not modify any source file.

suppressPackageStartupMessages(library(emmeans))

args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1L) args[[1L]] else file.path("data")
output_dir <- if (length(args) >= 2L) args[[2L]] else file.path("analysis", "results_descriptive_v4")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- c(
  "clean_plot_data_public.csv",
  "cluster_organ_data_public.csv",
  "retained_program_codebook_public.csv"
)
missing_files <- required_files[!file.exists(file.path(data_dir, required_files))]
if (length(missing_files) > 0L) {
  stop("Missing required public data file(s): ", paste(missing_files, collapse = ", "))
}

read_public_csv <- function(filename) {
  read.csv(
    file.path(data_dir, filename),
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA")
  )
}

plot_data <- read_public_csv("clean_plot_data_public.csv")
organ_data <- read_public_csv("cluster_organ_data_public.csv")
codebook <- read_public_csv("retained_program_codebook_public.csv")

years <- 2021:2025
class_order <- c("Untreated", "Biological only", "Mixed", "Conventional only")
architecture_order <- c(
  "Alternating rotation",
  "Complex sequence",
  "Formulated hybrid product",
  "Repeated tank mix"
)
metal_order <- c("Untreated", "Non-metal biological", "Metal-containing biological")
unresolved_flag <- "2022_T2_T3_block_assignment_unresolved"

inv_ha_logit <- function(x) {
  pmin(100, pmax(0, 101 * plogis(x) - 0.5))
}

mean_ci <- function(x) {
  x <- as.numeric(x[!is.na(x)])
  n <- length(x)
  if (n == 0L) {
    return(c(N = 0, Mean = NA_real_, Lower95 = NA_real_, Upper95 = NA_real_))
  }
  mean_x <- mean(x)
  if (n == 1L) {
    return(c(N = 1, Mean = mean_x, Lower95 = mean_x, Upper95 = mean_x))
  }
  half_width <- qt(0.975, df = n - 1L) * stats::sd(x) / sqrt(n)
  c(
    N = n,
    Mean = mean_x,
    Lower95 = max(0, mean_x - half_width),
    Upper95 = min(100, mean_x + half_width)
  )
}

collapse_flags <- function(x) {
  values <- sort(unique(x[!is.na(x) & nzchar(x)]))
  if (length(values) == 0L) NA_character_ else paste(values, collapse = "; ")
}

single_value <- function(x, label) {
  values <- unique(x[!is.na(x)])
  if (length(values) != 1L) {
    stop("Expected one ", label, " value per group, found ", length(values))
  }
  values[[1L]]
}

add_validation <- local({
  rows <- list()
  function(check = NULL, observed = NULL, expected = NULL, status = NULL, note = NULL, finish = FALSE) {
    if (!finish) {
      if (is.null(status)) {
        status <- if (identical(as.character(observed), as.character(expected))) "PASS" else "FAIL"
      }
      rows[[length(rows) + 1L]] <<- data.frame(
        Check = as.character(check),
        Observed = as.character(observed),
        Expected = as.character(expected),
        Status = as.character(status),
        Note = if (is.null(note)) "" else as.character(note),
        stringsAsFactors = FALSE
      )
      invisible(NULL)
    } else {
      do.call(rbind, rows)
    }
  }
})

# Basic source validations.
add_validation("R version", paste(R.version$major, R.version$minor, sep = "."), "4.5.3")
add_validation("Plot rows", nrow(plot_data), 264)
add_validation("Unique plots", length(unique(plot_data$Public_Plot_ID)), nrow(plot_data))
add_validation("Retained program rows", nrow(codebook), 66)
add_validation("Unique retained programs", length(unique(plot_data$Public_Program_ID)), nrow(codebook))
add_validation("Expected organ rows", nrow(organ_data), sum(plot_data$N_Clusters_Expected))
add_validation("Analyzed organ ratings", sum(!is.na(organ_data$Analysis_Rating_Pct)), sum(plot_data$N_Clusters_Rated))

required_plot_columns <- c(
  "Public_Plot_ID", "Public_Program_ID", "Year", "Block", "Program_Label",
  "Program_Class", "Biological_Metal_Group", "Architecture",
  "N_Clusters_Rated", "Cluster_Positive_Count", "Cluster_Incidence_Pct",
  "Cluster_Severity_Pct", "HA_Logit_Incidence", "HA_Logit_Severity", "QC_Flags"
)
if (!all(required_plot_columns %in% names(plot_data))) {
  stop("The plot data are missing required columns: ", paste(setdiff(required_plot_columns, names(plot_data)), collapse = ", "))
}

# Validate the organ-level public file against every plot-level cluster endpoint.
analyzed_organ_data <- organ_data[!is.na(organ_data$Analysis_Rating_Pct), , drop = FALSE]
organ_by_plot <- aggregate(
  analyzed_organ_data$Analysis_Rating_Pct,
  by = list(Public_Plot_ID = analyzed_organ_data$Public_Plot_ID),
  FUN = function(x) c(N = length(x), Positive = sum(x > 0), Severity = mean(x))
)
organ_check <- data.frame(
  Public_Plot_ID = organ_by_plot$Public_Plot_ID,
  N = organ_by_plot$x[, "N"],
  Positive = organ_by_plot$x[, "Positive"],
  Severity = organ_by_plot$x[, "Severity"],
  stringsAsFactors = FALSE
)
organ_check <- merge(
  organ_check,
  plot_data[, c("Public_Plot_ID", "N_Clusters_Rated", "Cluster_Positive_Count", "Cluster_Severity_Pct")],
  by = "Public_Plot_ID",
  all = TRUE,
  sort = FALSE
)
organ_endpoint_match <- nrow(organ_check) == nrow(plot_data) &&
  all(organ_check$N == organ_check$N_Clusters_Rated) &&
  all(organ_check$Positive == organ_check$Cluster_Positive_Count) &&
  max(abs(organ_check$Severity - organ_check$Cluster_Severity_Pct), na.rm = TRUE) < 1e-10
add_validation(
  "Organ-to-plot endpoint reconstruction",
  organ_endpoint_match,
  TRUE,
  note = "Rated count, positive count, and arithmetic mean severity agree for every retained plot."
)

# Validate the plot-to-codebook classification crosswalk.
codebook_xwalk <- codebook[, c(
  "Public_Program_ID", "Year", "Program_Label_Clean", "Program_Class_Clean",
  "Biological_Metal_Group", "Architecture_Clean"
)]
plot_xwalk <- unique(plot_data[, c(
  "Public_Program_ID", "Year", "Program_Label", "Program_Class",
  "Biological_Metal_Group", "Architecture"
)])
xwalk <- merge(plot_xwalk, codebook_xwalk, by = c("Public_Program_ID", "Year"), all = TRUE)
xwalk_match <- nrow(xwalk) == nrow(codebook) &&
  all(xwalk$Program_Label == xwalk$Program_Label_Clean) &&
  all(xwalk$Program_Class == xwalk$Program_Class_Clean) &&
  all(xwalk$Biological_Metal_Group.x == xwalk$Biological_Metal_Group.y) &&
  all(xwalk$Architecture == xwalk$Architecture_Clean)
add_validation("Plot-to-codebook crosswalk", xwalk_match, TRUE)

# Program-versus-untreated comparisons. These are intentionally separate,
# within-year fixed-effects models, retained only as secondary program-level
# context. The unresolved 2022 P02/P03 records are not used for these claims.
unresolved_rows <- !is.na(plot_data$QC_Flags) & grepl(unresolved_flag, plot_data$QC_Flags, fixed = TRUE)
unresolved_ids <- sort(unique(plot_data$Public_Program_ID[unresolved_rows]))
program_model_data <- plot_data[!plot_data$Public_Program_ID %in% unresolved_ids, , drop = FALSE]

add_validation("Unresolved 2022 plot rows excluded from program models", sum(unresolved_rows), 8)
add_validation(
  "Unresolved 2022 program IDs excluded from program models",
  paste(unresolved_ids, collapse = ";"),
  "Y2022_P02;Y2022_P03"
)
add_validation("Plot rows used in program models", nrow(program_model_data), 256)

endpoint_map <- c(Incidence = "HA_Logit_Incidence", Severity = "HA_Logit_Severity")
program_contrast_rows <- list()
program_model_support_rows <- list()

for (endpoint in names(endpoint_map)) {
  response <- unname(endpoint_map[[endpoint]])
  for (year in years) {
    year_data <- program_model_data[program_model_data$Year == year, , drop = FALSE]
    utc_ids <- unique(year_data$Public_Program_ID[year_data$Program_Class == "Untreated"])
    if (length(utc_ids) != 1L) stop("Expected exactly one untreated program in ", year)
    utc_id <- utc_ids[[1L]]
    program_levels <- c(utc_id, setdiff(sort(unique(year_data$Public_Program_ID)), utc_id))
    year_data$Program_Model <- factor(year_data$Public_Program_ID, levels = program_levels)
    year_data$Block_Model <- factor(year_data$Block, levels = sort(unique(year_data$Block)))

    model_formula <- reformulate(c("Program_Model", "Block_Model"), response = response)
    fit <- lm(model_formula, data = year_data)
    if (qr(model.matrix(fit))$rank != ncol(model.matrix(fit))) {
      stop("Rank-deficient program/block model for ", endpoint, " in ", year)
    }

    emm <- emmeans(fit, specs = "Program_Model", weights = "equal")
    emm_ci <- as.data.frame(confint(emm, level = 0.95, adjust = "none"))
    treatment_contrasts <- as.data.frame(
      summary(
        contrast(emm, method = "trt.vs.ctrl", ref = 1L, adjust = "none"),
        infer = c(TRUE, TRUE),
        adjust = "none"
      )
    )
    treatment_contrasts$p_Holm <- p.adjust(treatment_contrasts$p.value, method = "holm")

    utc_emm <- emm_ci[emm_ci$Program_Model == utc_id, , drop = FALSE]
    treatment_ids <- program_levels[-1L]
    treatment_emm <- emm_ci[match(treatment_ids, emm_ci$Program_Model), , drop = FALSE]
    treatment_tests <- treatment_contrasts[match(
      paste(treatment_ids, "-", utc_id),
      treatment_contrasts$contrast
    ), , drop = FALSE]

    metadata <- year_data[match(treatment_ids, year_data$Public_Program_ID), c(
      "Public_Program_ID", "Program_Label", "Program_Class"
    )]
    n_plots <- as.integer(table(factor(year_data$Public_Program_ID, levels = program_levels))[treatment_ids])

    direction <- ifelse(
      abs(treatment_tests$estimate) < 1e-12,
      "No numerical difference",
      ifelse(treatment_tests$estimate < 0, "Lower disease", "Higher disease")
    )

    program_contrast_rows[[length(program_contrast_rows) + 1L]] <- data.frame(
      Endpoint = endpoint,
      Year = year,
      Public_Program_ID = treatment_ids,
      Program_Label = metadata$Program_Label,
      Program_Class = metadata$Program_Class,
      N_Plots = n_plots,
      EMM_Pct = inv_ha_logit(treatment_emm$emmean),
      Lower95_Pct = inv_ha_logit(treatment_emm$lower.CL),
      Upper95_Pct = inv_ha_logit(treatment_emm$upper.CL),
      Difference_vs_UTC_Logit = treatment_tests$estimate,
      SE = treatment_tests$SE,
      t = treatment_tests$t.ratio,
      Residual_df = treatment_tests$df,
      p_raw = treatment_tests$p.value,
      p_Holm = treatment_tests$p_Holm,
      Direction_vs_UTC = direction,
      Untreated_EMM_Pct = rep(inv_ha_logit(utc_emm$emmean), length(treatment_ids)),
      Untreated_Lower95_Pct = rep(inv_ha_logit(utc_emm$lower.CL), length(treatment_ids)),
      Untreated_Upper95_Pct = rep(inv_ha_logit(utc_emm$upper.CL), length(treatment_ids)),
      stringsAsFactors = FALSE
    )

    program_model_support_rows[[length(program_model_support_rows) + 1L]] <- data.frame(
      Endpoint = endpoint,
      Year = year,
      N_Plots = nrow(year_data),
      N_Programs = length(program_levels),
      N_Treatment_Contrasts = length(treatment_ids),
      Residual_df = stats::df.residual(fit),
      Residual_SD_Logit = summary(fit)$sigma,
      Unresolved_2022_IDs_Excluded = if (year == 2022) paste(unresolved_ids, collapse = ";") else "",
      Model = paste(response, "~ program + block"),
      Multiplicity = "Holm within endpoint and year",
      stringsAsFactors = FALSE
    )
  }
}

program_contrasts <- do.call(rbind, program_contrast_rows)
program_model_support <- do.call(rbind, program_model_support_rows)
program_contrasts <- program_contrasts[
  order(match(program_contrasts$Endpoint, names(endpoint_map)), program_contrasts$Year, program_contrasts$Public_Program_ID),
]
program_model_support <- program_model_support[
  order(match(program_model_support$Endpoint, names(endpoint_map)), program_model_support$Year),
]

write.csv(program_contrasts, file.path(output_dir, "program_vs_utc_contrasts.csv"), row.names = FALSE, na = "")
write.csv(program_model_support, file.path(output_dir, "program_fixed_effect_model_support.csv"), row.names = FALSE, na = "")

expected_contrasts_by_year <- c(`2021` = 12L, `2022` = 13L, `2023` = 14L, `2024` = 9L, `2025` = 11L)
observed_contrasts <- table(program_contrasts$Endpoint, program_contrasts$Year)
contrast_counts_match <- all(observed_contrasts["Incidence", ] == expected_contrasts_by_year) &&
  all(observed_contrasts["Severity", ] == expected_contrasts_by_year)
add_validation("Program-versus-UTC contrast rows", nrow(program_contrasts), 118)
add_validation("Program contrast counts by endpoint/year", contrast_counts_match, TRUE)

# Unadjusted program means and t-based 95% confidence intervals for Figure S1.
program_keys <- unique(plot_data[, c("Year", "Public_Program_ID")])
program_keys <- program_keys[order(program_keys$Year, program_keys$Public_Program_ID), ]
program_summary_rows <- lapply(seq_len(nrow(program_keys)), function(i) {
  key <- program_keys[i, ]
  group <- plot_data[
    plot_data$Year == key$Year & plot_data$Public_Program_ID == key$Public_Program_ID,
    ,
    drop = FALSE
  ]
  inc <- mean_ci(group$Cluster_Incidence_Pct)
  sev <- mean_ci(group$Cluster_Severity_Pct)
  data.frame(
    Year = key$Year,
    Public_Program_ID = key$Public_Program_ID,
    Program_Label = single_value(group$Program_Label, "program label"),
    Program_Class = single_value(group$Program_Class, "program class"),
    Biological_Metal_Group = single_value(group$Biological_Metal_Group, "metal group"),
    Architecture = single_value(group$Architecture, "architecture"),
    N_Plots = nrow(group),
    Incidence_Mean_Pct = inc[["Mean"]],
    Incidence_Lower95_Pct = inc[["Lower95"]],
    Incidence_Upper95_Pct = inc[["Upper95"]],
    Severity_Mean_Pct = sev[["Mean"]],
    Severity_Lower95_Pct = sev[["Lower95"]],
    Severity_Upper95_Pct = sev[["Upper95"]],
    QC_Flags = collapse_flags(group$QC_Flags),
    stringsAsFactors = FALSE
  )
})
program_summary <- do.call(rbind, program_summary_rows)
write.csv(program_summary, file.path(output_dir, "program_summary.csv"), row.names = FALSE, na = "")
add_validation("Unadjusted program-summary rows", nrow(program_summary), 66)
add_validation("Program-summary plot-count total", sum(program_summary$N_Plots), nrow(plot_data))

# Mixed-program architecture summaries. These are unadjusted and descriptive
# because architecture categories are sparse and structurally imbalanced.
mixed_data <- plot_data[plot_data$Program_Class == "Mixed", , drop = FALSE]
unexpected_architectures <- setdiff(unique(mixed_data$Architecture), architecture_order)
if (length(unexpected_architectures) > 0L) {
  stop("Unexpected mixed-program architecture(s): ", paste(unexpected_architectures, collapse = ", "))
}
architecture_keys <- unique(mixed_data[, c("Year", "Architecture")])
architecture_keys$Architecture <- factor(architecture_keys$Architecture, levels = architecture_order)
architecture_keys <- architecture_keys[order(architecture_keys$Year, architecture_keys$Architecture), ]
architecture_keys$Architecture <- as.character(architecture_keys$Architecture)

architecture_summary_rows <- lapply(seq_len(nrow(architecture_keys)), function(i) {
  key <- architecture_keys[i, ]
  group <- mixed_data[
    mixed_data$Year == key$Year & mixed_data$Architecture == key$Architecture,
    ,
    drop = FALSE
  ]
  inc <- mean_ci(group$Cluster_Incidence_Pct)
  sev <- mean_ci(group$Cluster_Severity_Pct)
  data.frame(
    Year = key$Year,
    Architecture = key$Architecture,
    N_Programs = length(unique(group$Public_Program_ID)),
    N_Plots = nrow(group),
    N_Clusters_Rated = sum(group$N_Clusters_Rated),
    Incidence_Mean_Pct = inc[["Mean"]],
    Incidence_Lower95_Pct = inc[["Lower95"]],
    Incidence_Upper95_Pct = inc[["Upper95"]],
    Severity_Mean_Pct = sev[["Mean"]],
    Severity_Lower95_Pct = sev[["Lower95"]],
    Severity_Upper95_Pct = sev[["Upper95"]],
    Inference_Status = "Descriptive only; architecture categories are sparse and structurally imbalanced across years",
    stringsAsFactors = FALSE
  )
})
architecture_summary <- do.call(rbind, architecture_summary_rows)
write.csv(architecture_summary, file.path(output_dir, "mixed_architecture_descriptive.csv"), row.names = FALSE, na = "")
add_validation("Observed mixed architecture cells", nrow(architecture_summary), 10)
add_validation("Mixed programs represented", sum(architecture_summary$N_Programs), 15)
add_validation("Mixed plot rows represented", sum(architecture_summary$N_Plots), nrow(mixed_data))
add_validation(
  "Mixed architecture labels",
  paste(sort(unique(architecture_summary$Architecture)), collapse = ";"),
  paste(sort(architecture_order), collapse = ";")
)

# Support by year and program class.
class_keys <- expand.grid(
  Year = years,
  Program_Class = class_order,
  stringsAsFactors = FALSE
)
class_support_rows <- lapply(seq_len(nrow(class_keys)), function(i) {
  key <- class_keys[i, ]
  group <- plot_data[
    plot_data$Year == key$Year & plot_data$Program_Class == key$Program_Class,
    ,
    drop = FALSE
  ]
  data.frame(
    Year = key$Year,
    Program_Class = key$Program_Class,
    N_Programs = length(unique(group$Public_Program_ID)),
    N_Plots = nrow(group),
    N_Clusters_Rated = sum(group$N_Clusters_Rated),
    N_Cluster_Positive = sum(group$Cluster_Positive_Count),
    stringsAsFactors = FALSE
  )
})
class_support <- do.call(rbind, class_support_rows)
class_support$Program_Class <- factor(class_support$Program_Class, levels = class_order)
class_support <- class_support[order(class_support$Year, class_support$Program_Class), ]
class_support$Program_Class <- as.character(class_support$Program_Class)
write.csv(class_support, file.path(output_dir, "year_class_support.csv"), row.names = FALSE, na = "")
add_validation("Year-by-class support rows", nrow(class_support), 20)
add_validation("Year-by-class programs total", sum(class_support$N_Programs), 66)
add_validation("Year-by-class plots total", sum(class_support$N_Plots), 264)
add_validation("Year-by-class clusters total", sum(class_support$N_Clusters_Rated), 5277)

# Support for the 2023-2025 biological metal-status analysis, with untreated
# plots included as the reference group.
metal_data <- plot_data[
  plot_data$Year %in% 2023:2025 & plot_data$Program_Class %in% c("Untreated", "Biological only"),
  ,
  drop = FALSE
]
metal_data$Metal_Analysis_Group <- ifelse(
  metal_data$Program_Class == "Untreated",
  "Untreated",
  metal_data$Biological_Metal_Group
)
metal_keys <- expand.grid(
  Year = 2023:2025,
  Metal_Analysis_Group = metal_order,
  stringsAsFactors = FALSE
)
metal_support_rows <- lapply(seq_len(nrow(metal_keys)), function(i) {
  key <- metal_keys[i, ]
  group <- metal_data[
    metal_data$Year == key$Year & metal_data$Metal_Analysis_Group == key$Metal_Analysis_Group,
    ,
    drop = FALSE
  ]
  data.frame(
    Year = key$Year,
    Metal_Analysis_Group = key$Metal_Analysis_Group,
    N_Plots = nrow(group),
    N_Programs = length(unique(group$Public_Program_ID)),
    N_Clusters_Rated = sum(group$N_Clusters_Rated),
    N_Cluster_Positive = sum(group$Cluster_Positive_Count),
    stringsAsFactors = FALSE
  )
})
metal_support <- do.call(rbind, metal_support_rows)
metal_support$Metal_Analysis_Group <- factor(metal_support$Metal_Analysis_Group, levels = metal_order)
metal_support <- metal_support[order(metal_support$Year, metal_support$Metal_Analysis_Group), ]
metal_support$Metal_Analysis_Group <- as.character(metal_support$Metal_Analysis_Group)
write.csv(metal_support, file.path(output_dir, "metal_support.csv"), row.names = FALSE, na = "")
add_validation("Year-by-metal-group support rows", nrow(metal_support), 9)
add_validation("Metal-analysis programs across year-specific cells", sum(metal_support$N_Programs), 23)
add_validation("Metal-analysis plot rows", sum(metal_support$N_Plots), nrow(metal_data))

validation_summary <- add_validation(finish = TRUE)
if (any(validation_summary$Status != "PASS")) {
  write.csv(validation_summary, file.path(output_dir, "validation_summary.csv"), row.names = FALSE, na = "")
  print(validation_summary, row.names = FALSE)
  stop("One or more validation checks failed; see validation_summary.csv")
}
write.csv(validation_summary, file.path(output_dir, "validation_summary.csv"), row.names = FALSE, na = "")

cat("Black rot v4 descriptive outputs regenerated successfully.\n")
cat("R version:", R.version.string, "\n")
cat("Output directory:", normalizePath(output_dir, mustWork = TRUE), "\n")
cat("Program-versus-UTC contrasts:", nrow(program_contrasts), "rows; 8 unresolved 2022 plot rows (2 program IDs) excluded.\n")
cat("Program summary:", nrow(program_summary), "year-specific programs; unresolved rows retained descriptively and flagged.\n")
cat("Mixed architecture summary:", nrow(architecture_summary), "observed year-architecture cells across 15 mixed programs.\n")
cat("Validation checks:", nrow(validation_summary), "passed; no endpoint reconstruction or classification discrepancies detected.\n")
