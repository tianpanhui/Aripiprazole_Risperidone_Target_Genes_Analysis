#!/usr/bin/env Rscript

# Rare variant carrier phenotype association analysis
#
# This script creates a gene-level rare-variant carrier status from a PLINK
# .raw dosage file, merges SPARK phenotype/ancestry files, performs phenotype
# association analyses, and generates carrier vs non-carrier violin plots.
#
# Example:
# Rscript rare_variant_phenotype_analysis.R \
#   --gene DRD3 \
#   --raw DRD3.rare.raw.raw.txt \
#   --master SPARK.iWES_v2.mastertable.2023_01.xlsx \
#   --core core_descriptive_variables_production-2023-07-21.csv \
#   --ancestry SPARK.iWES_v2.ancestry.2023_01.xlsx \
#   --dcdq dcdq_production-2023-07-21.csv \
#   --outdir results/DRD3

suppressPackageStartupMessages({
  library(optparse)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(broom)
  library(readr)
})

option_list <- list(
  make_option(c("--gene"), type = "character", default = "DRD3",
              help = "Gene name used in output labels [default: %default]"),
  make_option(c("--raw"), type = "character", default = "DRD3.rare.raw.raw.txt",
              help = "PLINK .raw dosage file for rare variants"),
  make_option(c("--master"), type = "character", default = "SPARK.iWES_v2.mastertable.2023_01.xlsx",
              help = "SPARK master table"),
  make_option(c("--core"), type = "character", default = "core_descriptive_variables_production-2023-07-21.csv",
              help = "SPARK core phenotype CSV"),
  make_option(c("--ancestry"), type = "character", default = "SPARK.iWES_v2.ancestry.2023_01.xlsx",
              help = "SPARK ancestry XLSX"),
  make_option(c("--dcdq"), type = "character", default = "dcdq_production-2023-07-21.csv",
              help = "DCDQ phenotype CSV; optional but recommended"),
  make_option(c("--outdir"), type = "character", default = "results/DRD3",
              help = "Output directory [default: %default]"),
  make_option(c("--carrier-threshold"), type = "double", default = 0.9,
              help = "Dosage sum threshold for carrier status [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

standardize_spid <- function(df) {
  names(df)[1] <- "spid"
  df
}

safe_glm <- function(formula, data, family = gaussian()) {
  tryCatch(
    glm(formula, data = data, family = family),
    error = function(e) {
      message("Model failed: ", deparse(formula), " | ", e$message)
      NULL
    }
  )
}

extract_model_term <- function(model, term = "mutation") {
  if (is.null(model)) return(NULL)
  broom::tidy(model, conf.int = TRUE) %>%
    filter(.data$term == term)
}

format_p_fdr <- function(p) {
  paste0("p-FDR = ", formatC(p, format = "f", digits = 3))
}

# -----------------------------------------------------------------------------
# 1. Create gene-level carrier status from PLINK .raw file
# -----------------------------------------------------------------------------

raw_df <- read.table(opt$raw, header = TRUE, check.names = FALSE)
if (ncol(raw_df) <= 6) {
  stop("The .raw file must contain the first six PLINK columns plus variant dosage columns.")
}

carrier_df <- raw_df %>%
  mutate(variant_dosage_sum = rowSums(across(-(1:6)), na.rm = TRUE),
         mutation = if_else(variant_dosage_sum > opt$carrier_threshold, 1L, 0L)) %>%
  transmute(spid = .data[[2]], variant_dosage_sum, mutation)

write_csv(carrier_df, file.path(opt$outdir, paste0(opt$gene, "_rare_variant_carrier_status.csv")))

# -----------------------------------------------------------------------------
# 2. Load and merge phenotype/ancestry files
# -----------------------------------------------------------------------------

master <- read_excel(opt$master)
core <- read.csv(opt$core, header = TRUE) %>% standardize_spid()
ancestry <- read_excel(opt$ancestry)

merged_with_asd <- master %>%
  inner_join(carrier_df, by = "spid")

merged_df <- core %>%
  inner_join(merged_with_asd, by = "spid") %>%
  inner_join(ancestry, by = "spid")

if (file.exists(opt$dcdq)) {
  dcdq <- read.csv(opt$dcdq, header = TRUE) %>% standardize_spid()
  merged_df <- merged_df %>% inner_join(dcdq, by = "spid")
}

# Recode variables when present. In SPARK tables used here, sex/asd may be coded 1/2.
if ("asd" %in% names(merged_with_asd)) {
  merged_with_asd <- merged_with_asd %>% mutate(asd_binary = if_else(asd == 1, 0L, 1L))
}
if ("sex" %in% names(merged_with_asd)) {
  merged_with_asd <- merged_with_asd %>% mutate(sex_binary = if_else(sex == 1, 0L, 1L))
}
if ("asd.y" %in% names(merged_df)) {
  merged_df <- merged_df %>% mutate(asd_binary = if_else(asd.y == 1, 0L, 1L))
}
if ("sex.y" %in% names(merged_df)) {
  merged_df <- merged_df %>% mutate(sex_binary = if_else(sex.y == 1, 0L, 1L))
} else if ("sex" %in% names(merged_df)) {
  merged_df <- merged_df %>% mutate(sex_binary = if_else(sex == 1, 0L, 1L))
}

merged_df <- merged_df %>%
  mutate(mutation_fac = factor(mutation, levels = c(0, 1), labels = c("Non-carrier", "Carrier")))

write_csv(merged_df, file.path(opt$outdir, paste0(opt$gene, "_merged_analysis_dataset.csv")))

# -----------------------------------------------------------------------------
# 3. Association analyses
# -----------------------------------------------------------------------------

continuous_phenotypes <- c(
  "fsiq", "viq", "nviq",
  "scq_total_final_score", "rbsr_total_final_score",
  "vineland_abc_ss_latest", "used_words_age_mos", "walked_age_mos", "diagnosis_age",
  "control_during_movement", "fine_motor_handwriting", "general_coordination",
  "final_score", "q07_printing_writing_drawing_fast"
)
continuous_phenotypes <- intersect(continuous_phenotypes, names(merged_df))

binary_phenotypes <- c(
  "dcdq_dcd", "regress_lang_y_n", "regress_other_y_n",
  "current_depend_adult", "cognitive_impairment_latest"
)
binary_phenotypes <- intersect(binary_phenotypes, names(merged_df))

covariates_full <- c("sex_binary", "age_m", "superclass")
covariates_basic <- c("sex_binary", "age_m")

run_models <- function(data, phenotypes, family, covariates) {
  results <- lapply(phenotypes, function(y) {
    covs <- covariates[covariates %in% names(data)]
    formula_txt <- paste(y, "~ mutation", if (length(covs) > 0) paste("+", paste(covs, collapse = " + ")) else "")
    model <- safe_glm(as.formula(formula_txt), data = data, family = family)
    res <- extract_model_term(model, "mutation")
    if (is.null(res) || nrow(res) == 0) return(NULL)
    res %>% mutate(phenotype = y, model = formula_txt, .before = 1)
  })
  bind_rows(results) %>%
    mutate(p_FDR = p.adjust(p.value, method = "fdr"))
}

continuous_results <- run_models(merged_df, continuous_phenotypes, gaussian(), covariates_full)
binary_results <- run_models(merged_df, binary_phenotypes, binomial(), covariates_full) %>%
  mutate(OR = exp(estimate), OR_low = exp(conf.low), OR_high = exp(conf.high))

write_csv(continuous_results, file.path(opt$outdir, paste0(opt$gene, "_continuous_trait_associations.csv")))
write_csv(binary_results, file.path(opt$outdir, paste0(opt$gene, "_binary_trait_associations.csv")))

if (all(c("asd_binary", "sex_binary", "age_m") %in% names(merged_with_asd))) {
  asd_model <- safe_glm(asd_binary ~ mutation + sex_binary + age_m,
                        data = merged_with_asd, family = binomial())
  asd_result <- broom::tidy(asd_model, conf.int = TRUE) %>%
    mutate(OR = exp(estimate), OR_low = exp(conf.low), OR_high = exp(conf.high))
  write_csv(asd_result, file.path(opt$outdir, paste0(opt$gene, "_ASD_status_logistic_model.csv")))
}

# -----------------------------------------------------------------------------
# 4. Figure: five core ASD-related phenotypes
# -----------------------------------------------------------------------------

plot_phenotypes <- c("fsiq", "viq", "nviq", "scq_total_final_score", "rbsr_total_final_score")
plot_phenotypes <- intersect(plot_phenotypes, names(merged_df))

if (length(plot_phenotypes) > 0) {
  long_df <- merged_df %>%
    select(spid, mutation_fac, all_of(plot_phenotypes)) %>%
    pivot_longer(cols = all_of(plot_phenotypes), names_to = "Phenotype", values_to = "Score") %>%
    filter(!is.na(Score), !is.na(mutation_fac)) %>%
    mutate(Phenotype = factor(Phenotype, levels = plot_phenotypes))

  stats_long <- long_df %>%
    group_by(Phenotype, mutation_fac) %>%
    summarise(
      mean = mean(Score, na.rm = TRUE),
      sd = sd(Score, na.rm = TRUE),
      n = sum(!is.na(Score)),
      se = sd / sqrt(n),
      ci_low = mean - 1.96 * se,
      ci_high = mean + 1.96 * se,
      .groups = "drop"
    )

  p_label_df <- long_df %>%
    group_by(Phenotype) %>%
    summarise(
      p_raw = tryCatch(t.test(Score ~ mutation_fac)$p.value, error = function(e) NA_real_),
      y_max = max(Score, na.rm = TRUE),
      y_min = min(Score, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      p_FDR = p.adjust(p_raw, method = "fdr"),
      stars = case_when(
        p_FDR < 0.001 ~ "***",
        p_FDR < 0.01 ~ "**",
        p_FDR < 0.05 ~ "*",
        TRUE ~ ""
      ),
      label = if_else(stars == "", format_p_fdr(p_FDR), paste(format_p_fdr(p_FDR), stars)),
      range = if_else(y_max > y_min, y_max - y_min, abs(y_max) * 0.1 + 1),
      y_pos = y_max + 0.10 * range,
      y_pos_line = y_max + 0.05 * range,
      sig = p_FDR < 0.05
    )

  p_multi <- ggplot(long_df, aes(x = mutation_fac, y = Score)) +
    geom_violin(aes(fill = mutation_fac), trim = TRUE, alpha = 0.4,
                linewidth = 0.3, color = NA) +
    geom_boxplot(width = 0.15, fill = "white", color = "black",
                 alpha = 0.9, outlier.shape = NA, linewidth = 0.4) +
    geom_segment(data = p_label_df,
                 aes(x = 1, xend = 2, y = y_pos_line, yend = y_pos_line),
                 inherit.aes = FALSE, linewidth = 0.7, colour = "black") +
    geom_text(data = filter(p_label_df, !sig),
              aes(x = 1.5, y = y_pos, label = label), inherit.aes = FALSE,
              size = 4.2, fontface = "bold", colour = "black") +
    geom_text(data = filter(p_label_df, sig),
              aes(x = 1.5, y = y_pos, label = label), inherit.aes = FALSE,
              size = 4.2, fontface = "bold", colour = "red3") +
    scale_fill_manual(values = c("#9AD0A5", "#7DA2D5")) +
    labs(
      x = "",
      y = "Score",
      title = paste0("Phenotypic Differences Between Carriers and Non-carriers of Rare ", opt$gene, " Mutations")
    ) +
    facet_wrap(~ Phenotype, scales = "free_y", nrow = 1) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",
      strip.text = element_text(size = 14, face = "bold"),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 15),
      axis.text = element_text(size = 12, color = "black")
    )

  ggsave(file.path(opt$outdir, paste0(opt$gene, "_carrier_phenotype_violin.pdf")),
         p_multi, width = 14, height = 4.5)
  ggsave(file.path(opt$outdir, paste0(opt$gene, "_carrier_phenotype_violin.png")),
         p_multi, width = 14, height = 4.5, dpi = 300)

  write_csv(stats_long, file.path(opt$outdir, paste0(opt$gene, "_plot_summary_statistics.csv")))
  write_csv(p_label_df, file.path(opt$outdir, paste0(opt$gene, "_plot_t_test_FDR.csv")))
}

message("Analysis complete. Results written to: ", opt$outdir)
