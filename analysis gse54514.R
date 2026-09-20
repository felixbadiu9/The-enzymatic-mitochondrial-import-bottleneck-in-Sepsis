# GSE54514 LONGITUDINAL ANALYSIS 


library(GEOquery)
library(Biobase)
library(tidyverse)
library(lme4)
library(lmerTest)   
library(ggplot2)

options(na.print = "NA")


safe_print <- function(x, ...) {
  tryCatch(
    print(x, ...),
    error = function(e) {
      message("  NOTE: console print failed (", e$message,
              ") -- this is a DISPLAY issue only, the data/file was still saved correctly.")
      invisible(NULL)
    }
  )
}

GENES <- c("GCLC", "GCLM", "NFE2L2", "PINK1", "PRKN", "NLRP3", "GSDMD", "BNIP3L", "FUNDC1")


output_dir <- "C:/Users/andre/Downloads/output"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
}
setwd(output_dir)

#LOGGING
log_file <- file("gse54514_longitudinal_execution.log", open = "wt")
sink(log_file, split = TRUE)
sink(log_file, type = "message")
on.exit({
  sink(type = "message")
  sink()
  close(log_file)
  message("Logging safely closed. Outputs should be in: ", getwd())
})

message("Logging started at: ", Sys.time())
message("All outputs will be saved to: ", output_dir)


# DOWNLOAD + MAP GENES


message("Downloading GSE54514...")
gse  <- tryCatch(getGEO("GSE54514", GSEMatrix = TRUE, AnnotGPL = TRUE, destdir = "."),
                 error = function(e) stop("Download failed: ", e$message))
eset <- if (is.list(gse)) gse[[1]] else gse
message("  Dimensions: ", nrow(eset), " probes x ", ncol(eset), " samples")

map_probes_to_genes <- function(eset, genes) {
  fdata <- fData(eset)
  symbol_cols <- c("Gene Symbol", "Gene symbol", "GENE_SYMBOL", "gene_assignment",
                   "Symbol", "GeneSymbol", "Gene_Symbol", "SYMBOL", "gene_symbol")
  sym_col <- intersect(symbol_cols, colnames(fdata))[1]
  if (is.na(sym_col)) stop("No gene symbol column found in fData(eset).")
  message("  Using gene symbol column: '", sym_col, "'")
  
  fdata$gene_symbol_clean <- trimws(gsub(" ///.*", "", fdata[[sym_col]]))
  expr_mat <- as.data.frame(exprs(eset))
  expr_mat$gene_symbol <- fdata$gene_symbol_clean[match(rownames(expr_mat), rownames(fdata))]
  
  # Ensure max() is used for mapping Illumina probes
  expr_mat %>%
    dplyr::filter(gene_symbol %in% genes) %>%
    dplyr::group_by(gene_symbol) %>%
    dplyr::summarise(dplyr::across(dplyr::where(is.numeric), \(x) max(x, na.rm = TRUE)), .groups = "drop") %>%
    tibble::column_to_rownames("gene_symbol")
}

expr_genes <- map_probes_to_genes(eset, GENES)
found_genes   <- intersect(GENES, rownames(expr_genes))
missing_genes <- setdiff(GENES, rownames(expr_genes))
message("  Genes found:     ", paste(found_genes, collapse = ", "))
if (length(missing_genes) > 0) message("  Genes NOT found: ", paste(missing_genes, collapse = ", "))

if (length(found_genes) == 0) stop("FATAL ERROR: No target genes found. Halting script.")


# extract patient_id, day, condition for every sample


pdata     <- pData(eset)
raw_title <- as.character(pdata$title)
char_cols <- grep("characteristics", colnames(pdata), value = TRUE)
search_text <- raw_title
if (length(char_cols) > 0) {
  for (cc in char_cols) search_text <- paste(search_text, as.character(pdata[[cc]]), sep = " | ")
}

condition_std <- dplyr::case_when(
  grepl("control|healthy|normal|volunteer", raw_title, ignore.case = TRUE) ~ "Control",
  grepl("sepsis|septic|shock|survivor|nonsurvivor", raw_title, ignore.case = TRUE) ~ "Sepsis",
  TRUE ~ "Other"
)

day <- as.integer(stringr::str_match(search_text, "(?i)(?:day|time\\s*point)[_\\s:=]*(\\d+)")[, 2])
patient_id <- stringr::str_match(search_text, "(?i)(?:patient|subject|\\bid)[_\\s:=]*(\\d+)")[, 2]
missing_idx <- which(is.na(patient_id))

if (length(missing_idx) > 0) {
  fallback <- trimws(gsub("(?i)[\\s_,-]*(day|time\\s*point)[_\\s:=]*\\d+.*$", "", raw_title[missing_idx]))
  fallback[fallback == ""] <- NA
  patient_id[missing_idx] <- fallback
}

sample_meta <- data.frame(
  sample_id     = colnames(expr_genes),
  condition_std = condition_std,
  day           = day,
  patient_id    = patient_id,
  stringsAsFactors = FALSE
) %>%
  dplyr::filter(condition_std %in% c("Sepsis", "Control"))

# Drop unresolved samples to prevent phantom patients
sample_meta <- sample_meta %>%
  dplyr::filter(!is.na(patient_id), patient_id != "", !is.na(day))

ROGUE_SAMPLES <- c("GSM1317938")  # sepsis_nonsurvivor, Day_4, ID=20 -- orphan, contradicts documented n=9

n_before <- nrow(sample_meta)
sample_meta <- sample_meta %>% dplyr::filter(!sample_id %in% ROGUE_SAMPLES)
n_removed <- n_before - nrow(sample_meta)

message("Rogue-sample exclusion:")
message("  Removed ", n_removed, " sample(s): ", paste(ROGUE_SAMPLES, collapse = ", "))
if (n_removed == 0) {
  warning("Expected to remove 1 rogue sample (GSM1317938) but removed 0. ",
          "Check that sample IDs in this GEO pull match the accessions used above ",
          "(re-verify against the current GEO record for GSE54514).")
}

# CONFIRMATION
n_sepsis_check  <- length(unique(sample_meta$patient_id[sample_meta$condition_std == "Sepsis"]))
n_control_check <- length(unique(sample_meta$patient_id[sample_meta$condition_std == "Control"]))

message("  Post-exclusion unique patients -> Sepsis: ", n_sepsis_check,
        " | Control: ", n_control_check,
        "  (expected: Sepsis = 35, Control = 18)")

if (n_sepsis_check != 35 || n_control_check != 18) {
  warning("Patient counts (Sepsis=", n_sepsis_check, ", Control=", n_control_check,
          ") do NOT match the documented GSE54514 design (35/18). ",
          "Do not proceed to publication-facing tables/figures until this is resolved -- ",
          "check patient_id extraction and the rogue-sample exclusion list above.")
} else {
  message("  OK: cohort matches documented design (35 septic / 18 control).")
}

message("Samples retained for longitudinal model: ", nrow(sample_meta))
message("Day distribution:")
print(table(sample_meta$condition_std, sample_meta$day))

if (nrow(sample_meta) == 0) stop("FATAL ERROR: No valid samples left after metadata filtering.")


# SECTION 3: LONG-FORMAT DATA + MIXED-EFFECTS MODEL PER GENE


long_df <- purrr::map_dfr(found_genes, function(g) {
  sample_meta %>%
    dplyr::mutate(
      gene       = g,
      expression = as.numeric(expr_genes[g, sample_id])
    )
})

results_list <- list()
model_objects <- list()


for (g in found_genes) {
  df_g <- long_df %>% dplyr::filter(gene == g)
  df_g$condition_std <- factor(df_g$condition_std, levels = c("Control", "Sepsis"))
  
  fit <- try(lmerTest::lmer(expression ~ condition_std + day + (1 | patient_id), data = df_g), silent = TRUE)
  
  if (inherits(fit, "try-error")) {
    message("  ", g, ": Model failed to fit -- ", fit)
    next
  }
  
  model_objects[[g]] <- fit
  s <- summary(fit)$coefficients
  cond_row <- grep("^condition_std", rownames(s))[1]
  
  if (is.na(cond_row)) next
  
  singular <- lme4::isSingular(fit)
  
  vc <- as.data.frame(lme4::VarCorr(fit))
  var_intercept <- vc$vcov[vc$grp == "patient_id" & is.na(vc$var2)]
  var_residual  <- vc$vcov[vc$grp == "Residual"]
  icc <- if (length(var_intercept) == 1 && length(var_residual) == 1 && (var_intercept + var_residual) > 0) {
    var_intercept / (var_intercept + var_residual)
  } else NA_real_
  
  n_patients_total <- length(unique(df_g$patient_id))
  avg_cluster_size <- nrow(df_g) / n_patients_total
  design_effect <- if (!is.na(icc)) 1 + (avg_cluster_size - 1) * icc else NA_real_
  effective_n   <- if (!is.na(design_effect)) round(nrow(df_g) / design_effect, 1) else NA_real_
  
  fit_int <- try(lmerTest::lmer(expression ~ condition_std * day + (1 | patient_id), data = df_g), silent = TRUE)
  interaction_p <- NA_real_
  if (!inherits(fit_int, "try-error")) {
    s_int <- summary(fit_int)$coefficients
    int_row <- grep(":", rownames(s_int))[1]
    if (!is.na(int_row)) interaction_p <- s_int[int_row, "Pr(>|t|)"]
  }
  
  results_list[[g]] <- data.frame(
    gene = g,
    n_obs = nrow(df_g),
    n_patients_sepsis  = length(unique(df_g$patient_id[df_g$condition_std == "Sepsis"])),
    n_patients_control = length(unique(df_g$patient_id[df_g$condition_std == "Control"])),
    estimate = round(s[cond_row, "Estimate"], 4),      
    se       = round(s[cond_row, "Std. Error"], 4),
    df       = round(s[cond_row, "df"], 1),
    t_value  = round(s[cond_row, "t value"], 3),
    p_value  = s[cond_row, "Pr(>|t|)"],
    converged = TRUE,
    singular_fit = singular,
    icc = round(icc, 3),
    avg_days_per_patient = round(avg_cluster_size, 2),
    effective_n = effective_n,
    interaction_p_value = interaction_p,
    stringsAsFactors = FALSE
  )
}

if (length(results_list) == 0) stop("FATAL ERROR: No models converged.")

results_df <- dplyr::bind_rows(results_list) %>%
  dplyr::mutate(
    p_adjusted = p.adjust(p_value, method = "BH"),
    sig_label  = dplyr::case_when(
      is.na(p_adjusted)  ~ "NA",
      p_adjusted < 0.001 ~ "***",
      p_adjusted < 0.01  ~ "**",
      p_adjusted < 0.05  ~ "*",
      TRUE               ~ "ns"
    ),
    interaction_sig = dplyr::case_when(
      is.na(interaction_p_value) ~ NA_character_,
      interaction_p_value < 0.05 ~ "Trajectories diverge over time (p<0.05, uncorrected)",
      TRUE ~ "No evidence trajectories diverge over time"
    ),
    direction = ifelse(estimate > 0, "UP in Sepsis", "DOWN in Sepsis")
  )

# Save the CSV immediately before printing
readr::write_csv(results_df, "gse54514_longitudinal_results.csv")
message("Saved: gse54514_longitudinal_results.csv")

message("LONGITUDINAL MIXED-MODEL RESULTS:")
safe_print(tibble::as_tibble(results_df) %>%
             dplyr::select(gene, n_obs, n_patients_sepsis, n_patients_control,
                           estimate, p_value, p_adjusted, sig_label, direction,
                           icc, avg_days_per_patient, effective_n,
                           singular_fit, interaction_p_value, interaction_sig),
           n = Inf, width = Inf)


# figures


message("\Generating figures...")

if (nrow(long_df) > 0) {
  tryCatch({
    p1 <- ggplot(long_df, aes(x = day, y = expression, group = patient_id, color = condition_std)) +
      geom_line(alpha = 0.25, linewidth = 0.4) +
      geom_point(alpha = 0.35, size = 1) +
      geom_smooth(aes(group = condition_std), method = "loess", se = TRUE, linewidth = 1.1) +
      facet_wrap(~ gene, scales = "free_y") +
      scale_color_manual(values = c("Sepsis" = "#e74c3c", "Control" = "#3498db")) +
      scale_x_continuous(breaks = sort(unique(long_df$day))) +
      labs(title = "GSE54514: Per-Patient Expression Trajectories Over Time",
           subtitle = paste0("Thin lines = individual patients; thick line = LOESS trend per group. ",
                             nrow(sample_meta), " samples used (rogue orphan sample excluded)."),
           x = "Day", y = "Expression (log2 units)", color = "Condition") +
      theme_bw(base_size = 11) +
      theme(plot.title = element_text(face = "bold"))
    ggsave("figure1_longitudinal_trajectories.png", p1, width = 12, height = 8, dpi = 150)
    message("  Saved: figure1_longitudinal_trajectories.png")
  }, error = function(e) message("  Figure 1 FAILED: ", e$message))
}

if (nrow(results_df) > 0) {
  tryCatch({
    forest_df <- results_df %>%
      dplyr::filter(!is.na(estimate)) %>%
      dplyr::mutate(
        ci_lo = estimate - 1.96 * se,
        ci_hi = estimate + 1.96 * se,
        pathway = dplyr::case_when(
          gene %in% c("GCLC","GCLM") ~ "GSH Synthesis",
          gene == "NFE2L2"            ~ "Nrf2",
          gene %in% c("PINK1","PRKN", "BNIP3L", "FUNDC1") ~ "Mitophagy",
          gene %in% c("NLRP3","GSDMD") ~ "Inflammasome/Pyroptosis"
        )
      )
    
    p2 <- ggplot(forest_df, aes(x = estimate, y = reorder(gene, estimate), color = pathway)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.15, linewidth = 0.9) +
      geom_point(size = 3) +
      labs(title    = "GSE54514 Longitudinal Model: Sepsis vs Control Effect (adjusted for day)",
           subtitle = "Point = mixed-model estimate; bars = 95% CI. Patient random intercept; rogue orphan sample excluded.",
           x = "Estimate (Sepsis - Control, log2 units)", y = NULL, color = "Pathway") +
      theme_bw(base_size = 11) +
      theme(plot.title = element_text(face = "bold"))
    ggsave("figure2_longitudinal_forest.png", p2, width = 8, height = 5, dpi = 150)
    message("  Saved: figure2_longitudinal_forest.png")
  }, error = function(e) message("  Figure 2 FAILED: ", e$message))
}


# SECTION 5: SENSITIVITY CHECK -- ROGUE ENTRY ON DAY 4 


message("Running leave-one-patient-out sensitivity check...")

lopo_results <- list()

for (g in found_genes) {
  df_g <- long_df %>% dplyr::filter(gene == g)
  df_g$condition_std <- factor(df_g$condition_std, levels = c("Control", "Sepsis"))
  all_patients <- unique(df_g$patient_id)
  
  for (pid in all_patients) {
    df_sub <- df_g %>% dplyr::filter(patient_id != pid)
    if (length(unique(df_sub$condition_std)) < 2) next
    
    fit_sub <- try(lmerTest::lmer(expression ~ condition_std + day + (1 | patient_id), data = df_sub), silent = TRUE)
    if (inherits(fit_sub, "try-error")) next
    
    s_sub <- summary(fit_sub)$coefficients
    cond_row_sub <- grep("^condition_std", rownames(s_sub))[1]
    if (is.na(cond_row_sub)) next
    
    lopo_results[[paste(g, pid, sep = "_")]] <- data.frame(
      gene = g,
      excluded_patient = pid,
      n_days_excluded  = sum(df_g$patient_id == pid),
      estimate_without = round(s_sub[cond_row_sub, "Estimate"], 4),
      p_value_without  = s_sub[cond_row_sub, "Pr(>|t|)"],
      stringsAsFactors = FALSE
    )
  }
}

if (length(lopo_results) > 0) {
  lopo_df <- dplyr::bind_rows(lopo_results) %>%
    dplyr::left_join(
      results_df %>% dplyr::select(gene, estimate_full = estimate, p_value_full = p_value),
      by = "gene"
    ) %>%
    dplyr::mutate(
      estimate_shift     = round(estimate_without - estimate_full, 4),
      flips_significance = (p_value_full < 0.05) != (p_value_without < 0.05)
    )
  
  readr::write_csv(lopo_df, "gse54514_longitudinal_LOPO_sensitivity.csv")
  message("Saved: gse54514_longitudinal_LOPO_sensitivity.csv")
} else {
  message("No LOPO models converged.")
}

message("Done. Check the log file if anything failed.")