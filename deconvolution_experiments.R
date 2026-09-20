

library(xCell)
library(GEOquery)
library(Biobase)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)
library(writexl)

GENES <- c("GCLC", "GCLM", "NFE2L2", "NLRP3", "GSDMD", "PINK1", "PRKN", "BNIP3L", "BNIP3", "FUNDC1","SLC25A39","SLC25A40","AFG3L2")

datasets <- c("GSE54514", "GSE13904", "GSE28750")

CELLS_OF_INTEREST <- c(
  "Monocytes", "Neutrophils", "NK cells", "CD8+ T-cells", "CD4+ T-cells",
  "B-cells", "Dendritic cells", "Macrophages M1", "Macrophages M2"
)

# Output goes straight to Downloads.
output_dir <- "C:/Users/andre/Downloads"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
setwd(output_dir)
message("All outputs will be saved to: ", output_dir)


# SECTION 2: HELPER FUNCTIONS


map_probes_to_genes <- function(eset) {
  fdata <- fData(eset)
  symbol_cols <- c("Gene Symbol", "Gene symbol", "GENE_SYMBOL", "gene_assignment",
                   "Symbol", "GeneSymbol", "Gene_Symbol", "SYMBOL", "gene_symbol")
  sym_col <- intersect(symbol_cols, colnames(fdata))[1]
  
  if (is.na(sym_col)) {
    message("  WARNING: No gene symbol column found. Columns: ",
            paste(colnames(fdata), collapse = ", "))
    return(NULL)
  }
  
  message("  Using gene symbol column: '", sym_col, "'")
  fdata$gene_symbol_clean <- gsub(" ///.*", "", fdata[[sym_col]])
  fdata$gene_symbol_clean <- trimws(fdata$gene_symbol_clean)
  
  expr_mat <- as.data.frame(exprs(eset))
  expr_mat$gene_symbol <- fdata$gene_symbol_clean[match(rownames(expr_mat), rownames(fdata))]
  
  expr_all <- expr_mat %>%
    dplyr::filter(!is.na(gene_symbol), gene_symbol != "") %>%
    dplyr::group_by(gene_symbol) %>%
    dplyr::summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
    tibble::column_to_rownames("gene_symbol")
  
  return(expr_all)
}

# ---- explicitly flag unresolved patient IDs instead of letting them
#      silently collapse into one NA group downstream --------------------
extract_patient_day_gse54514 <- function(eset) {
  pdata <- pData(eset)
  raw_title <- as.character(pdata$title)
  
  char_cols <- grep("characteristics", colnames(pdata), value = TRUE)
  search_text <- raw_title
  if (length(char_cols) > 0) {
    for (cc in char_cols) {
      search_text <- paste(search_text, as.character(pdata[[cc]]), sep = " | ")
    }
  }
  
  day <- as.integer(stringr::str_match(search_text, "(?i)(?:day|time\\s*point)[_\\s:=]*(\\d+)")[, 2])
  
  patient_id <- stringr::str_match(search_text, "(?i)(?:patient|subject|\\bid)[_\\s:=]*(\\d+)")[, 2]
  
  missing_idx <- which(is.na(patient_id))
  if (length(missing_idx) > 0) {
    fallback <- gsub("(?i)[\\s_,-]*(day|time\\s*point)[_\\s:=]*\\d+.*$", "", raw_title[missing_idx])
    fallback <- trimws(fallback)
    fallback[fallback == ""] <- NA
    patient_id[missing_idx] <- fallback
  }
  
  still_missing <- is.na(patient_id) | patient_id == ""
  if (any(still_missing)) {
    message("  WARNING: could not resolve patient_id for ", sum(still_missing),
            " sample(s) -- these will be DROPPED before deduplication to avoid",
            " forming a phantom patient group:")
    message("    ", paste(raw_title[still_missing], collapse = " | "))
  }
  
  data.frame(
    condition_std = dplyr::case_when(
      grepl("control|healthy|normal|volunteer", raw_title, ignore.case = TRUE) ~ "Control",
      grepl("sepsis|septic|shock|survivor|nonsurvivor", raw_title, ignore.case = TRUE) ~ "Sepsis",
      TRUE ~ "Other"
    ),
    day = day,
    patient_id = patient_id,
    still_missing_id = still_missing,
    stringsAsFactors = FALSE
  )
}

extract_groups <- function(eset, dataset_id) {
  pdata <- pData(eset)
  raw_title <- as.character(pdata$title)
  
  condition_std <- dplyr::case_when(
    grepl("control|healthy|normal|volunteer", raw_title, ignore.case = TRUE) ~ "Control",
    grepl("sepsis|septic|shock|survivor|nonsurvivor", raw_title, ignore.case = TRUE) ~ "Sepsis",
    TRUE ~ "Other"
  )
  
  char_cols <- grep("characteristics", colnames(pdata), value = TRUE)
  search_text <- raw_title
  if (length(char_cols) > 0) {
    for (cc in char_cols) {
      search_text <- paste(search_text, as.character(pdata[[cc]]), sep = " | ")
    }
  }
  
  day <- as.integer(stringr::str_match(search_text, "(?i)(?:day|time\\s*point|time)[_\\s:=]*(\\d+)")[, 2])
  
  # Force Healthy Controls to Day 1
  day[condition_std == "Control" & is.na(day)] <- 1
  
  patient_id <- stringr::str_match(search_text, "(?i)(?:patient|subject|\\bid)[_\\s:=]*(\\w+)")[, 2]
  
  # Fallback for missing IDs
  missing_idx <- which(is.na(patient_id))
  if (length(missing_idx) > 0) {
    fallback <- gsub("(?i)[\\s_,-]*(day|time\\s*point|24h|48h|72h|96h|120h).*$", "", raw_title[missing_idx])
    fallback <- trimws(fallback)
    fallback[fallback == ""] <- NA
    patient_id[missing_idx] <- fallback
  }
  
  still_missing <- is.na(patient_id) | patient_id == ""
  
  if (any(still_missing)) {
    message("  WARNING: could not resolve patient_id for ", sum(still_missing), " sample(s) in ", dataset_id)
  }
  
  df <- data.frame(
    condition_std = condition_std,
    day = day,
    patient_id = patient_id,
    still_missing_id = still_missing,
    stringsAsFactors = FALSE
  )
  
  message("  Condition distribution (pre-filtering) for ", dataset_id, ":")
  print(table(df$condition_std))
  return(df)
}


#day-1 filter
filter_to_day1 <- function(expr_mat, groups_df, dataset_id) {
  n_before <- ncol(expr_mat)
  df <- groups_df %>% dplyr::mutate(idx = seq_len(n_before))
  
  n_unresolved <- sum(df$still_missing_id, na.rm = TRUE)
  if (n_unresolved > 0) {
    message("  ", dataset_id, ": dropping ", n_unresolved,
            " sample(s) with unresolved patient_id before deduplication")
    df <- df %>% dplyr::filter(!still_missing_id)
  }
  
  # Dedup Ctrls
  df_ctrl_kept <- df %>%
    dplyr::filter(condition_std == "Control") %>%
    dplyr::filter(day == 1 | is.na(day)) %>%
    dplyr::group_by(patient_id) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
  
  # Dedup Sepsis
  df_sepsis_kept <- df %>%
    dplyr::filter(condition_std == "Sepsis") %>%
    dplyr::filter(day == 1 | is.na(day)) %>% # Accepts NA day for GSE13904 if missing, but slices by ID
    dplyr::group_by(patient_id) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
  
  df_keep <- dplyr::bind_rows(df_ctrl_kept, df_sepsis_kept) %>% dplyr::arrange(idx)
  keep_idx <- seq_len(n_before) %in% df_keep$idx
  
  message("  ", dataset_id, " day-1 filter applied:")
  message("    Control: -> ", nrow(df_ctrl_kept), " unique cross-sectional controls")
  message("    Sepsis:  -> ", nrow(df_sepsis_kept), " unique day-1 patients")
  
  list(expr = expr_mat[, keep_idx, drop = FALSE], groups = groups_df$condition_std[keep_idx])
}


zscore_within_dataset <- function(expr_mat) {
  z <- t(scale(t(as.matrix(expr_mat))))
  z[is.nan(z)] <- NA   # genes with zero variance in this dataset -> NA, not NaN
  z
}

# SECTION 3: DOWNLOAD, PROCESS, AND DECONVOLVE EACH DATASET


results_list     <- list()
xcell_scores_all <- list()
gene_expr_all    <- list()
groups_all       <- list()

for (gse_id in datasets) {
  message("", strrep("=", 60))
  message("Processing ", gse_id, "...")
  message(strrep("=", 60))
  
  tryCatch({
    gse  <- getGEO(gse_id, GSEMatrix = TRUE, AnnotGPL = TRUE, destdir = ".")
    eset <- if (is.list(gse)) gse[[1]] else gse
    message("  Dimensions: ", nrow(eset), " probes x ", ncol(eset), " samples")
    
    expr_all <- map_probes_to_genes(eset)
    if (is.null(expr_all)) {
      message("  ERROR: Gene mapping failed. Skipping.")
      next
    }
    message("  Unique genes after mapping: ", nrow(expr_all))
    
    groups_df <- extract_groups(eset, gse_id)
    
    filt <- filter_to_day1(expr_all, groups_df, gse_id)
    expr_all <- filt$expr
    groups   <- filt$groups
    
    groups_all[[gse_id]] <- groups
    
    keep_idx  <- which(groups %in% c("Sepsis", "Control"))
    expr_filt <- expr_all[, keep_idx, drop = FALSE]
    grp_filt  <- groups[keep_idx]
    
    message("  Samples retained (Sepsis + Control): ", length(keep_idx))
    message("  [Platform note] This dataset's samples/scores are never pooled",
            " with any other dataset in any statistical test below --",
            " every wilcox.test()/cor.test() runs within ", gse_id, " only.")
    
    message("  Running xCell deconvolution...")
    expr_matrix <- as.matrix(expr_filt)   # raw/platform-normalized values -> xCell (NOT z-scored)
    
    xcell_raw <- xCellAnalysis(expr_matrix, parallel.sz = 1)
    
    xcell_df <- as.data.frame(t(xcell_raw))
    xcell_df$sample_id <- colnames(expr_matrix)
    xcell_df$condition <- grp_filt
    xcell_df$dataset   <- gse_id
    
    xcell_scores_all[[gse_id]] <- xcell_df
    message("  xCell complete. Cell types scored: ", ncol(xcell_raw))
    
    found_genes <- intersect(GENES, rownames(expr_filt))
    message("  Target genes found: ", paste(found_genes, collapse = ", "))
    missing_target_genes <- setdiff(GENES, rownames(expr_filt))
    if (length(missing_target_genes) > 0) {
      message("  Target genes NOT found on this platform: ", paste(missing_target_genes, collapse = ", "))
    }
    
    gene_expr_df <- as.data.frame(t(expr_filt[found_genes, , drop = FALSE]))
    gene_expr_df$sample_id <- colnames(expr_filt)
    gene_expr_df$condition  <- grp_filt
    gene_expr_df$dataset    <- gse_id
    
    # NEW: within-dataset z-scored copies of the same genes, suffixed "_z".
    # Only used for the Spearman-invariance check in Section 5.
    gene_expr_z <- zscore_within_dataset(expr_filt[found_genes, , drop = FALSE])
    gene_expr_z_df <- as.data.frame(t(gene_expr_z))
    colnames(gene_expr_z_df) <- paste0(colnames(gene_expr_z_df), "_z")
    gene_expr_z_df$sample_id <- colnames(expr_filt)
    gene_expr_df <- dplyr::left_join(gene_expr_df, gene_expr_z_df, by = "sample_id")
    
    gene_expr_all[[gse_id]] <- gene_expr_df
    
    message("  Computing Sepsis vs Control cell-type differences...")
    available_cells <- intersect(CELLS_OF_INTEREST, colnames(xcell_df))
    
    for (cell in available_cells) {
      sep_scores  <- xcell_df[[cell]][xcell_df$condition == "Sepsis"]
      ctrl_scores <- xcell_df[[cell]][xcell_df$condition == "Control"]
      
      if (length(sep_scores) < 3 || length(ctrl_scores) < 3) next
      
      wt <- wilcox.test(sep_scores, ctrl_scores, exact = FALSE)
      
      med_sep  <- median(sep_scores,  na.rm = TRUE)
      med_ctrl <- median(ctrl_scores, na.rm = TRUE)
      mean_sep  <- mean(sep_scores,  na.rm = TRUE)
      mean_ctrl <- mean(ctrl_scores, na.rm = TRUE)
      
      # Standardized (within-dataset z-scored) delta, for cross-platform-
      # comparable MAGNITUDE only -- never used for the test itself.
      z_all   <- scale(c(sep_scores, ctrl_scores))
      z_sep   <- z_all[seq_along(sep_scores)]
      z_ctrl  <- z_all[(length(sep_scores) + 1):length(z_all)]
      delta_z <- median(z_sep, na.rm = TRUE) - median(z_ctrl, na.rm = TRUE)
      
      results_list[[paste(gse_id, cell, sep = "_")]] <- data.frame(
        dataset        = gse_id,
        cell_type      = cell,
        n_sepsis       = length(sep_scores),
        n_control      = length(ctrl_scores),
        median_sepsis  = round(med_sep, 4),
        median_control = round(med_ctrl, 4),
        delta_score    = round(med_sep - med_ctrl, 4),     # PRIMARY (median-based, matches Wilcoxon)
        mean_sepsis    = round(mean_sep, 4),                # secondary, for transparency only
        mean_control   = round(mean_ctrl, 4),
        delta_mean     = round(mean_sep - mean_ctrl, 4),
        delta_zscore   = round(delta_z, 4),                 # standardized magnitude, cross-platform-comparable
        p_value        = wt$p.value,
        stringsAsFactors = FALSE
      )
    }
    
    message("  Done: ", gse_id)
    
  }, error = function(e) {
    message("  ERROR processing ", gse_id, ": ", e$message)
  })
}

# SECTION 3.5: ZERO-INFLATION DIAGNOSTIC (data-derived, not assumed)


message("a", strrep("=", 60))
message("Zero-inflation diagnostic (computed from actual xCell output)...")

compute_zero_inflation <- function(xcell_scores_all, cells, tol = 1e-8) {
  purrr::map_dfr(names(xcell_scores_all), function(gse_id) {
    df <- xcell_scores_all[[gse_id]]
    purrr::map_dfr(intersect(cells, colnames(df)), function(cell) {
      sep  <- df[[cell]][df$condition == "Sepsis"]
      ctrl <- df[[cell]][df$condition == "Control"]
      data.frame(
        dataset           = gse_id,
        cell_type         = cell,
        n_sepsis          = length(sep),
        n_zero_sepsis     = sum(abs(sep) < tol, na.rm = TRUE),
        pct_zero_sepsis   = round(100 * mean(abs(sep) < tol, na.rm = TRUE), 1),
        n_control         = length(ctrl),
        n_zero_control    = sum(abs(ctrl) < tol, na.rm = TRUE),
        pct_zero_control  = round(100 * mean(abs(ctrl) < tol, na.rm = TRUE), 1),
        pct_zero_overall  = round(100 * mean(abs(c(sep, ctrl)) < tol, na.rm = TRUE), 1),
        stringsAsFactors  = FALSE
      )
    })
  })
}

zero_inflation_df <- compute_zero_inflation(xcell_scores_all, CELLS_OF_INTEREST)

message("ZERO-INFLATION BY DATASET x CELL TYPE:")
print(tibble::as_tibble(zero_inflation_df) %>% dplyr::arrange(dataset, cell_type), n = Inf)
readr::write_csv(zero_inflation_df, "experiment3c_zero_inflation.csv")
message("Saved: experiment3c_zero_inflation.csv")


# compile cell-type results


message("", strrep("=", 60))
message("Compiling cell-type enrichment results...")

if (length(results_list) == 0) {
  stop("No results generated. Check errors above.")
}

celltype_df <- dplyr::bind_rows(results_list) %>%
  dplyr::group_by(dataset) %>%
  dplyr::mutate(p_adjusted = p.adjust(p_value, method = "BH")) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    significant  = p_adjusted < 0.05,
    sig_label    = dplyr::case_when(
      p_adjusted < 0.001 ~ "***",
      p_adjusted < 0.01  ~ "**",
      p_adjusted < 0.05  ~ "*",
      TRUE               ~ "ns"
    ),
    direction = ifelse(delta_score > 0, "Higher in Sepsis", "Lower in Sepsis")
  ) %>%
  dplyr::left_join(
    zero_inflation_df %>% dplyr::select(dataset, cell_type, pct_zero_overall),
    by = c("dataset", "cell_type")
  )

ZERO_INFLATION_THRESHOLD <- 15   # percent; adjust if you want a stricter/looser cutoff

concordance_flags <- sapply(seq_len(nrow(celltype_df)), function(i) {
  row    <- celltype_df[i, ]
  others <- celltype_df %>%
    dplyr::filter(cell_type == row$cell_type, dataset != row$dataset, significant)
  if (nrow(others) == 0) return(NA)
  any(sign(others$delta_score) == sign(row$delta_score))
})
celltype_df$concordant_with_significant_other_platform <- concordance_flags

celltype_df <- celltype_df %>%
  dplyr::mutate(
    power_interpretation = dplyr::case_when(
      dataset == "GSE54514" & !significant & isTRUE(concordant_with_significant_other_platform) &
        pct_zero_overall >= ZERO_INFLATION_THRESHOLD ~
        paste0("Power-limited, not a failed replication: n=", n_sepsis + n_control,
               ", ", pct_zero_overall, "% of this cell type's xCell scores are exactly ",
               "zero in this cohort, and the effect direction matches a significant ",
               "result in another (Affymetrix) platform."),
      dataset == "GSE54514" & !significant & isTRUE(concordant_with_significant_other_platform) ~
        paste0("Direction-concordant with a significant result in another platform; ",
               "small n (", n_sepsis + n_control, ") likely limits power here."),
      dataset == "GSE54514" & !significant & !isTRUE(concordant_with_significant_other_platform) ~
        "No concordant significant result in another platform for this cell type -- treated as a genuine null, not a power issue.",
      TRUE ~ NA_character_
    )
  )

message("CELL-TYPE ENRICHMENT RESULTS (median = primary; mean/z-score = secondary):")
print(tibble::as_tibble(celltype_df) %>%
        dplyr::select(dataset, cell_type, n_sepsis, n_control,
                      delta_score, delta_mean, delta_zscore, pct_zero_overall,
                      p_value, p_adjusted, sig_label, direction) %>%
        dplyr::arrange(cell_type, dataset), n = Inf)

message("GSE54514 'ns' RESULTS -- REINTERPRETED (power limitation vs. genuine null):")
print(tibble::as_tibble(celltype_df) %>%
        dplyr::filter(dataset == "GSE54514", !significant) %>%
        dplyr::select(cell_type, n_sepsis, n_control, pct_zero_overall,
                      delta_score, power_interpretation),
      n = Inf, width = Inf)

readr::write_csv(celltype_df, "experiment3c_celltype_results.csv")
message("Saved: experiment3c_celltype_results.csv")

# gene vs cell-type


message("Computing gene-celltype correlations...")

corr_list   <- list()
corr_z_list <- list()   # NEW: parallel z-scored version, for the invariance check

for (gse_id in names(xcell_scores_all)) {
  xcell_df    <- xcell_scores_all[[gse_id]]
  gene_df     <- gene_expr_all[[gse_id]]
  
  merged <- dplyr::inner_join(
    gene_df  %>% dplyr::select(-condition, -dataset),
    xcell_df %>% dplyr::select(-condition, -dataset),
    by = "sample_id"
  )
  
  available_cells <- intersect(CELLS_OF_INTEREST, colnames(xcell_df))
  available_genes <- intersect(GENES, colnames(gene_df))
  
  for (gene in available_genes) {
    gene_z <- paste0(gene, "_z")
    for (cell in available_cells) {
      if (!gene %in% colnames(merged) || !cell %in% colnames(merged)) next
      x <- merged[[gene]]
      y <- merged[[cell]]
      if (sum(!is.na(x) & !is.na(y)) < 5) next
      
      ct <- cor.test(x, y, method = "spearman", exact = FALSE)
      corr_list[[paste(gse_id, gene, cell, sep = "_")]] <- data.frame(
        dataset   = gse_id,
        gene      = gene,
        cell_type = cell,
        rho       = round(unname(ct$estimate), 3),
        p_value   = ct$p.value,
        stringsAsFactors = FALSE
      )
      
      # NEW: same correlation on z-scored expression -- should match
      # the raw-value rho exactly (Spearman is rank-based), demonstrating
      # the result isn't an artifact of platform intensity scale.
      if (gene_z %in% colnames(merged)) {
        xz <- merged[[gene_z]]
        if (sum(!is.na(xz) & !is.na(y)) >= 5) {
          ct_z <- cor.test(xz, y, method = "spearman", exact = FALSE)
          corr_z_list[[paste(gse_id, gene, cell, sep = "_")]] <- data.frame(
            dataset = gse_id, gene = gene, cell_type = cell,
            rho_z = round(unname(ct_z$estimate), 3),
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
}

corr_df <- dplyr::bind_rows(corr_list) %>%
  dplyr::mutate(p_adjusted = p.adjust(p_value, method = "BH"),
                sig_label  = dplyr::case_when(
                  p_adjusted < 0.001 ~ "***",
                  p_adjusted < 0.01  ~ "**",
                  p_adjusted < 0.05  ~ "*",
                  TRUE               ~ "ns"
                ))

corr_z_df <- dplyr::bind_rows(corr_z_list)
if (nrow(corr_z_df) > 0) {
  check <- dplyr::inner_join(corr_df, corr_z_df, by = c("dataset", "gene", "cell_type")) %>%
    dplyr::mutate(abs_diff = abs(rho - rho_z))
  max_diff <- max(check$abs_diff, na.rm = TRUE)
  message("[PLATFORM BATCH-EFFECT CHECK] Max |raw rho - z-scored rho| across all ",
          nrow(check), " gene-celltype pairs: ", round(max_diff, 6))
  if (max_diff < 1e-6) {
    message("  Confirmed: Spearman rho is IDENTICAL under within-dataset z-score",
            " standardization (as expected, since rho depends only on rank order).",
            " This demonstrates these correlations are not driven by differences",
            " in absolute intensity scale between the Illumina (GSE54514) and",
            " Affymetrix (GSE13904/GSE28750) platforms.")
  } else {
    message("  NOTE: non-zero difference detected -- investigate before citing",
            " rank-invariance in the rebuttal (should be ~0 for Spearman).")
  }
}

message("GENE-CELL TYPE CORRELATIONS (significant only):")
print(tibble::as_tibble(corr_df) %>%
        dplyr::filter(p_adjusted < 0.05) %>%
        dplyr::arrange(gene, cell_type, dataset), n = Inf)
readr::write_csv(corr_df, "experiment3c_correlations.csv")
message("Saved: experiment3c_correlations.csv")
message("Working directory: ", getwd())

# figures 

message("Generating figures...")

# ---- Figure 1: cell-type delta heatmap (median-based, primary) -----------
delta_wide <- celltype_df %>%
  dplyr::select(dataset, cell_type, delta_score) %>%
  tidyr::pivot_wider(names_from = dataset, values_from = delta_score) %>%
  tibble::column_to_rownames("cell_type")

sig_wide <- celltype_df %>%
  dplyr::select(dataset, cell_type, sig_label) %>%
  tidyr::pivot_wider(names_from = dataset, values_from = sig_label) %>%
  tibble::column_to_rownames("cell_type")

# FIX #5: force identical row/column order before pheatmap sees them, so the
# significance-star matrix can never silently misalign with the score matrix.
common_rows <- intersect(rownames(delta_wide), rownames(sig_wide))
common_cols <- intersect(colnames(delta_wide), colnames(sig_wide))
delta_wide  <- delta_wide[common_rows, common_cols, drop = FALSE]
sig_wide    <- sig_wide[common_rows, common_cols, drop = FALSE]
sig_wide[is.na(sig_wide)] <- ""   # pheatmap display_numbers can't take NA

if (nrow(delta_wide) > 0) {
  max_val <- max(abs(delta_wide), na.rm = TRUE)
  png("figure1_celltype_heatmap.png", width = 900, height = 700, res = 120)
  pheatmap(as.matrix(delta_wide),
           color           = colorRampPalette(c("#2166AC", "white", "#D6604D"))(100),
           breaks          = seq(-max_val, max_val, length.out = 101),
           cluster_rows    = TRUE,
           cluster_cols    = FALSE,
           display_numbers = as.matrix(sig_wide),
           fontsize        = 11,
           fontsize_number = 10,
           main            = "Cell-Type Enrichment: Sepsis vs Control (median delta)(*** p.adj<0.001, ** <0.01, * <0.05; platforms never pooled)",
           na_col          = "grey90",
           border_color    = "white")
  dev.off()
  message("  Saved: figure1_celltype_heatmap.png")
}

# ---- Figure 1b: SUPPLEMENTARY standardized (z-scored) delta heatmap ------
# For readers/reviewers who want to compare RELATIVE magnitude across the
# two platforms; the significance stars are IDENTICAL to Figure 1 (same
# underlying test) -- only the color scale changes to SD units.
zdelta_wide <- celltype_df %>%
  dplyr::select(dataset, cell_type, delta_zscore) %>%
  tidyr::pivot_wider(names_from = dataset, values_from = delta_zscore) %>%
  tibble::column_to_rownames("cell_type")
zdelta_wide <- zdelta_wide[common_rows, intersect(colnames(zdelta_wide), common_cols), drop = FALSE]

if (nrow(zdelta_wide) > 0) {
  sig_wide_z <- sig_wide[rownames(zdelta_wide), colnames(zdelta_wide), drop = FALSE]
  max_val_z  <- max(abs(zdelta_wide), na.rm = TRUE)
  png("figure1b_celltype_heatmap_standardized.png", width = 900, height = 700, res = 120)
  pheatmap(as.matrix(zdelta_wide),
           color           = colorRampPalette(c("#2166AC", "white", "#D6604D"))(100),
           breaks          = seq(-max_val_z, max_val_z, length.out = 101),
           cluster_rows    = TRUE,
           cluster_cols    = FALSE,
           display_numbers = as.matrix(sig_wide_z),
           fontsize        = 11,
           fontsize_number = 10,
           main            = "Cell-Type Enrichment: Standardized delta (SD units, within-platform z-score)Supplementary -- for relative cross-platform magnitude comparison only",
           na_col          = "grey90",
           border_color    = "white")
  dev.off()
  message("  Saved: figure1b_celltype_heatmap_standardized.png")
}

# ---- Figure 2: NLRP3 vs Monocytes, annotated with the EXACT stored rho ---
nlrp3_mono_list <- list()
for (gse_id in names(xcell_scores_all)) {
  xcell_df <- xcell_scores_all[[gse_id]]
  gene_df  <- gene_expr_all[[gse_id]]
  if (!"NLRP3" %in% colnames(gene_df)) next
  if (!"Monocytes" %in% colnames(xcell_df)) next
  
  merged <- dplyr::inner_join(
    gene_df  %>% dplyr::select(sample_id, NLRP3, condition),
    xcell_df %>% dplyr::select(sample_id, Monocytes),
    by = "sample_id"
  ) %>% dplyr::filter(condition %in% c("Sepsis", "Control"))
  
  merged$dataset <- gse_id
  nlrp3_mono_list[[gse_id]] <- merged
}

if (length(nlrp3_mono_list) > 0) {
  nlrp3_mono <- dplyr::bind_rows(nlrp3_mono_list)
  
  # Annotation positions computed from the actual data range per facet,
  # labels pulled directly from corr_df (BH-adjusted) -- not recomputed.
  ann_nlrp3 <- corr_df %>%
    dplyr::filter(gene == "NLRP3", cell_type == "Monocytes") %>%
    dplyr::inner_join(
      nlrp3_mono %>%
        dplyr::group_by(dataset) %>%
        dplyr::summarise(
          x = min(Monocytes, na.rm = TRUE) + 0.05 * diff(range(Monocytes, na.rm = TRUE)),
          y = max(NLRP3, na.rm = TRUE) - 0.05 * diff(range(NLRP3, na.rm = TRUE)),
          .groups = "drop"
        ),
      by = "dataset"
    ) %>%
    dplyr::mutate(label = paste0("rho = ", rho, "  ", sig_label, " (p.adj = ", signif(p_adjusted, 2), ")"))
  
  p2 <- ggplot(nlrp3_mono, aes(x = Monocytes, y = NLRP3, color = condition)) +
    geom_point(alpha = 0.6, size = 1.8) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.8, color = "grey40") +
    geom_text(data = ann_nlrp3, aes(x = x, y = y, label = label),
              inherit.aes = FALSE, hjust = 0, size = 3.1, fontface = "bold") +
    facet_wrap(~ dataset, scales = "free") +
    scale_color_manual(values = c("Sepsis" = "#e74c3c", "Control" = "#3498db")) +
    labs(title    = "NLRP3 Expression vs Monocyte Enrichment Score",
         subtitle = "Annotated rho/p.adj = exact Spearman values from experiment3c_correlations.csv (BH-adjusted); trend line is an illustrative OLS fit only",
         x = "xCell Monocyte Score", y = "NLRP3 Expression") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))
  ggsave("figure2_NLRP3_vs_monocytes.png", p2, width = 10, height = 4, dpi = 150)
  message("  Saved: figure2_NLRP3_vs_monocytes.png")
}

gclm_nk_list <- list()
for (gse_id in names(xcell_scores_all)) {
  xcell_df <- xcell_scores_all[[gse_id]]
  gene_df  <- gene_expr_all[[gse_id]]
  if (!"GCLM" %in% colnames(gene_df)) next
  if (!"NK cells" %in% colnames(xcell_df)) next
  
  merged <- dplyr::inner_join(
    gene_df  %>% dplyr::select(sample_id, GCLM, condition),
    xcell_df %>% dplyr::select(sample_id, `NK cells`),
    by = "sample_id"
  ) %>% dplyr::filter(condition %in% c("Sepsis", "Control"))
  
  merged$dataset <- gse_id
  gclm_nk_list[[gse_id]] <- merged
}

if (length(gclm_nk_list) > 0) {
  gclm_nk <- dplyr::bind_rows(gclm_nk_list)
  
  ann_gclm <- corr_df %>%
    dplyr::filter(gene == "GCLM", cell_type == "NK cells") %>%
    dplyr::inner_join(
      gclm_nk %>%
        dplyr::group_by(dataset) %>%
        dplyr::summarise(
          x = min(`NK cells`, na.rm = TRUE) + 0.05 * diff(range(`NK cells`, na.rm = TRUE)),
          y = max(GCLM, na.rm = TRUE) - 0.05 * diff(range(GCLM, na.rm = TRUE)),
          .groups = "drop"
        ),
      by = "dataset"
    ) %>%
    dplyr::mutate(label = paste0("rho = ", rho, "  ", sig_label, " (p.adj = ", signif(p_adjusted, 2), ")"))
  
  p3 <- ggplot(gclm_nk, aes(x = `NK cells`, y = GCLM, color = condition)) +
    geom_point(alpha = 0.6, size = 1.8) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.8, color = "grey40") +
    geom_text(data = ann_gclm, aes(x = x, y = y, label = label),
              inherit.aes = FALSE, hjust = 0, size = 3.1, fontface = "bold") +
    facet_wrap(~ dataset, scales = "free") +
    scale_color_manual(values = c("Sepsis" = "#e74c3c", "Control" = "#3498db")) +
    labs(title    = "GCLM Expression vs NK Cell Enrichment Score",
         subtitle = "Annotated rho/p.adj = exact Spearman values from experiment3c_correlations.csv (BH-adjusted); trend line is an illustrative OLS fit only",
         x = "xCell NK Cell Score", y = "GCLM Expression") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))
  ggsave("figure3_GCLM_vs_NKcells.png", p3, width = 10, height = 4, dpi = 150)
  message("  Saved: figure3_GCLM_vs_NKcells.png")
}

plot_cells <- c("Monocytes", "NK cells", "Neutrophils", "CD8+ T-cells")

all_xcell <- dplyr::bind_rows(xcell_scores_all)
available_plot_cells <- intersect(plot_cells, colnames(all_xcell))

if (length(available_plot_cells) > 0) {

  y_pos <- all_xcell %>%
    dplyr::filter(condition %in% c("Sepsis", "Control")) %>%
    dplyr::select(dataset, dplyr::all_of(available_plot_cells)) %>%
    tidyr::pivot_longer(-dataset, names_to = "cell_type", values_to = "score") %>%
    dplyr::group_by(dataset, cell_type) %>%
    dplyr::summarise(y.position = max(score, na.rm = TRUE) * 1.08, .groups = "drop")
  
  stats_df <- celltype_df %>%
    dplyr::filter(cell_type %in% available_plot_cells) %>%
    dplyr::transmute(dataset, cell_type,
                     group1 = "Control", group2 = "Sepsis",
                     p.adj = p_adjusted, sig_label) %>%
    dplyr::left_join(y_pos, by = c("dataset", "cell_type"))
  
  plot_list <- list()
  for (cell in available_plot_cells) {
    df_plot <- all_xcell %>%
      dplyr::filter(condition %in% c("Sepsis", "Control")) %>%
      dplyr::select(dplyr::all_of(c("condition", "dataset", cell))) %>%
      dplyr::rename(score = dplyr::all_of(cell))
    
    stats_cell <- stats_df %>% dplyr::filter(cell_type == cell)
    
    p <- ggplot(df_plot, aes(x = condition, y = score, fill = condition)) +
      geom_boxplot(outlier.shape = 21, outlier.size = 1.2, alpha = 0.8) +
      geom_jitter(width = 0.15, size = 0.6, alpha = 0.3) +
      scale_fill_manual(values = c("Sepsis" = "#e74c3c", "Control" = "#3498db")) +
      ggpubr::stat_pvalue_manual(
        data = stats_cell, label = "sig_label",
        xmin = "group1", xmax = "group2", y.position = "y.position",
        tip.length = 0.01
      ) +
      facet_wrap(~ dataset, scales = "free_y", nrow = 1) +
      labs(title = cell, x = NULL, y = "xCell Score") +
      theme_bw(base_size = 9) +
      theme(legend.position = "none",
            plot.title = element_text(face = "bold", size = 10))
    plot_list[[cell]] <- p
  }
  
  combined <- ggarrange(plotlist = plot_list,
                        ncol = 1, nrow = length(plot_list))
  ggsave("figure4_celltype_boxplots.png", combined,
         width = 10, height = 4 * length(plot_list), dpi = 150, limitsize = FALSE)
  message("  Saved: figure4_celltype_boxplots.png (significance stars = BH-adjusted, matches results table exactly)")
}

message("", strrep("=", 60))
message("EXPERIMENT 3c COMPLETE (v3)")
message(strrep("=", 60))
message("Output files in ", output_dir, ":")
message("  experiment3c_zero_inflation.csv          - % exactly-zero xCell scores per dataset x cell type (data-derived)")
message("  experiment3c_celltype_results.csv       - Sepsis vs Control per cell type (median=primary, mean/z-score=secondary; includes power_interpretation for GSE54514 'ns' rows)")
message("  experiment3c_correlations.csv           - Gene ~ cell type Spearman rho (BH-adjusted)")
message("  figure1_celltype_heatmap.png            - Delta score heatmap (median-based, primary)")
message("  figure1b_celltype_heatmap_standardized.png - Same data, z-scored for cross-platform magnitude comparison")
message("  figure2_NLRP3_vs_monocytes.png          - NLRP3 vs monocyte scatter (annotated with exact stored rho)")
message("  figure3_GCLM_vs_NKcells.png             - GCLM vs NK cell scatter (annotated with exact stored rho)")
message("  figure4_celltype_boxplots.png           - Cell type boxplots (stars = BH-adjusted, matches table)")
message("All comparisons and correlations above were computed strictly WITHIN",
        " one platform at a time -- GSE54514 (Illumina) and GSE13904/GSE28750",
        " (Affymetrix) values are never pooled or merged in any single test.")