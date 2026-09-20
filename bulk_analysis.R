# Datasets: GSE54514, GSE13904, GSE28750

library(GEOquery)
library(limma)
library(Biobase)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(pheatmap)
library(RColorBrewer)
library(writexl)

GENES <- c("GCLC", "GCLM", "NFE2L2", "PINK1", "PRKN", "NLRP3", "GSDMD", "SLC25A39", "SLC25A40", "AFG3L2")

datasets <- c("GSE54514", "GSE13904", "GSE28750")

dir.create("experiment3b_output_final", showWarnings = FALSE)
setwd("experiment3b_output_final")

# Helper functions

map_probes_to_genes <- function(eset) {
  fdata <- fData(eset)
  
  symbol_cols <- c("Gene Symbol", "Gene symbol", "GENE_SYMBOL", "gene_assignment",
                   "Symbol", "GeneSymbol", "Gene_Symbol", "SYMBOL",
                   "gene_symbol", "GeneName")
  
  sym_col <- intersect(symbol_cols, colnames(fdata))[1]
  
  if (is.na(sym_col)) {
    message("  WARNING: Could not find gene symbol column.")
    return(NULL)
  }
  
  message("  Using column '", sym_col, "' for gene symbols")
  fdata$gene_symbol_clean <- fdata[[sym_col]]
  fdata$gene_symbol_clean <- gsub(" ///.*", "", fdata$gene_symbol_clean)
  fdata$gene_symbol_clean <- trimws(fdata$gene_symbol_clean)
  
  expr_mat <- exprs(eset)
  expr_mat <- as.data.frame(expr_mat)
  expr_mat$gene_symbol <- fdata$gene_symbol_clean[match(rownames(expr_mat), rownames(fdata))]
  
  expr_genes <- expr_mat %>%
    filter(gene_symbol %in% GENES) %>%
    group_by(gene_symbol) %>%
    summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
    column_to_rownames("gene_symbol")
  
  return(expr_genes)
}

extract_groups <- function(eset, dataset_id) {
  pdata <- pData(eset)
  groups <- list()
  raw_title <- as.character(pdata$title)
  
  groups$condition_std <- case_when(
    grepl("control|healthy|normal", raw_title, ignore.case = TRUE) ~ "Control",
    grepl("sepsis|septic|shock", raw_title, ignore.case = TRUE)    ~ "Sepsis",
    TRUE ~ "Other"
  )
  
  # Search all characteristics columns for Day and Patient ID mapping
  char_cols <- grep("characteristics", colnames(pdata), value = TRUE)
  search_text <- raw_title
  if (length(char_cols) > 0) {
    for (cc in char_cols) {
      search_text <- paste(search_text, as.character(pdata[[cc]]), sep = " | ")
    }
  }
  
  day_match <- stringr::str_match(search_text, "(?i)(?:day|time\\s*point|time)[_\\s:=]*(\\d+)")[, 2]
  groups$day <- as.integer(day_match)
  
  # Force Healthy Controls to Day 1
  groups$day[groups$condition_std == "Control" & is.na(groups$day)] <- 1
  
  pid_match <- stringr::str_match(search_text, "(?i)(?:patient|subject|id)[_\\s:=]*(\\w+)")[, 2]
  groups$patient_id <- pid_match
  
  # fallback - GSM accession number for faulty or no ID's
  missing_idx <- which(is.na(groups$patient_id))
  if (length(missing_idx) > 0) {
    groups$patient_id[missing_idx] <- sampleNames(eset)[missing_idx]
  }
  
  message("  Condition distribution (all samples, pre-filtering) for ", dataset_id, ":")
  print(table(groups$condition_std))
  return(groups)
}

filter_to_day1 <- function(expr_genes, groups, dataset_id) {
  n_before <- ncol(expr_genes)
  
  df <- data.frame(
    idx        = seq_len(n_before),
    condition  = groups$condition_std,
    day        = groups$day,
    patient_id = groups$patient_id,
    stringsAsFactors = FALSE
  )
  
  df <- df %>% dplyr::filter(!is.na(patient_id), patient_id != "")
  
  df_ctrl_kept <- df %>%
    dplyr::filter(condition == "Control") %>%
    dplyr::filter(day == 1 | is.na(day)) %>% 
    dplyr::group_by(patient_id) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
  
  # Dedup (Strictly Day 1, take first available per ID)
  df_sepsis_kept <- df %>%
    dplyr::filter(condition == "Sepsis") %>%
    dplyr::filter(day == 1 | is.na(day)) %>% 
    dplyr::group_by(patient_id) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
  
  df_keep <- dplyr::bind_rows(df_ctrl_kept, df_sepsis_kept) %>% dplyr::arrange(idx)
  keep_idx <- seq_len(n_before) %in% df_keep$idx
  
  expr_filtered <- expr_genes[, keep_idx, drop = FALSE]
  groups_filtered <- lapply(groups, function(g) {
    if (length(g) == n_before) g[keep_idx] else g
  })
  
  message("  Day-1 deduplication filter applied to ", dataset_id, ":")
  message("    Sepsis:  -> ", nrow(df_sepsis_kept), " unique day-1 patients")
  message("    Control: -> ", nrow(df_ctrl_kept), " unique cross-sectional controls")
  
  list(expr = expr_filtered, groups = groups_filtered,
       n_sepsis_samples_used = nrow(df_sepsis_kept), n_sepsis_patients = nrow(df_sepsis_kept))
}

run_mwu <- function(expr_vec, group_vec, gene, dataset) {
  df <- data.frame(expr = expr_vec, group = group_vec) %>%
    filter(group %in% c("Sepsis", "Control"), !is.na(expr))
  
  if (nrow(df) < 6 || length(unique(df$group)) < 2) {
    return(data.frame(dataset = dataset, gene = gene,
                      n_sepsis = NA, n_control = NA,
                      median_sepsis = NA, median_control = NA,
                      fold_change = NA, log2FC = NA,
                      W_statistic = NA, p_value = NA,
                      p_adjusted = NA, direction = NA,
                      note = "Insufficient samples"))
  }
  
  sep  <- df$expr[df$group == "Sepsis"]
  ctrl <- df$expr[df$group == "Control"]
  
  test <- wilcox.test(sep, ctrl, exact = FALSE)
  
  med_sep  <- median(sep,  na.rm = TRUE)
  med_ctrl <- median(ctrl, na.rm = TRUE)
  
  data_range <- range(df$expr, na.rm = TRUE)
  is_log <- data_range[2] < 25
  
  if (is_log) {
    log2FC     <- med_sep - med_ctrl
    fold_change <- 2^log2FC
  } else {
    fold_change <- med_sep / med_ctrl
    log2FC      <- log2(fold_change)
  }
  
  data.frame(
    dataset      = dataset,
    gene         = gene,
    n_sepsis     = length(sep),
    n_control    = length(ctrl),
    median_sepsis  = round(med_sep,  3),
    median_control = round(med_ctrl, 3),
    fold_change  = round(fold_change, 3),
    log2FC       = round(log2FC, 3),
    W_statistic  = test$statistic,
    p_value      = test$p.value,
    p_adjusted   = NA,
    direction    = ifelse(log2FC > 0, "UP in Sepsis", "DOWN in Sepsis"),
    note         = ""
  )
}

# Download and processing

results_list  <- list()
expr_data_all <- list()
sample_size_log <- list()

for (gse_id in datasets) {
  message("", strrep("=", 60))
  message("Processing ", gse_id, "...")
  message(strrep("=", 60))
  
  tryCatch({
    gse  <- getGEO(gse_id, GSEMatrix = TRUE, AnnotGPL = TRUE, destdir = ".")
    eset <- if (is.list(gse)) gse[[1]] else gse
    
    message("  Dimensions: ", nrow(eset), " probes x ", ncol(eset), " samples")
    
    expr_genes <- map_probes_to_genes(eset)
    
    if (is.null(expr_genes) || nrow(expr_genes) == 0) {
      message("  ERROR: Could not map genes. Skipping ", gse_id)
      next
    }
    
    found_genes   <- intersect(GENES, rownames(expr_genes))
    missing_genes <- setdiff(GENES, rownames(expr_genes))
    message("  Genes found:     ", paste(found_genes,  collapse = ", "))
    
    groups <- extract_groups(eset, gse_id)
    
    filt <- filter_to_day1(expr_genes, groups, gse_id)
    expr_genes <- filt$expr
    groups     <- filt$groups
    
    ssu <- if (is.null(filt$n_sepsis_samples_used)) sum(groups$condition_std == "Sepsis") else filt$n_sepsis_samples_used
    snp <- if (is.null(filt$n_sepsis_patients)) NA else filt$n_sepsis_patients
    sample_size_log[[gse_id]] <- list(
      sepsis_samples_used = ssu,
      sepsis_patients      = snp,
      control_samples      = sum(groups$condition_std == "Control")
    )
    
    expr_data_all[[gse_id]] <- list(expr = expr_genes, groups = groups)
    
    for (gene in found_genes) {
      expr_vec  <- as.numeric(expr_genes[gene, ])
      group_vec <- groups$condition_std
      result    <- run_mwu(expr_vec, group_vec, gene, gse_id)
      results_list[[paste(gse_id, gene, sep = "_")]] <- result
    }
    
    message("  Done: ", gse_id)
    
  }, error = function(e) {
    message("  ERROR processing ", gse_id, ": ", e$message)
  })
}

#Compiling + FDR

message("", strrep("=", 60))
message("Compiling results...")

results_df <- bind_rows(results_list) %>%
  group_by(dataset) %>%
  mutate(p_adjusted = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  mutate(
    significant = p_adjusted < 0.05,
    sig_label   = case_when(
      p_adjusted < 0.001 ~ "***",
      p_adjusted < 0.01  ~ "**",
      p_adjusted < 0.05  ~ "*",
      TRUE               ~ "ns"
    )
  )

message("RESULTS TABLE:")
print(tibble::as_tibble(results_df) %>%
        select(dataset, gene, n_sepsis, n_control, log2FC,
               p_value, p_adjusted, sig_label, direction) %>%
        arrange(gene, dataset), n = Inf)

write_csv(results_df, "experiment3b_results.csv")
message("Saved: experiment3b_results.csv")

if (length(sample_size_log) > 0) {
  ss_df <- bind_rows(lapply(names(sample_size_log), function(id) {
    x <- sample_size_log[[id]]
    data.frame(dataset = id,
               sepsis_samples_used = x$sepsis_samples_used,
               sepsis_patients     = x$sepsis_patients,
               control_samples     = x$control_samples)
  }))
  write_csv(ss_df, "experiment3b_sample_sizes.csv")
  message("Saved: experiment3b_sample_sizes.csv")
}

consistency <- results_df %>%
  filter(!is.na(log2FC)) %>%
  group_by(gene) %>%
  summarise(
    n_datasets_tested      = n(),
    n_datasets_significant = sum(significant, na.rm = TRUE),
    n_upregulated          = sum(log2FC > 0, na.rm = TRUE),
    n_downregulated        = sum(log2FC < 0, na.rm = TRUE),
    mean_log2FC            = round(mean(log2FC, na.rm = TRUE), 3),
    sd_log2FC              = round(sd(log2FC,   na.rm = TRUE), 3),
    consistent_direction   = n_upregulated == n() | n_downregulated == n(),
    consensus_direction    = ifelse(n_upregulated > n_downregulated,
                                    "Consistently UP", "Consistently DOWN"),
    .groups = "drop"
  ) %>%
  arrange(gene)

write_csv(consistency, "experiment3b_consistency.csv")

#Figures (supplementary materials)

message("Generating figures...")

fc_matrix <- results_df %>%
  filter(!is.na(log2FC)) %>%
  select(gene, dataset, log2FC) %>%
  pivot_wider(names_from = dataset, values_from = log2FC) %>%
  column_to_rownames("gene")

gene_order <- intersect(GENES, rownames(fc_matrix))
fc_matrix  <- fc_matrix[gene_order, , drop = FALSE]
max_val    <- max(abs(fc_matrix), na.rm = TRUE)

png("figure1_log2FC_heatmap.png", width = 800, height = 600, res = 120)
pheatmap(fc_matrix,
         color         = colorRampPalette(c("#2166AC","white","#D6604D"))(100),
         breaks        = seq(-max_val, max_val, length.out = 101),
         cluster_rows  = FALSE,
         cluster_cols  = FALSE,
         display_numbers = TRUE,
         number_format = "%.2f",
         fontsize      = 12,
         main          = "Log2 Fold Change: Sepsis vs Control",
         na_col        = "grey90",
         border_color  = "white")
dev.off()
message("  Saved: figure1_log2FC_heatmap.png")

forest_data <- results_df %>%
  filter(!is.na(log2FC)) %>%
  mutate(
    label   = paste0(gene, " (", dataset, ")"),
    pathway = case_when(
      gene %in% c("GCLC","GCLM") ~ "GSH Synthesis",
      gene == "NFE2L2"            ~ "Nrf2",
      gene %in% c("PINK1","PRKN") ~ "Mitophagy",
      gene %in% c("NLRP3","GSDMD") ~ "Inflammasome/Pyroptosis",
      gene %in% c("SLC25A39", "SLC25A40", "AFG3L2") ~ "Mitochondrial Transport & Regulation",
      TRUE ~ "Other"
    )
  )

p2 <- ggplot(forest_data,
             aes(x = log2FC, y = reorder(label, log2FC),
                 color = pathway, shape = significant)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 3.5, stroke = 1.2) +
  scale_shape_manual(values = c("FALSE" = 1, "TRUE" = 16),
                     labels = c("FALSE" = "p.adj >= 0.05", "TRUE" = "p.adj < 0.05"),
                     name   = "Significance") +
  scale_color_manual(values = c(
    "GSH Synthesis"                        = "#1a7abf",
    "Nrf2"                                 = "#3aacda",
    "Mitophagy"                            = "#9b59b6",
    "Inflammasome/Pyroptosis"              = "#e74c3c",
    "Mitochondrial Transport & Regulation" = "#f39c12",
    "Other"                                = "grey40"), name = "Pathway") +
  labs(title    = "Gene Expression Changes: Sepsis vs Healthy Controls",
       subtitle = "Filled circle = p.adj < 0.05; strict deduplication to 1 sample/patient",
       x = "Log2 Fold Change", y = NULL) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

ggsave("figure2_forest_plot.png", p2, width = 9, height = 7, dpi = 150)
message("  Saved: figure2_forest_plot.png")