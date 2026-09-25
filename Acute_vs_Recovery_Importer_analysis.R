
library(GEOquery)
library(Biobase)
library(tidyverse)
library(lmerTest)


message("Downloading GSE46955 to temporary directory...")
temp_dir <- tempdir()
gse_data <- getGEO("GSE46955", GSEMatrix = TRUE, AnnotGPL = TRUE, destdir = temp_dir)
eset <- gse_data[[1]]


fdata <- fData(eset)
expr_mat <- as.data.frame(exprs(eset))


gene_sym_col <- "Gene symbol" 
expr_mat$gene_symbol <- fdata[[gene_sym_col]][match(rownames(expr_mat), rownames(fdata))]


target_expr <- expr_mat %>%
  dplyr::filter(grepl("AFG3L2|SLC25A39", gene_symbol)) %>%
  dplyr::group_by(gene_symbol) %>%
  dplyr::summarise(dplyr::across(dplyr::where(is.numeric), \(x) max(x, na.rm = TRUE)), .groups = "drop") %>%
  tibble::column_to_rownames("gene_symbol") %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("sample_id")

# 3. Extract Clinical Metadata (Corrected Regex & Filtering)
pdata <- pData(eset)

meta <- data.frame(
  sample_id = rownames(pdata),
  title = pdata$title,
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(
    # Extract just the patient number from the end of the title
    patient_id = stringr::str_extract(title, "\\d+$"),
    
    # Map timepoints based on exact title text
    timepoint = dplyr::case_when(
      grepl("sepsis", title, ignore.case = TRUE) ~ "Acute",
      grepl("recovery", title, ignore.case = TRUE) ~ "Recovery",
      grepl("healthy", title, ignore.case = TRUE) ~ "Healthy",
      TRUE ~ "Unknown"
    ),
    
    # Identify whether the sample was basal or LPS-stimulated
    stimulation = ifelse(grepl("lps", title, ignore.case = TRUE), "LPS", "Basal")
  ) %>%
  # Keep only the basal samples (reduces the 44 rows back down to the 22 independent samples)
  dplyr::filter(stimulation == "Basal")

long_df <- dplyr::inner_join(meta, target_expr, by = "sample_id")

message("Data merged. Final row count: ", nrow(long_df), " (Expected: 22)")

Run the Linear Mixed-Effects Model
message("--- Running Linear Mixed-Effects Model ---")
lmm_fit <- lmerTest::lmer(SLC25A39 ~ AFG3L2 + (1 | patient_id), data = long_df)
print(summary(lmm_fit))


message("Running Spearman (For Comparison)")
cor_test <- cor.test(long_df$AFG3L2, long_df$SLC25A39, method = "spearman")
print(cor_test)
