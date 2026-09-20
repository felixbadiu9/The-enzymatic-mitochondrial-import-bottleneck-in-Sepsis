
#AFG3L2 ~ SLC25A39 correlation in Acute Sepsis vs Recovery GSE46955

library(GEOquery)
library(limma)
library(dplyr)
library(ggplot2)

gse46955 <- getGEO("GSE46955", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]


expr_mat <- exprs(gse46955)
fdata <- fData(gse46955)

sym_col <- intersect(c("Gene symbol", "Gene Symbol", "Symbol", "SYMBOL"), colnames(fdata))[1]
gene_symbols <- as.character(fdata[[sym_col]])
gene_symbols <- gsub(" ///.*", "", gene_symbols) 

# 4. Collapse Probes 
expr_collapsed <- avereps(expr_mat, ID = gene_symbols)

targets <- c("AFG3L2", "SLC25A39")
df_expr <- as.data.frame(t(expr_collapsed[targets, ]))
df_expr$sample_id <- rownames(df_expr)

# 6. Extract and Stratify Metadata
pdata <- pData(gse46955)
df_meta <- data.frame(
  sample_id = rownames(pdata),
  title = as.character(pdata$title),
  stringsAsFactors = FALSE
)

df_meta <- df_meta %>%
  dplyr::filter(!grepl("LPS", title, ignore.case = TRUE)) %>%
  dplyr::mutate(
    clinical_phase = dplyr::case_when(
      grepl("acute|sepsis_acute|d1", title, ignore.case = TRUE) ~ "Acute Sepsis",
      grepl("recovery|convalescent", title, ignore.case = TRUE) ~ "Recovery",
      grepl("healthy|control", title, ignore.case = TRUE) ~ "Control",
      TRUE ~ "Unknown"
    )
  )


df_final <- inner_join(df_expr, df_meta, by = "sample_id")


run_cor <- function(data_subset, label) {
  if(nrow(data_subset) < 3) return(NULL)
  
  ct <- cor.test(data_subset$AFG3L2, data_subset$SLC25A39, 
                 method = "spearman", exact = FALSE)
  
  data.frame(
    Phase = label,
    n_samples = nrow(data_subset),
    Spearman_rho = round(ct$estimate, 3),
    p_value = ct$p.value
  )
}

results_table <- bind_rows(
  run_cor(df_final, "Aggregate (Composite)"),
  run_cor(df_final %>% filter(clinical_phase == "Acute Sepsis"), "Acute Sepsis"),
  run_cor(df_final %>% filter(clinical_phase == "Recovery"), "Recovery")
) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    Sig = case_when(
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

message("AFG3L2 ~ SLC25A39 co-transcription results")
print(results_table)