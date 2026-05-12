############################################################
# Bulk EWAS: Cancer vs Normal DNA methylation analysis
#
# Author:
# Barbara Karakyriakou
#
# Description:
# limma-based epigenome-wide association study (EWAS)
# comparing tumor and normal breast tissue DNA methylation
# profiles with optional adjustment for estimated cell-type
# proportions.
#
# Associated manuscript:
# "Tumor endothelial cell-specific DNA methylation
# alterations in breast cancer"
#
# Date: 2026
############################################################

library(minfi)
library(limma)
library(dplyr)
library(tibble)
library(EnhancedVolcano)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# -----------------------------
# 1. Load input data
# -----------------------------

# beta_mat:
#   Matrix of DNA methylation beta values
#   Rows = CpG probes
#   Columns = samples
#   Column names must match pheno$Sample_Name

# pheno:
#   Sample metadata table
#   Required columns:
#   Sample_Name, Status, Age
#   For cell-type-adjusted model, also include:
#   Myeloid, Lymphocyte, Endothelial, Epithelial, Stromal

beta_mat <- readRDS("beta_matrix.rds")
pheno <- read.csv("sample_metadata.csv")

# Align phenotype rows to beta matrix columns
pheno <- pheno[match(colnames(beta_mat), pheno$Sample_Name), ]

# Sanity check
stopifnot(all(colnames(beta_mat) == pheno$Sample_Name))

# -----------------------------
# 2. Select EWAS model
# -----------------------------

# Choose one:
#   "age"          = adjusted for age only
#   "age_celltype" = adjusted for age and estimated cell-type proportions

model_type <- "age"
# model_type <- "age_celltype"

pheno$Status <- factor(pheno$Status, levels = c("Normal", "Cancer"))

if (model_type == "age") {
  
  design <- model.matrix(
    ~ 0 + Status + Age,
    data = pheno
  )
  
  colnames(design) <- c("Normal", "Cancer", "Age")
  model_label <- "age_adjusted"
  plot_subtitle <- "Age-adjusted bulk EWAS"
  
} else if (model_type == "age_celltype") {
  
  design <- model.matrix(
    ~ 0 + Status + Age +
      Myeloid + Lymphocyte + Endothelial + Epithelial + Stromal,
    data = pheno
  )
  
  colnames(design) <- c(
    "Normal", "Cancer", "Age",
    "Myeloid", "Lymphocyte", "Endothelial", "Epithelial", "Stromal"
  )
  
  model_label <- "age_celltype_adjusted"
  plot_subtitle <- "Age- and cell-type-adjusted bulk EWAS"
  
} else {
  stop("model_type must be either 'age' or 'age_celltype'")
}

# -----------------------------
# 3. Run limma EWAS
# -----------------------------

fit <- lmFit(beta_mat, design)

contrast_matrix <- makeContrasts(
  Cancer_vs_Normal = Cancer - Normal,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# -----------------------------
# 4. Annotate EWAS results
# -----------------------------

ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

ann450k_sub <- ann450k[
  match(rownames(beta_mat), ann450k$Name),
  c(1:4, 12:19, 24:ncol(ann450k))
]

ewas_results <- topTable(
  fit2,
  coef = "Cancer_vs_Normal",
  num = Inf,
  genelist = ann450k_sub,
  sort.by = "P"
)

ewas_results <- ewas_results %>%
  rownames_to_column("CpG") %>%
  rename(Delta_Beta = logFC) %>%
  mutate(
    absDelta_Beta = abs(Delta_Beta),
    significant = adj.P.Val < 0.05 & absDelta_Beta >= 0.2
  )

# -----------------------------
# 5. Save result tables
# -----------------------------

write.csv(
  ewas_results,
  paste0("EWAS_Cancer_vs_Normal_", model_label, "_all_CpGs.csv"),
  row.names = FALSE
)

ewas_sig <- ewas_results %>%
  filter(significant)

write.csv(
  ewas_sig,
  paste0("EWAS_Cancer_vs_Normal_", model_label, "_significant_CpGs.csv"),
  row.names = FALSE
)

# Summary:
# Delta_Beta > 0 = hypermethylated in cancer
# Delta_Beta < 0 = hypomethylated in cancer

table(ewas_sig$Delta_Beta > 0)
table(ewas_sig$Delta_Beta < 0)

# -----------------------------
# 6. Volcano plot
# -----------------------------

p_value_cutoff <- max(
  ewas_results$P.Value[ewas_results$adj.P.Val < 0.05],
  na.rm = TRUE
)

volcano_plot <- EnhancedVolcano(
  ewas_results,
  x = "Delta_Beta",
  y = "P.Value",
  lab = NA,
  FCcutoff = 0.2,
  pCutoff = p_value_cutoff,
  xlim = c(-0.6, 0.6),
  ylim = c(0, 330),
  title = "Cancer vs Normal",
  subtitle = plot_subtitle,
  caption = paste0("Total = ", nrow(ewas_results), " CpGs"),
  xlab = expression(Delta~Beta),
  ylab = expression(-log[10](P)),
  legendPosition = "none",
  gridlines.major = FALSE,
  gridlines.minor = FALSE,
  pointSize = 0.5,
  col = c("#999999", "#e699ff", "#8080ff", "#ff0066")
)

volcano_plot

ggsave(
  filename = paste0("EWAS_Cancer_vs_Normal_", model_label, "_volcano.pdf"),
  plot = volcano_plot,
  width = 6,
  height = 5
)