############################################################
# CellDMC analysis
#
# Authors:
# Barbara Karakyriakou

# Description:
# Identifies cell-type-specific differentially methylated
# CpGs using CellDMC with HiTIMED-derived cell-type
# proportions and sample-level covariates.
############################################################

library(EpiDISH)
library(dplyr)
library(tibble)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# -----------------------------
# 1. Required input objects
# -----------------------------

# beta_mat:
#   DNA methylation beta matrix
#   Rows = CpGs
#   Columns = samples
#
# pheno:
#   Sample metadata dataframe
#   Required columns:
#   Sample_Name, Status, Age, Sex, Cohort
#
# cell_fractions:
#   HiTIMED deconvolution output
#   Rows = samples
#   Columns = cell-type proportions

beta_mat <- readRDS("processed_beta_matrix.rds")
pheno <- read.csv("sample_metadata.csv")
cell_fractions <- read.csv(
  "HiTIMED_deconvolution_for_CellDMC.csv",
  row.names = 1
)

# -----------------------------
# 2. Align inputs
# -----------------------------

pheno <- pheno[match(colnames(beta_mat), pheno$Sample_Name), ]
cell_fractions <- cell_fractions[
  match(pheno$Sample_Name, rownames(cell_fractions)),
]

stopifnot(all(colnames(beta_mat) == pheno$Sample_Name))
stopifnot(all(rownames(cell_fractions) == pheno$Sample_Name))

pheno$Status <- factor(pheno$Status, levels = c("Normal", "Cancer"))

# -----------------------------
# 3. Run CellDMC
# -----------------------------

covariate_model <- model.matrix(
  ~ Age,
  data = pheno
)

# Optional:
# Additional covariates such as cohort or batch variables
# may be included if appropriate for the study design.
#
# Example:
# covariate_model <- model.matrix(
#   ~ Age + Sex,
#   data = pheno
# )

cdmc_results <- CellDMC(
  beta.m = beta_mat,
  pheno.v = pheno$Status,
  frac.m = as.matrix(cell_fractions),
  cov.mod = covariate_model
)

write.csv(
  cdmc_results,
  "CellDMC_all_CpGs.csv"
)

# -----------------------------
# 4. Annotate CellDMC results
# -----------------------------

ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

ann450k_sub <- ann450k[
  match(rownames(cdmc_results), ann450k$Name),
]

cdmc_annotated <- cbind(
  cdmc_results,
  ann450k_sub
)

write.csv(
  cdmc_annotated,
  "CellDMC_all_CpGs_annotated.csv"
)

# -----------------------------
# 5. Extract significant CellDMC CpGs
# -----------------------------

cdmc_dmc <- cdmc_annotated %>%
  filter(dmct.DMC == 1)

write.csv(
  cdmc_dmc,
  "CellDMC_significant_DMCs_annotated.csv"
)

# -----------------------------
# 6. Extract endothelial-specific DMCs
# -----------------------------

endothelial_dmc <- cdmc_dmc %>%
  filter(dmct.Endothelial != 0)

write.csv(
  endothelial_dmc,
  "CellDMC_endothelial_specific_DMCs_annotated.csv"
)

# -----------------------------
# 7. Optional summary counts
# -----------------------------

table(cdmc_dmc$dmct.Endothelial != 0)

table(endothelial_dmc$dmct.Endothelial)