############################################################
# HiTIMED deconvolution
#
# Authors:
# Barbara Karakyriakou
#
# Description:
# Estimates cell-type proportions from bulk DNA methylation
# beta values using HiTIMED. Outputs deconvolution estimates
# for downstream visualization and CellDMC analysis.
############################################################

library(HiTIMED)
library(EpiDISH)
library(dplyr)

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
#
# Required columns:
#   Sample_Name
#   Status
#
# Status should indicate:
#   Cancer
#   Normal

beta_mat <- readRDS("processed_beta_matrix.rds")
pheno <- read.csv("sample_metadata.csv")

# Align phenotype rows to beta matrix columns
pheno <- pheno[match(colnames(beta_mat), pheno$Sample_Name), ]

stopifnot(all(colnames(beta_mat) == pheno$Sample_Name))

# -----------------------------
# 2. Optional: load HiTIMED reference library
# -----------------------------

# Load the HiTIMED reference library if required by your
# local installation or workflow.
#
# Example:
# load("HiTIMED_Library.RData")

# -----------------------------
# 3. Split tumor and normal samples
# -----------------------------

normal_samples <- pheno$Sample_Name[pheno$Status == "Normal"]
cancer_samples <- pheno$Sample_Name[pheno$Status == "Cancer"]

beta_normal <- beta_mat[, colnames(beta_mat) %in% normal_samples]
beta_cancer <- beta_mat[, colnames(beta_mat) %in% cancer_samples]

# -----------------------------
# 4. Run HiTIMED separately for visualization
# -----------------------------

# These estimates may be used for summary plots comparing
# deconvolved cell-type composition between cancer and normal
# samples.

hitimed_normal <- HiTIMED_deconvolution(
  tumor_beta = beta_normal,
  tumor_type = "BRCA",
  tissue_type = "normal",
  h = 3
)

hitimed_cancer <- HiTIMED_deconvolution(
  tumor_beta = beta_cancer,
  tumor_type = "BRCA",
  tissue_type = "tumor",
  h = 3
)

hitimed_visualization <- rbind(
  hitimed_cancer,
  hitimed_normal
)

write.csv(
  hitimed_visualization,
  "HiTIMED_deconvolution_cancer_normal_visualization.csv"
)