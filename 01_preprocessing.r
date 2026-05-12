############################################################
# DNA methylation preprocessing from IDAT files
#
# Authors:
# Barbara Karakyriakou
#
# Description:
# Preprocesses Illumina 450K methylation IDAT files using minfi,
# performs functional normalization, removes low-quality samples
# and CpGs, filters SNP/cross-reactive probes using the HiTIMED
# annotation, and exports a beta-value matrix.
############################################################

library(minfi)
library(ENmix)
library(HiTIMED)

# -----------------------------
# 1. Required input files
# -----------------------------

# idat_dir:
#   Directory containing IDAT files and sample sheet.
#
# Sample sheet should include:
#   Sample_Name
#   TissueType or other phenotype variables
#
# HM450_annotation:
#   Annotation object containing probe filtering information.
#   Required column:
#   MASK_general

idat_dir <- "path_to_idat_directory"

# Load annotation object used for probe filtering
# Example:
# load("HM450_annotation.RData")

# -----------------------------
# 2. Read IDAT files
# -----------------------------

targets <- read.metharray.sheet(idat_dir)

rgSet <- read.metharray.exp(
  targets = targets,
  extended = TRUE,
  force = TRUE
)

sampleNames(rgSet) <- targets$Sample_Name

# Add phenotype information
pheno <- pData(rgSet)
pheno$TissueType <- targets$TissueType

# -----------------------------
# 3. Quality control
# -----------------------------

# Generate ENmix control probe QC plots
plotCtrl(rgSet)

# Extract QC information
qc <- QCinfo(rgSet)

# Optional: save QC object
saveRDS(qc, "QCinfo_ENmix.rds")

# -----------------------------
# 4. Functional normalization
# -----------------------------

mSet_norm <- preprocessFunnorm(rgSet)

# -----------------------------
# 5. Remove SNP/cross-reactive/failed probes
# -----------------------------

filtered_probes <- HM450_annotation$probeID[
  HM450_annotation$MASK_general == TRUE
]

mSet_filt <- mSet_norm[
  !(featureNames(mSet_norm) %in% filtered_probes),
]

# -----------------------------
# 6. Extract beta values
# -----------------------------

beta_mat <- getBeta(mSet_filt)

# -----------------------------
# 7. Remove low-quality samples and CpGs
# -----------------------------

beta_clean <- beta_mat[
  !(rownames(beta_mat) %in% qc$badCpG),
  !(colnames(beta_mat) %in% qc$badSample)
]

# -----------------------------
# 8. Save processed beta matrix
# -----------------------------

saveRDS(beta_clean, "processed_beta_matrix.rds")