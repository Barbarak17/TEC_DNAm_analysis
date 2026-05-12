############################################################
# Batch correction and UMAP visualization
#
# Authors:
# Barbara Karakyriakou
# Description:
# Performs ComBat batch correction across methylation
# cohorts and visualizes sample clustering before and
# after correction using UMAP.
############################################################

library(sva)
library(umap)
library(plotly)
library(ENmix)
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
#   Status (Cancer/Normal)
#   Cohort

beta_mat <- readRDS("merged_beta_matrix.rds")
pheno <- read.csv("sample_metadata.csv")

# -----------------------------
# 2. Align metadata
# -----------------------------

pheno <- pheno[
  match(colnames(beta_mat), pheno$Sample_Name),
]

stopifnot(all(colnames(beta_mat) == pheno$Sample_Name))

# -----------------------------
# 3. UMAP before batch correction
# -----------------------------

set.seed(123)

umap_before <- umap(
  t(beta_mat),
  n_components = 2
)

umap_before_df <- data.frame(
  UMAP1 = umap_before$layout[,1],
  UMAP2 = umap_before$layout[,2],
  Status = pheno$Status,
  Cohort = pheno$Cohort
)

fig_before <- plot_ly(
  umap_before_df,
  x = ~UMAP1,
  y = ~UMAP2,
  color = ~Status,
  symbol = ~Cohort,
  type = "scatter",
  mode = "markers"
) %>%
  layout(
    title = "Before ComBat batch correction",
    plot_bgcolor = "#FFFFFF",
    xaxis = list(showgrid = FALSE),
    yaxis = list(showgrid = FALSE)
  )

fig_before

# -----------------------------
# 4. ComBat batch correction
# -----------------------------

# Convert beta values to M-values
m_values <- B2M(beta_mat)

# Preserve biological variable of interest
combat_model <- model.matrix(
  ~ Status,
  data = pheno
)

m_values_combat <- ComBat(
  dat = m_values,
  batch = pheno$Cohort,
  mod = combat_model,
  mean.only = FALSE
)

# Convert back to beta values
beta_combat <- M2B(m_values_combat)

saveRDS(
  beta_combat,
  "beta_matrix_combat_corrected.rds"
)

# -----------------------------
# 5. UMAP after batch correction
# -----------------------------

set.seed(123)

umap_after <- umap(
  t(beta_combat),
  n_components = 2
)

umap_after_df <- data.frame(
  UMAP1 = umap_after$layout[,1],
  UMAP2 = umap_after$layout[,2],
  Status = pheno$Status,
  Cohort = pheno$Cohort
)

fig_after <- plot_ly(
  umap_after_df,
  x = ~UMAP1,
  y = ~UMAP2,
  color = ~Status,
  symbol = ~Cohort,
  type = "scatter",
  mode = "markers"
) %>%
  layout(
    title = "After ComBat batch correction",
    plot_bgcolor = "#FFFFFF",
    xaxis = list(showgrid = FALSE),
    yaxis = list(showgrid = FALSE)
  )

fig_after