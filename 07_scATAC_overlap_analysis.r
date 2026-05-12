############################################################
# scATAC-seq overlap analysis
#
# Authors:
# Barbara Karakyriakou, Brock C. Christensen
#
# Description:
# Compares TEC-specific methylation alterations with
# endothelial cell differential chromatin accessibility peaks.
# CpG coordinates are lifted from hg19 to hg38 and overlapped
# with ATAC peaks. Gene-level overlap is assessed using a
# permutation test.
############################################################

library(rtracklayer)
library(GenomicRanges)
library(dplyr)
library(readr)
library(stringr)
library(tibble)
library(ggplot2)

# -----------------------------
# 1. Required input files
# -----------------------------

# tec_cpgs:
#   Annotated TEC-specific CpGs.
#   Required columns:
#   chr, pos, UCSC_RefGene_Name
#
# atac_peaks:
#   Differential ATAC peak results.
#   Required columns:
#   peak, p_val_fdr, avg_log2FC, significant, direction, gene_name
#   Peak format should be: chr-start-end
#
# chain_file:
#   UCSC hg19ToHg38.over.chain file for coordinate liftover.

tec_cpgs <- read.csv(
  "TEC_specific_CpGs_annotated.csv",
  row.names = 1
)

atac_peaks <- read.csv(
  "endothelial_scATAC_differential_peaks.csv"
)

chain_file <- "hg19ToHg38.over.chain"

# -----------------------------
# 2. Prepare ATAC peak GRanges, hg38
# -----------------------------

atac_peaks <- atac_peaks %>%
  mutate(
    chr = sub("-.*", "", peak),
    start = as.integer(sub("^[^-]+-([0-9]+)-.*", "\\1", peak)),
    end = as.integer(sub("^[^-]+-[0-9]+-([0-9]+)$", "\\1", peak))
  )

atac_gr <- makeGRangesFromDataFrame(
  atac_peaks,
  keep.extra.columns = TRUE,
  seqnames.field = "chr",
  start.field = "start",
  end.field = "end"
)

# -----------------------------
# 3. Prepare CpG GRanges, hg19
# -----------------------------

cpg_hg19 <- makeGRangesFromDataFrame(
  tec_cpgs %>% mutate(end = pos),
  keep.extra.columns = TRUE,
  seqnames.field = "chr",
  start.field = "pos",
  end.field = "end"
)

# -----------------------------
# 4. Liftover CpGs from hg19 to hg38
# -----------------------------

chain <- import.chain(chain_file)

cpg_hg38_list <- liftOver(cpg_hg19, chain)

mapped_idx <- lengths(cpg_hg38_list) == 1

cpg_hg38 <- unlist(cpg_hg38_list[mapped_idx])
tec_cpgs_mapped <- tec_cpgs[mapped_idx, ]

cat(
  "CpGs successfully lifted:",
  length(cpg_hg38),
  "of",
  nrow(tec_cpgs),
  "\n"
)

# -----------------------------
# 5. Direct CpG–ATAC peak overlap
# -----------------------------

hits_direct <- findOverlaps(
  atac_gr,
  cpg_hg38,
  ignore.strand = TRUE
)

direct_overlap_df <- data.frame(
  atac_peaks[queryHits(hits_direct), ],
  CpG = rownames(tec_cpgs_mapped)[subjectHits(hits_direct)],
  tec_cpgs_mapped[subjectHits(hits_direct), ]
)

write_csv(
  direct_overlap_df,
  "ATAC_TEC_CpG_direct_overlaps_hg38.csv"
)

# -----------------------------
# 6. Window-based CpG–ATAC peak overlap
# -----------------------------

window_bp <- 2000L

atac_window_gr <- atac_gr

ranges(atac_window_gr) <- IRanges(
  start = pmax(1, start(atac_gr) - window_bp),
  end = end(atac_gr) + window_bp
)

hits_window <- findOverlaps(
  atac_window_gr,
  cpg_hg38,
  ignore.strand = TRUE
)

window_overlap_df <- data.frame(
  atac_peaks[queryHits(hits_window), ],
  CpG = rownames(tec_cpgs_mapped)[subjectHits(hits_window)],
  tec_cpgs_mapped[subjectHits(hits_window), ]
)

write_csv(
  window_overlap_df,
  "ATAC_TEC_CpG_window_overlaps_hg38.csv"
)

# -----------------------------
# 7. CpG-level overlap summary
# -----------------------------

overlap_summary <- data.frame(
  total_TEC_CpGs = nrow(tec_cpgs),
  lifted_TEC_CpGs = length(cpg_hg38),
  direct_overlap_peaks = length(unique(queryHits(hits_direct))),
  direct_overlap_CpGs = length(unique(subjectHits(hits_direct))),
  window_bp = window_bp,
  window_overlap_peaks = length(unique(queryHits(hits_window))),
  window_overlap_CpGs = length(unique(subjectHits(hits_window)))
)

write_csv(
  overlap_summary,
  "ATAC_TEC_CpG_overlap_summary.csv"
)

print(overlap_summary)

# -----------------------------
# 8. Gene-level overlap
# -----------------------------

extract_genes <- function(x) {
  x %>%
    strsplit(";") %>%
    unlist() %>%
    trimws() %>%
    unique() %>%
    .[!is.na(.) & . != ""]
}

tec_genes <- extract_genes(tec_cpgs$UCSC_RefGene_Name)

atac_genes <- atac_peaks %>%
  filter(significant == TRUE) %>%
  pull(gene_name) %>%
  unique() %>%
  na.omit()

shared_genes <- intersect(tec_genes, atac_genes)

write_csv(
  data.frame(shared_gene = shared_genes),
  "TEC_ATAC_shared_genes.csv"
)

cat("TEC genes:", length(tec_genes), "\n")
cat("ATAC genes:", length(atac_genes), "\n")
cat("Shared genes:", length(shared_genes), "\n")

# -----------------------------
# 9. Gene universe for enrichment
# -----------------------------

# all_cellDMC_results:
#   Annotated CellDMC result table used to define the tested
#   methylation gene universe.

all_cellDMC_results <- read.csv(
  "CellDMC_all_CpGs_annotated.csv",
  row.names = 1
)

gene_universe <- extract_genes(all_cellDMC_results$UCSC_RefGene_Name)

tec_genes_u <- intersect(tec_genes, gene_universe)
atac_genes_u <- intersect(atac_genes, gene_universe)

observed_overlap <- length(intersect(tec_genes_u, atac_genes_u))

# -----------------------------
# 10. Permutation test
# -----------------------------

set.seed(123)

n_perm <- 10000

perm_overlaps <- replicate(n_perm, {
  random_genes <- sample(
    gene_universe,
    size = length(tec_genes_u),
    replace = FALSE
  )
  length(intersect(random_genes, atac_genes_u))
})

p_empirical <- (sum(perm_overlaps >= observed_overlap) + 1) / (n_perm + 1)

expected_overlap <- mean(perm_overlaps)
fold_enrichment <- observed_overlap / expected_overlap
z_score <- (observed_overlap - mean(perm_overlaps)) / sd(perm_overlaps)

perm_summary <- data.frame(
  universe_size = length(gene_universe),
  tec_gene_count = length(tec_genes_u),
  atac_gene_count = length(atac_genes_u),
  observed_overlap = observed_overlap,
  expected_overlap = expected_overlap,
  sd_expected = sd(perm_overlaps),
  fold_enrichment = fold_enrichment,
  z_score = z_score,
  empirical_p_value = p_empirical,
  n_permutations = n_perm
)

write_csv(
  perm_summary,
  "TEC_ATAC_gene_overlap_permutation_summary.csv"
)

write_csv(
  data.frame(perm_overlap = perm_overlaps),
  "TEC_ATAC_gene_overlap_permutations.csv"
)

print(perm_summary)

# -----------------------------
# 11. Optional null distribution plot
# -----------------------------

perm_plot_df <- data.frame(overlap = perm_overlaps)

p_perm <- ggplot(
  perm_plot_df,
  aes(x = overlap)
) +
  geom_histogram(
    bins = 50,
    fill = "grey80",
    color = "black"
  ) +
  geom_vline(
    xintercept = observed_overlap,
    linetype = "dashed",
    linewidth = 1
  ) +
  labs(
    title = "Permutation test for TEC gene overlap with ATAC genes",
    x = "Overlap with ATAC genes in random gene sets",
    y = "Count"
  ) +
  theme_classic(base_size = 14)

ggsave(
  "TEC_ATAC_gene_overlap_permutation_plot.pdf",
  p_perm,
  width = 6,
  height = 5
)