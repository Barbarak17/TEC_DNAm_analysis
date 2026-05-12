############################################################
# Targeted GO enrichment analysis
#
# Authors:
# Barbara Karakyriakou, Brock C. Christensen
#
# Description:
# Performs GO enrichment analysis for TEC-specific CpGs using
# missMethyl, extracts angiogenesis/endothelial/VEGF-related
# biological process terms, and generates the final gene-count
# dot plot used for visualization.
############################################################

library(missMethyl)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(tibble)
library(readr)

# -----------------------------
# 1. Required input files
# -----------------------------

# tec_cpgs:
#   Dataframe containing significant TEC-specific CpGs.
#   Row names should be CpG probe IDs.
#
# all_cpgs:
#   Dataframe containing all CpGs tested in CellDMC.
#   Row names should be CpG probe IDs.

tec_cpgs <- read.csv(
  "TEC_specific_CpGs.csv",
  row.names = 1
)

all_cpgs <- read.csv(
  "CellDMC_all_CpGs_annotated.csv",
  row.names = 1
)

sig_cpg <- rownames(tec_cpgs)
all_cpg <- rownames(all_cpgs)

# -----------------------------
# 2. Load 450K annotation
# -----------------------------

ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# -----------------------------
# 3. Run GO enrichment using gometh
# -----------------------------

go_res <- gometh(
  sig.cpg = sig_cpg,
  all.cpg = all_cpg,
  collection = "GO",
  plot.bias = FALSE,
  prior.prob = TRUE,
  anno = ann450k,
  equiv.cpg = TRUE,
  fract.counts = TRUE,
  genomic.features = c(
    "ALL",
    "TSS200",
    "TSS1500",
    "Body",
    "1stExon",
    "3'UTR",
    "5'UTR",
    "ExonBnd"
  ),
  sig.genes = TRUE
) %>%
  as_tibble()

write_csv(
  go_res,
  "GO_enrichment_all_terms.csv"
)

# -----------------------------
# 4. Filter angiogenesis/endothelial/VEGF-related GO terms
# -----------------------------

term_pattern <- paste(
  "ANGIOGEN",
  "VASCULOGEN",
  "VEGF",
  "VASCULAR ENDOTHELIAL GROWTH",
  "BLOOD VESSEL",
  "ENDOTHELIAL",
  "ENDOTHELIUM",
  "TUBE",
  "SPROUT",
  "ENDOTHELIAL CELL MIGRATION",
  "SHEAR STRESS",
  "BLOOD",
  "ARTERY",
  sep = "|"
)

go_targeted <- go_res %>%
  filter(ONTOLOGY == "BP") %>%
  filter(!is.na(FDR)) %>%
  filter(str_detect(toupper(TERM), term_pattern)) %>%
  arrange(FDR)

write_csv(
  go_targeted,
  "GO_targeted_angiogenesis_endothelial_VEGF_terms.csv"
)

# -----------------------------
# 5. Extract genes from targeted GO terms
# -----------------------------

gene_column <- if ("SigGenesInSet" %in% colnames(go_targeted)) {
  "SigGenesInSet"
} else {
  "DEgenes"
}

go_term_gene_long <- go_targeted %>%
  select(TERM, genes = all_of(gene_column)) %>%
  mutate(genes = as.character(genes)) %>%
  separate_rows(genes, sep = "[,;]") %>%
  mutate(gene = str_trim(genes)) %>%
  filter(!is.na(gene), gene != "") %>%
  distinct(TERM, gene)

genes_per_term <- go_term_gene_long %>%
  group_by(TERM) %>%
  summarise(
    genes = paste(sort(unique(gene)), collapse = "; "),
    .groups = "drop"
  )

# -----------------------------
# 6. Convert gene symbols to Entrez IDs
# -----------------------------

go_sets_symbol <- split(
  go_term_gene_long$gene,
  go_term_gene_long$TERM
)

go_sets_symbol <- lapply(go_sets_symbol, unique)

all_symbols <- unique(unlist(go_sets_symbol))

symbol_to_entrez <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = all_symbols,
  keytype = "SYMBOL",
  columns = "ENTREZID"
)

go_sets_entrez <- lapply(go_sets_symbol, function(gene_set) {
  unique(symbol_to_entrez$ENTREZID[
    match(gene_set, symbol_to_entrez$SYMBOL)
  ])
})

go_sets_entrez <- lapply(go_sets_entrez, function(x) x[!is.na(x)])
go_sets_entrez <- go_sets_entrez[lengths(go_sets_entrez) >= 5]

# -----------------------------
# 7. Run targeted gene set enrichment using gsameth
# -----------------------------

targeted_go_res <- gsameth(
  sig.cpg = sig_cpg,
  all.cpg = all_cpg,
  collection = go_sets_entrez
)

targeted_go_df <- targeted_go_res %>%
  as.data.frame() %>%
  rownames_to_column("TERM") %>%
  left_join(genes_per_term, by = "TERM") %>%
  arrange(FDR)

write_csv(
  targeted_go_df,
  "targeted_GO_angiogenesis_endothelial_VEGF_results.csv"
)

# -----------------------------
# 8. Generate final gene-count dot plot
# -----------------------------

plot_df <- targeted_go_df %>%
  filter(!is.na(FDR)) %>%
  slice_head(n = 20) %>%
  mutate(
    neglog10FDR = -log10(FDR + 1e-300),
    pathway_label = str_trunc(TERM, 70)
  )

p_gene_count <- ggplot(
  plot_df,
  aes(
    x = DE,
    y = reorder(pathway_label, DE),
    color = neglog10FDR,
    size = DE
  )
) +
  geom_point(alpha = 0.85) +
  scale_color_gradient(
    low = "blue",
    high = "red"
  ) +
  labs(
    x = "Gene Counts",
    y = "Angiogenesis, Endothelial & VEGF GO terms",
    color = expression(-log[10]("FDR")),
    size = "# DM genes",
    title = "Targeted GO enrichment of TEC-specific CpGs"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.y = element_text(size = 9),
    legend.position = "right"
  )

p_gene_count

ggsave(
  filename = "targeted_GO_gene_count_dotplot.pdf",
  plot = p_gene_count,
  width = 7,
  height = 6
)