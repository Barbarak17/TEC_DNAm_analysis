# TEC_DNAm_analysis

## Description

This repository contains analysis scripts associated with the manuscript:

**Tumor endothelial cell-specific DNA methylation alterations in breast cancer**

The workflow reconstructs lineage-specific DNA methylation alterations from bulk Illumina HumanMethylation450K array data using HiTIMED deconvolution and CellDMC interaction modeling. Additional analyses include bulk EWAS, pathway enrichment, batch correction, and integration with endothelial single-cell ATAC-seq datasets.

---

# Repository Structure

```text
01_preprocessing.R
02_batch_correction_and_UMAP.R
03_HiTIMED_deconvolution.R
04_bulk_EWAS.R
05_CellDMC_analysis.R
06_pathway_enrichment.R
07_scATAC_overlap_analysis.R
```

---

# Analysis Overview

## 1. Preprocessing

`01_preprocessing.R`

* Reads IDAT files
* Performs functional normalization using `preprocessFunnorm`
* Removes low-quality samples and probes
* Filters SNP and cross-reactive probes
* Exports processed beta matrices

---

## 2. Batch correction and visualization

`02_batch_correction_and_UMAP.R`

* Performs ComBat batch correction
* Preserves biological status variables during harmonization
* Generates UMAP visualizations before and after correction

---

## 3. HiTIMED deconvolution

`03_HiTIMED_deconvolution.R`

* Estimates cell-type proportions from bulk methylation data
* Generates CellDMC-compatible deconvolution matrices

---

## 4. Bulk EWAS

`05_bulk_EWAS.R`

* Performs limma-based differential methylation analysis
* Supports:

  * age-adjusted models
  * age and cell-type-adjusted models
* Generates volcano plots and summary tables

---

## 5. CellDMC analysis

`04_CellDMC_analysis.R`

* Identifies cell-type-specific differentially methylated CpGs
* Integrates HiTIMED-derived cell fractions
* Exports annotated endothelial-specific DMCs

---

## 6. Pathway enrichment analysis

`06_pathway_enrichment.R`

* Performs GO enrichment using `missMethyl`
* Focuses on angiogenesis, endothelial, and VEGF-related pathways
* Generates targeted gene-count enrichment plots

---

## 7. scATAC-seq overlap analysis

`07_scATAC_overlap_analysis.R`

* Lifts CpG coordinates from hg19 to hg38
* Identifies overlaps between TEC-specific CpGs and endothelial differential ATAC peaks
* Performs gene-level overlap and permutation enrichment analyses

---

# Input Data Requirements

Scripts assume the presence of:

* DNA methylation beta matrices
* sample metadata tables
* HiTIMED deconvolution outputs
* annotated CellDMC result tables
* differential scATAC-seq peak tables

Input file paths are provided as placeholders and should be adapted locally.

---

# Software Requirements

Analyses were performed in R using packages including:

* minfi
* limma
* ENmix
* sva
* HiTIMED
* EpiDISH
* missMethyl
* GenomicRanges
* rtracklayer
* ggplot2
* EnhancedVolcano

---

# Data Availability

Public datasets used in this study were obtained from:

* TCGA
* GEO
* GTEx

Accession numbers are provided in the associated manuscript.

---

# Citation

If using these scripts or derivative workflows, please cite the associated manuscript.
