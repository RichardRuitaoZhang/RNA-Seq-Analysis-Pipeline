# 🧬 RNA-Seq Analysis Pipeline

This repository contains a complete, modular pipeline for **Bulk** and **Single-Cell RNA-Seq** analysis, integrating Shell-based workflow control and R-based statistical analysis.
It provides an end-to-end workflow from **raw sequencing data** to **differential expression** and **functional enrichment**, supporting both **bulk RNA-Seq** and **10x Genomics single-cell RNA-Seq** data.

---

### Repository Structure

```
RNA-Seq-Analysis-Pipeline/
├── RFunctions/         # Contains custom reusable R function definitions supporting both bulk and single-cell workflows
├── RScripts/      # Core R functions used by pipeline scripts
├── bulk/          # Complete bulk RNA-Seq pipeline (FASTQ → DE → enrichment)
├── single cell/   # Single-cell RNA-Seq pipeline (10x Genomics + Seurat)
├── paper/         # Example datasets from published studies
└── README.md
```

---

### Project Overview

This pipeline automates **RNA-Seq data analysis** for both bulk and single-cell experiments.
It separates computational steps into modular scripts and reusable R functions, allowing flexible execution on local or HPC environments.

**Key Features:**

* Supports both **bulk** and **single-cell** RNA-Seq workflows
* Modular **Shell + R** integration for clarity and reproducibility
* Fully reproducible with all outputs serialized as `.RData`
* Includes **example datasets** for demonstration under `paper/`

---

## Bulk RNA-Seq Analysis

### 🔹 Overview

The **Bulk RNA-Seq module** provides a complete end-to-end workflow for standard transcriptome profiling — from **raw FASTQ files** to **differential expression** and **functional enrichment analysis**.

This workflow is fully modular, allowing users to execute or debug each step independently.

---

### 🔹 Workflow Summary

| Step            | Script                                                                                     | Description                                       |
| --------------- | ------------------------------------------------------------------------------------------ | ------------------------------------------------- |
| 1.1             | `1.1_download.bulk.command.sh`                                                             | Download raw FASTQ data                           |
| 1.2             | `1.2_quality.control.sh`                                                                   | Quality control using FastQC                      |
| 1.2.1           | `1.2.1_trimming.to.remove.lq.and.adapters.sh`                                              | Adapter trimming and removal of low-quality reads |
| 1.3             | `1.3_alignment.via.STAR.sh`                                                                | Alignment to reference genome via STAR            |
| 1.3.0 / 1.3.0.1 | Reference genome download and indexing (hg38)                                              |                                                   |
| 1.3.1           | `1.3.1_post_alignment.qc.sh`                                                               | Post-alignment QC                                 |
| 1.4             | `1.4_samtools.sh`                                                                          | BAM sorting, indexing, and filtering              |
| 1.5             | `1.5_RSEM.sh`, `1.5.1_construct.expression.matrix*.sh`, `1.5.2_import.exp.results.merge.R` | Expression quantification and matrix construction |
| 1.6             | `1.6.1_Generate.colData.for.DESeq2.Analysis.R`                                             | Generate metadata (colData) for DESeq2            |
| 1.7             | `1.7_DESeq2.R`, `1.7_DESeq2.Rmd`                                                           | Differential expression analysis                  |
| 1.8             | `1.8.1_gene.ontology.R`, `1.8.2_KEGG.R`, `1.8.3_gene.set.enrichment.analysis.R`            | Functional enrichment analysis (GO/KEGG/GSEA)     |
| 1.8.4           | `1.8.4_heatmap.R`                                                                          | Expression heatmap visualization                  |
| 1.8.5           | `1.8.5_tri_conditions.R`                                                                   | Multi-condition comparison analysis               |

> ⚙️ **Note:** Variants like `1.5.1_construct.expression.matrix_PC9.sh` and `..._SKMEL28.sh` handle cell line–specific expression construction.

---

### 🔹 Dependencies

**External tools:** `FastQC`, `Trim Galore`, `STAR`, `Samtools`, `RSEM`
**R packages:** `DESeq2`, `clusterProfiler`, `org.Hs.eg.db`, `ggplot2`, `enrichplot`, `fgsea`, `dplyr`, `tidyr`

---

### 🔹 Example Data

The `paper/cancers-16-01001-with-cover.pdf` dataset demonstrates the **bulk RNA-Seq** workflow, serving as a test and educational dataset for this pipeline.

---

## Single-Cell RNA-Seq Analysis

### 🔹 Overview

The **Single-Cell RNA-Seq module** automates 10x Genomics data processing and downstream Seurat analysis.
It provides a complete workflow from **FASTQ** to **cell clustering and dataset integration**.

---

### 🔹 Workflow Summary

| Step          | Script                                                | Description                                                         |
| ------------- | ----------------------------------------------------- | ------------------------------------------------------------------- |
| 2.1           | `2.1_download.sc.command.sh`                          | Download single-cell FASTQ files                                    |
| 2.2.1 / 2.2.2 | `construct.gtf.for.ref.sh`, `construct.ref.genome.sh` | Build CellRanger reference genome and GTF annotation                |
| 2.3           | `2.3_cellranger.count.sh`                             | Perform read alignment and quantification using CellRanger          |
| 2.4           | `2.4_cluster.by.seurat.R`                             | Clustering, dimensionality reduction, and visualization with Seurat |
| 2.5           | `2.5_scRNA_integration.R`                             | Multi-sample integration and batch correction                       |

---

### 🔹 Dependencies

**External tools:** `CellRanger`
**R packages:** `Seurat`, `ggplot2`, `dplyr`, `tidyr`

---

### 🔹 Example Data

The dataset `paper/cancers-2874917-supplementary.pdf` demonstrates the **single-cell RNA-Seq** pipeline.

---

## 🔢 RScripts and RData

### `RScripts/`

Contains all **core R functions** used by both bulk and single-cell pipelines. These include:

* Differential expression (DESeq2 wrappers)
* Enrichment analysis (GO, KEGG, GSEA)
* Visualization (MA, Volcano, Heatmap, Bar, Bubble)

### `RData/`

Contains serialized `.RData` objects for every functions wrote in `RScripts/`, enabling reproducibility and fast reload.

---

## Paper Directory

Contains published example datasets used for demonstration and validation.

| File                                | Dataset Type        | Description                            |
| ----------------------------------- | ------------------- | -------------------------------------- |
| `cancers-16-01001-with-cover.pdf`   | Bulk RNA-Seq        | Example dataset for bulk analysis      |
| `cancers-2874917-supplementary.pdf` | Single-cell RNA-Seq | Example dataset for scRNA-seq workflow |

---

## 🚀 Quick Start

```bash
# Clone repository
git clone https://github.com/RichardRuitaoZhang/RNA-Seq-Analysis-Pipeline.git
cd RNA-Seq-Analysis-Pipeline

# Run bulk RNA-Seq pipeline
bash bulk/1.1_download.bulk.command.sh
bash bulk/1.3_alignment.via.STAR.sh
Rscript bulk/1.7_DESeq2.R

# Run single-cell RNA-Seq pipeline
bash single\ cell/2.3_cellranger.count.sh
Rscript single\ cell/2.4_cluster.by.seurat.R
```

---

## Design Philosophy

* **Modularized:** Each step is independent and reusable.
* **Reproducible:** All results stored as `.RData` objects.
* **Extensible:** Can integrate new tools or data types.
* **Cross-platform:** Designed for Linux and HPC environments.

---

# Author & Course Information

**Author:** Ruitao Zhang
**Course:** BME311 – Computational Genomics, Fall 2024
**Instructor:** Dr. Zhe Ji
Northwestern University, Evanston, IL
