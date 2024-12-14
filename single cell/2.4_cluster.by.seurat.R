# Title: Single-Cell RNA-Seq Analysis with Seurat
# Author:  Group 7
# Date: 2024-11-22
# Project: Final Project for BMD_ENG_311
# Description:
#     This script uses Seurat 5.0.1 to process single-cell RNA-Seq data.
#     Steps include quality control, normalization, dimensionality reduction,
#     clustering, and differential expression analysis.
#     
# Dependencies:
#     - R (>= 4.1.0)
#     - Seurat v5
#     - ggplot2, dplyr (for additional plotting and data manipulation)
#     
# Input:
#     - scRNA-Seq data matrix (e.g., counts or UMI matrix)
#     
# Output:
#     - Processed Seurat object
#     - Clustering results
#     - Plots for visualization (e.g., UMAP, feature expression)

# set working dir
setwd("/projects/e30836/project/group7/sc/counts/")

# load packages
library(dplyr)
library(Seurat)
library(patchwork)

cells <- read.table(file = 'GSM6939133_HFC1_barcodes.tsv.gz', sep = '\t', header = TRUE)
features <- read.table(file = 'GSM6939133_HFC1_features.tsv.gz', sep = '\t', header = TRUE)
gsm <- ReadMtx("GSM6939133_HFC1_matrix.mtx.gz",  cells = "GSM6939133_HFC1_barcodes.tsv.gz", features = "GSM6939133_HFC1_features.tsv.gz")
gsm_1 <- CreateSeuratObject(counts = gsm, project = "HFC1")
head(gsm_1)

gsm_1[["percent.mt"]] <- PercentageFeatureSet(gsm_1, pattern = "^MT-")
VlnPlot(gsm_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(gsm_1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(gsm_1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

hfc1_v <- subset(gsm_1, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

hfc1_v <- NormalizeData(hfc1_v, normalization.method = "LogNormalize", scale.factor = 10000)

hfc1_v <- FindVariableFeatures(hfc1_v, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(hfc1_v), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(hfc1_v)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(hfc1_v)
hfc1_v <- ScaleData(hfc1_v, features = all.genes)

hfc1_v <- RunPCA(hfc1_v, features = VariableFeatures(object = hfc1_v))
# Examine and visualize PCA results a few different ways
print(hfc1_v[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(hfc1_v, dims = 1:2, reduction = "pca")

DimPlot(hfc1_v, reduction = "pca") + NoLegend()
DimHeatmap(hfc1_v, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(hfc1_v, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(hfc1_v)

hfc1_v <- FindNeighbors(hfc1_v, dims = 1:10)
hfc1_v <- FindClusters(hfc1_v, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(hfc1_v), 5)
hfc1_v <- RunUMAP(hfc1_v, dims = 1:10)
DimPlot(hfc1_v, reduction = "umap")
