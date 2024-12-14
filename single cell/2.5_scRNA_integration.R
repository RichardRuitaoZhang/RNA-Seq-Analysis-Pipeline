library(Seurat)
library(ggplot2)
library(patchwork)

# Load and preprocess pc9_u
gsm_u <- ReadMtx("GSM4932159_sample1_matrix.mtx.gz", cells = "GSM4932159_sample1_barcodes.tsv.gz", features = "GSM4932159_sample1_genes.tsv.gz")
pc9_u <- CreateSeuratObject(counts = gsm_u, project = "PC9_Untreated")
pc9_u[["percent.mt"]] <- PercentageFeatureSet(pc9_u, pattern = "^MT-")
pc9_u <- subset(pc9_u, subset = nFeature_RNA > 100 & nFeature_RNA < 5000 & percent.mt < 10)
pc9_u <- NormalizeData(pc9_u, assay = "RNA")
pc9_u <- FindVariableFeatures(pc9_u, assay = "RNA")
pc9_u <- ScaleData(pc9_u, assay = "RNA")

# Load and preprocess pc9_24
gsm_24 <- ReadMtx("GSM4932160_sample2_matrix.mtx.gz", cells = "GSM4932160_sample2_barcodes.tsv.gz", features = "GSM4932160_sample2_genes.tsv.gz")
pc9_24 <- CreateSeuratObject(counts = gsm_24, project = "PC9_24")
pc9_24[["percent.mt"]] <- PercentageFeatureSet(pc9_24, pattern = "^MT-")
pc9_24 <- subset(pc9_24, subset = nFeature_RNA > 100 & nFeature_RNA < 5000 & percent.mt < 10)
pc9_24 <- NormalizeData(pc9_24, assay = "RNA")
pc9_24 <- FindVariableFeatures(pc9_24, assay = "RNA")
pc9_24 <- ScaleData(pc9_24, assay = "RNA")

# Load and preprocess pc9_48
gsm_48 <- ReadMtx("GSM4932161_sample3_matrix.mtx.gz", cells = "GSM4932161_sample3_barcodes.tsv.gz", features = "GSM4932161_sample3_genes.tsv.gz")
pc9_48 <- CreateSeuratObject(counts = gsm_48, project = "PC9_48")
pc9_48[["percent.mt"]] <- PercentageFeatureSet(pc9_48, pattern = "^MT-")
pc9_48 <- subset(pc9_48, subset = nFeature_RNA > 100 & nFeature_RNA < 5000 & percent.mt < 10)
pc9_48 <- NormalizeData(pc9_48, assay = "RNA")
pc9_48 <- FindVariableFeatures(pc9_48, assay = "RNA")
pc9_48 <- ScaleData(pc9_48, assay = "RNA")

# Load and preprocess pc9_72
gsm_72 <- ReadMtx("GSM4932162_sample4_matrix.mtx.gz", cells = "GSM4932162_sample4_barcodes.tsv.gz", features = "GSM4932162_sample4_genes.tsv.gz")
pc9_72 <- CreateSeuratObject(counts = gsm_72, project = "PC9_72")
pc9_72[["percent.mt"]] <- PercentageFeatureSet(pc9_72, pattern = "^MT-")
pc9_72 <- subset(pc9_72, subset = nFeature_RNA > 100 & nFeature_RNA < 5000 & percent.mt < 10)
pc9_72 <- NormalizeData(pc9_72, assay = "RNA")
pc9_72 <- FindVariableFeatures(pc9_72, assay = "RNA")
pc9_72 <- ScaleData(pc9_72, assay = "RNA")

# Integration pipeline
combined <- list(pc9_u, pc9_24, pc9_48, pc9_72)

# normalize and identify variable features for each dataset independently
combined <- lapply(X = combined, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst")
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = combined)

## Perform integration

radiation.anchors <- FindIntegrationAnchors(object.list = combined, anchor.features = features)
# this command creates an 'integrated' data assay
radiation.combined <- IntegrateData(anchorset = radiation.anchors)

## Perform an integrated analysis
DefaultAssay(radiation.combined) <- "RNA"
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
radiation.combined[["percent.mt"]] <- PercentageFeatureSet(ctrl, pattern = "^mt-")

# Visualize QC metrics as a violin plot
VlnPlot(radiation.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

#plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(ctrl, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#plot1 + plot2

ctrl <- subset(ctrl, subset = nFeature_RNA > 100 & nFeature_RNA < 1500)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(radiation.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
radiation.combined <- ScaleData(radiation.combined, verbose = FALSE)
radiation.combined <- RunPCA(radiation.combined, npcs = 30, verbose = FALSE)
radiation.combined <- RunUMAP(radiation.combined, reduction = "pca", dims = 1:30)
radiation.combined <- FindNeighbors(radiation.combined, reduction = "pca", dims = 1:30)
radiation.combined <- FindClusters(radiation.combined, resolution = 0.5)
# Visualization
ElbowPlot(radiation.combined)
DimPlot(radiation.combined, reduction = "umap")
p1 <- DimPlot(radiation.combined, reduction = "umap", group.by = "orig.ident")
p1_1 <- DimPlot(radiation.combined, reduction = "umap", split.by = "orig.ident", ncol = 2)
p3 <- DimPlot(radiation.combined, reduction = "umap", label = T)
ggsave(p1, file = "Original Identity Plot.pdf", width = 8, height = 7, limitsize = FALSE)
ggsave(p1_1, file = "Splitted Identity Plot.pdf", width =10, height = 10, limitsize = FALSE)
p2 <- DimPlot(radiation.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
ggsave(p3, file = "Original Identity Plot_clusters.jpg", width = 7, height = 7)
