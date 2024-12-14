# 3 conditions together

# Load Libraries
library("DESeq2") # package for DESeq run
library("combinat")
library("BiocParallel") # for parallel evaluation 
library("dplyr") # for data compilation
library("ggplot2") # plotting
library("vsn") #
library("pheatmap") # for heatmap plotting
library("RColorBrewer") # provides several color palette

# set working dir
setwd("/Users/richardzhang/Desktop/Northwestern/Y1/BME_311/Final Project/data/PC9/")

# read in counts dataset as a matrix
cts <- as.matrix(read.csv("PC9_expression_matrix.csv", row.names = 1))
# read in annotation dataset, and use the first column as rownames
coldata <- read.csv("colData_PC9.csv", row.names=1)
# subset 2 columns: condition and type
coldata <- coldata[,c("file_id", "condition", "read_type")]
# convert type of items in condition into 'factor' type.
coldata$condition <- factor(coldata$condition)
# convert type of items in type into 'factor' type.
coldata$read_type <- factor(coldata$read_type)
# remove the character 'fb' in rownames of coldata
rownames(coldata) <- coldata$file_id
# subset the cts by column names included in the row names of 'coldata', and reorder them in order of coldata.
cts <- cts[, rownames(coldata)]
cts <- round(cts)
# Create the DESeqDataSet
dds_all <- DESeqDataSetFromMatrix(
  countData = cts,
  colData = coldata,
  design = ~ condition
)
# Run the DESeq2 pipeline
dds_all <- DESeq(dds_all)

levels(dds_all$condition)
# Output: "day3" "day9" "untreated"

# Change the Reference Level (Optional)
dds_all$condition <- relevel(dds_all$condition, ref = "Untreated")
dds_all <- DESeq(dds_all)

# Perform variance stabilizing transformation
vsd <- vst(dds_all, blind = FALSE)

# 1. Generate PCA data
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)

# Extract percentage variance for labels
percentVar <- round(100 * attr(pca_data, "percentVar"))

# Create the PCA plot with a title
pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA Plot for PC9 Samples") +
  theme_bw()

# Save the plot to a file
ggsave("PCA_plot_PC9.png", plot = pca_plot, width = 8, height = 6, dpi = 300)


# 2. Sample-to-sample distances heatmap
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)

# Ensure the number of rows in coldata matches the number of samples in the matrix
coldata <- coldata[rownames(sampleDistMatrix), , drop = FALSE]
rownames(sampleDistMatrix) <- paste(coldata$condition, coldata$read_type, sep = "-")

colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

png(filename = file.path("sample_distance_heatmap_PC9.png"), 
    width = 3200, height = 3200, res = 300)
pheatmap(sampleDistMatrix, 
         clustering_distance_rows = sampleDists, 
         clustering_distance_cols = sampleDists, 
         col = colors, 
         main = paste("Sample Distance Heatmap - PC9"))
dev.off()


### SKMEL28
# set working dir
setwd("/Users/richardzhang/Desktop/Northwestern/Y1/BME_311/Final Project/data/SKMEL28/")

# read in counts dataset as a matrix
cts <- as.matrix(read.csv("SKMEL28_expression_matrix.csv", row.names = 1))
# read in annotation dataset, and use the first column as rownames
coldata <- read.csv("colData_SKMEL28.csv", row.names=1)
# subset 2 columns: condition and type
coldata <- coldata[,c("file_id", "condition", "read_type")]
# convert type of items in condition into 'factor' type.
coldata$condition <- factor(coldata$condition)
# convert type of items in type into 'factor' type.
coldata$read_type <- factor(coldata$read_type)
# remove the character 'fb' in rownames of coldata
rownames(coldata) <- coldata$file_id
# subset the cts by column names included in the row names of 'coldata', and reorder them in order of coldata.
cts <- cts[, rownames(coldata)]
cts <- round(cts)
# Create the DESeqDataSet
dds_all <- DESeqDataSetFromMatrix(
  countData = cts,
  colData = coldata,
  design = ~ condition
)
# Run the DESeq2 pipeline
dds_all <- DESeq(dds_all)

levels(dds_all$condition)
# Output: "day3" "day9" "untreated"

# Change the Reference Level (Optional)
dds_all$condition <- relevel(dds_all$condition, ref = "Untreated")
dds_all <- DESeq(dds_all)

# Perform variance stabilizing transformation
vsd <- vst(dds_all, blind = FALSE)

# 1. Generate PCA data
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)

# Extract percentage variance for labels
percentVar <- round(100 * attr(pca_data, "percentVar"))

# Create the PCA plot with a title
pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA Plot for SKMEL28 Samples") +
  theme_bw()

# Save the plot to a file
ggsave("PCA_plot_SKMEL28.png", plot = pca_plot, width = 8, height = 6, dpi = 300)

# 2. Sample-to-sample distances heatmap
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)

# Ensure the number of rows in coldata matches the number of samples in the matrix
coldata <- coldata[rownames(sampleDistMatrix), , drop = FALSE]
rownames(sampleDistMatrix) <- paste(coldata$condition, coldata$read_type, sep = "-")

colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

png(filename = file.path("sample_distance_heatmap_SKMEL28.png"), 
    width = 3200, height = 3200, res = 300)
pheatmap(sampleDistMatrix, 
         clustering_distance_rows = sampleDists, 
         clustering_distance_cols = sampleDists, 
         col = colors, 
         main = paste("Sample Distance Heatmap - SKMEL28"))
dev.off()
