library(pheatmap)
library(RColorBrewer) # provides several color palette

# Assume `normalized_counts` is your normalized expression matrix
# Rows are genes, columns are samples
# `de_genes` is the list of significant DEGs
# read in counts dataset as a matrix
# set working dir
setwd("/Users/richardzhang/Desktop/Northwestern/Y1/BME_311/Final Project/data/PC9/")
cts <- as.matrix(read.csv("PC9_expression_matrix.csv", row.names = 1))
d9_u <- read.csv("Day9_vs_Untreated_deg_05/deg_05.csv")

# Create a mapping of ENSEMBL to SYMBOL from the d9_u data
ensembl_to_symbol <- setNames(d9_u$SYMBOL, d9_u$ENSEMBL)

# Subset the matrix for DEGs and map row names to SYMBOL
rownames(cts) <- sub("\\.\\d+$", "", rownames(cts))
de_matrix <- cts[rownames(cts) %in% names(ensembl_to_symbol), ]

# Update row names to SYMBOL
rownames(de_matrix) <- ensembl_to_symbol[rownames(de_matrix)]

# Proceed with the rest of the analysis (e.g., log2-transform, scale, and heatmap)
de_matrix <- log2(de_matrix + 1)
scaled_matrix <- t(scale(t(de_matrix)))

# Define output dimensions for heatmap
num_genes <- nrow(scaled_matrix)
num_samples <- ncol(scaled_matrix)
heatmap_height <- num_genes * 0.15
heatmap_width <- num_samples * 1.5

# Save the heatmap to a PDF file
pdf("heatmap_pc9_long.pdf", width = heatmap_width, height = heatmap_height)
pheatmap(
  mat = scaled_matrix,
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "PiYG")))(100),  # Color gradient
  cluster_rows = TRUE,  # Cluster genes
  cluster_cols = TRUE,  # Cluster samples
  show_rownames = TRUE, # Show gene names
  show_colnames = TRUE  # Show sample names
)
dev.off()



# Set working directory
setwd("/Users/richardzhang/Desktop/Northwestern/Y1/BME_311/Final Project/data/SKMEL28/")

# Load expression matrix
cts <- as.matrix(read.csv("skmel28_expression_matrix.csv", row.names = 1))
d9_u <- read.csv("Day9_vs_Untreated_deg_05/deg_05.csv")

# Remove version numbers from ENSEMBL IDs
rownames(cts) <- sub("\\.\\d+$", "", rownames(cts))

# Create a mapping of ENSEMBL to SYMBOL from the d9_u data
ensembl_to_symbol <- setNames(d9_u$SYMBOL, d9_u$ENSEMBL)
# Subset `de_matrix` (assumes `d9_u` SYMBOL mapping is already handled)
de_matrix <- cts[rownames(cts) %in% names(ensembl_to_symbol), ]

# Map ENSEMBL to SYMBOL for row names
rownames(de_matrix) <- ensembl_to_symbol[rownames(de_matrix)]

# Log2-transform the matrix (to stabilize variance)
de_matrix <- log2(de_matrix + 1)

# Scale rows (genes) to mean=0 and standard deviation=1
scaled_matrix <- t(scale(t(de_matrix)))

# Define heatmap dimensions
num_genes <- nrow(scaled_matrix)
num_samples <- ncol(scaled_matrix)
heatmap_height <- num_genes * 0.15
heatmap_width <- num_samples * 1.5

# Save the heatmap as PDF
pdf("heatmap_skmel28_long.pdf", width = heatmap_width, height = heatmap_height)

# Generate the heatmap
pheatmap(
  mat = scaled_matrix,
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "PiYG")))(100),  # Color gradient
  cluster_rows = TRUE,  # Cluster genes
  cluster_cols = TRUE,  # Cluster samples
  show_rownames = TRUE, # Show gene names
  show_colnames = TRUE  # Show sample names
)

# Close the PDF device
dev.off()
print("PDF file generated successfully.")
