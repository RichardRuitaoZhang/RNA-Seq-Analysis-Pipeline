# Load libraries
library(clusterProfiler)
library(org.Hs.eg.db)  # Use org.Mm.eg.db for mouse, etc.
library(enrichplot)
library(DOSE)
library(GO.db)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)  # Fo

# Assume `normalized_counts` is your normalized expression matrix
# Rows are genes, columns are samples
# `de_genes` is the list of significant DEGs
# read in counts dataset as a matrix
# set working dir
setwd("/Users/richardzhang/Desktop/Northwestern/Y1/BME_311/Final Project/data/")
cts <- read.csv("gene_expression.csv")

cts$ID <- sub("\\.\\d+$", "", cts$ID)


# Map ENSEMBL IDs to Gene Symbols
gene_annotations <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = cts$ID,
  columns = c("SYMBOL"),      # Retrieve gene symbols
  keytype = "ENSEMBL"         # Input IDs are ENSEMBL
)

cts$ENSEMBL <- cts$ID
cts <- merge(cts, gene_annotations, by = "ENSEMBL")


# 加载数据
# expression_data <- read.csv("your_expression_data.csv", row.names = 1)

expression_data <- cts[,c(3:11)]
rownames(expression_data) <- make.names(cts$SYMBOL, unique=TRUE)

# Assign unique row names
rownames(expression_data) <- expression_data$SYMBOL
expression_data <- expression_data[ , -1]  # Remove SYMBOL column

# Remove rows where rownames contain "NA"
expression_data <- expression_data[!grepl("NA", rownames(expression_data)), ]

# 提取样本数和基因数
num_genes <- nrow(expression_data)
num_samples <- ncol(expression_data)

# 创建 GCT 文件头信息
header <- c("#1.2", paste(num_genes, num_samples))

# Add column names (NAME, DESCRIPTION, and colnames of expression_data)
column_names <- c("NAME", "DESCRIPTION", colnames(expression_data))

# Rename colnames to Sample1, Sample2, ..., Sample9
colnames(expression_data) <- paste0("Sample", seq_len(ncol(expression_data)))

# Prepare the header
header <- "#1.2"  # GCT format header
data_info <- paste(nrow(expression_data), ncol(expression_data), sep = "\t")  # Row and column count

# Add column names (NAME, DESCRIPTION, and colnames of expression_data)
column_names <- c("NAME", "DESCRIPTION", colnames(expression_data))
data_info <- paste(nrow(expression_data), ncol(expression_data), sep = "\t")  # Row and column count
# Write to GCT file
writeLines(
  c(
    header,
    data_info = data_info,
    paste(column_names, collapse = "\t"),  # Add column names
    paste(
      rownames(expression_data),           # NAME
      rep("na", nrow(expression_data)),    # DESCRIPTION
      apply(expression_data, 1, paste, collapse = "\t"),  # Expression data
      sep = "\t"
    )
  ),
  "expression_data.gct"
)

### SKMEL28
# set working dir
setwd("/Users/richardzhang/Desktop/Northwestern/Y1/BME_311/Final Project/data/SKMEL28/")
cts <- read.csv("skmel28_expression_matrix.csv")

cts$ID <- sub("\\.\\d+$", "", cts$X)


# Map ENSEMBL IDs to Gene Symbols
gene_annotations <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = cts$ID,
  columns = c("SYMBOL"),      # Retrieve gene symbols
  keytype = "ENSEMBL"         # Input IDs are ENSEMBL
)

cts$ENSEMBL <- cts$ID
cts <- merge(cts, gene_annotations, by = "ENSEMBL")


# 加载数据
# expression_data <- read.csv("your_expression_data.csv", row.names = 1)

expression_data <- cts[,c(3:11)]
rownames(expression_data) <- make.names(cts$SYMBOL, unique=TRUE)

# Assign unique row names
rownames(expression_data) <- expression_data$SYMBOL
#expression_data <- expression_data[ , -1]  # Remove SYMBOL column

# Remove rows where rownames contain "NA"
expression_data <- expression_data[!grepl("NA", rownames(expression_data)), ]

# 提取样本数和基因数
num_genes <- nrow(expression_data)
num_samples <- ncol(expression_data)

# 创建 GCT 文件头信息
header <- c("#1.2", paste(num_genes, num_samples))

# Add column names (NAME, DESCRIPTION, and colnames of expression_data)
column_names <- c("NAME", "DESCRIPTION", colnames(expression_data))

# Rename colnames to Sample1, Sample2, ..., Sample9
colnames(expression_data) <- paste0("Sample", seq_len(ncol(expression_data)))

# Prepare the header
header <- "#1.2"  # GCT format header
data_info <- paste(nrow(expression_data), ncol(expression_data), sep = "\t")  # Row and column count

# Add column names (NAME, DESCRIPTION, and colnames of expression_data)
column_names <- c("NAME", "DESCRIPTION", colnames(expression_data))
data_info <- paste(nrow(expression_data), ncol(expression_data), sep = "\t")  # Row and column count
# Write to GCT file
writeLines(
  c(
    header,
    data_info = data_info,
    paste(column_names, collapse = "\t"),  # Add column names
    paste(
      rownames(expression_data),           # NAME
      rep("na", nrow(expression_data)),    # DESCRIPTION
      apply(expression_data, 1, paste, collapse = "\t"),  # Expression data
      sep = "\t"
    )
  ),
  "expression_data_SKMEL28.gct"
)
