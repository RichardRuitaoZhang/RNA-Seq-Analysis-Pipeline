# Title: Gene Ontology Analysis Using DAVID Outputs
# Author: [Your Name or Team Name]
# Date: 2024-12-04
# Project: Gene Ontology Analysis for Differentially Expressed Genes
# Description:
#     This script automates the processing of gene ontology data obtained from DAVID.
#     It is designed to:
#         - Identify and organize up- and down-regulated gene lists for GO categories
#           (Biological Process, Cellular Component, Molecular Function).
#         - Read DAVID output files and load them into R as structured data.
#         - Print the top 5 terms for quick inspection from each category.
#         
#     The analysis includes:
#         - Reading GO terms from files with `_up.txt` and `_down.txt` suffixes.
#         - Organizing the terms into lists based on regulation status (up/down).
#         - Displaying summary results for quick inspection.


# Load Required Packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)

# Get lists of file names for up- and down-regulated terms
input_dir <- "/Users/richardzhang/Desktop/Northwestern/Y1/BME_311/Homework/Homework 3/"

# Get lists of file names for up- and down-regulated terms
up_files <- list.files(path = input_dir, pattern = "_up\\.txt$", full.names = TRUE)
down_files <- list.files(path = input_dir, pattern = "_down\\.txt$", full.names = TRUE)

# Read "up" files into a list with file names as identifiers
up_data <- lapply(up_files, function(file) {
  read.delim(file)
})

# Read "down" files into a separate list
down_data <- lapply(down_files, function(file) {
  read.delim(file)
})

# Assign names to the lists based on file names
names(up_data) <- sub("\\.txt$", "", up_files)
names(down_data) <- sub("\\.txt$", "", down_files)

# 2: Define a Function to Print Top 5 Rows
# Function to print the top 5 rows of each data frame
print_top5 <- function(data, title) {
  cat("\n", title, "\n")
  if (!is.null(data)) {
    print(head(data, 5))
  } else {
    cat("No data available for", title, "\n")
  }
}

# ==============================================================================
# Step 3: Display Top 5 Rows for Each GO Term Category
# ==============================================================================
# Up-regulated categories
print_top5(up_data$BP_up, "Top 5 Rows of Up-regulated: Biological Process")
print_top5(up_data$CC_up, "Top 5 Rows of Up-regulated: Cellular Component")
print_top5(up_data$MF_up, "Top 5 Rows of Up-regulated: Molecular Function")

# Down-regulated categories
print_top5(down_data$BP_down, "Top 5 Rows of Down-regulated: Biological Process")
print_top5(down_data$CC_down, "Top 5 Rows of Down-regulated: Cellular Component")
print_top5(down_data$MF_down, "Top 5 Rows of Down-regulated: Molecular Function")

category_list <- c("Biological Process", "Cellular Component", "Molecular Function")

# 4. Bubble Plots
data <- read.delim("/Users/richardzhang/Desktop/Northwestern/Y1/BME_311/Final Project/results/DAVID/chart_B97D0AA610CA1733440428480.txt")
data$logP <- -log10(data$PValue)
data <- data %>%
  separate(Term, into = c("GO_ID", "Term"), sep = "~")
data$Category <- rep(category_list[1], nrow(data))
data <- data[data$PValue < 0.05, ]

# Create the bubble plot
ggplot(data, aes(x = Fold.Enrichment, y = reorder(Term, Fold.Enrichment), size = Count, color = logP)) +
  geom_point(alpha = 0.7) +  # Add transparency
  scale_size(range = c(3, 10)) +  # Adjust bubble sizes
  scale_color_gradientn(colors = brewer.pal(9, "RdYlBu")) +
  theme_bw()+
  labs(
    title = "Bubble Plot of GO Terms",
    x = "Fold Enrichment",
    y = "GO Term",
    color = "-log10(p-value)",
    size = "Gene Count"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10))  # Adjust axis text size for clarity


ggplot(data, aes(x = reorder(Term, -Count), y = Count, fill = -log10(PValue))) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(
    title = "Bar Plot of GO Terms",
    x = "GO Term",
    y = "Gene Count",
    fill = "-log10(p-value)"
  ) +
  theme_minimal()

ggplot(data, aes(x = Fold.Enrichment, y = reorder(Term, Fold.Enrichment), color = -log10(PValue), size = Count)) +
  geom_point(alpha = 0.7) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(
    title = "Dot Plot of GO Terms",
    x = "Fold Enrichment",
    y = "GO Term",
    color = "-log10(p-value)",
    size = "Gene Count"
  ) +
  theme_minimal()

library(treemap)

treemap(
  data,
  index = "Term",
  vSize = "Count",
  vColor = "logP",
  type = "value",
  title = "Treemap of GO Terms",
  palette = "RdYlBu"
)
