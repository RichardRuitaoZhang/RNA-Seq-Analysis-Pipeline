# Title: Generate colData for DESeq2 Analysis
# Author: Group 7
# Date: 2024-11-25
# Project: Final Project for BMD_ENG_311
# Description:
#     This script generates a colData metadata file for RNA-Seq differential expression analysis using DESeq2.
#     It reads sample names from a provided list, assigns experimental annotations based on known conditions,
#     and saves the metadata in CSV format.
#     
#     Metadata includes:
#         - Sample names
#         - Cell type (CD9 or HE28)
#         - Condition (untreated, 3_day_treated, or 9_day_treated)
#         - Replicates (1, 2, or 3)
#
# Usage:
#     Place the sample names in a text file named 'SRR_Acc_List.txt'.
#     Update the file name in the script if needed.
#     Run the script in R to generate the 'colData.csv' file for DESeq2.

# set working dir
setwd("/Users/richardzhang/Desktop/Northwestern/Y1/BME_311/Final Project/data/")

# read sample names from the file
samples <- readLines("SRR_Acc_List.txt")

# refine metadata based on the sample information
running_id <- sprintf("GSM4932%03d", 167:184)
cell_type <- c(rep("PC9", 9), rep("SKMEL28", 9))
type <- rep("single-end", 18)
condition <- rep(c("Untreated", "Day3", "Day9"), each = 3)
replicate <- rep(1:3, times = 6)

# Combine into a data frame
colData <- data.frame(
  sample = samples[33:50], # First 18 samples
  running_id = running_id,
  read_type = type,
  cell_type = cell_type,
  condition = condition,
  replicate = replicate
)

# Print the metadata to check
print(colData)

# Save colData to a CSV file for DESeq2
write.csv(colData, "colData.csv", row.names = FALSE)

