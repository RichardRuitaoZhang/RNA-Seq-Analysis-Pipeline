# Title: Prepare RSEM Output for DESeq2 Using tximport
# Author: Group 7
# Date: 2024-11-25
# Project: Final Project for BMD_ENG_311
# Description:
#     This script uses the tximport package to organize gene-level expression
#     data generated by RSEM. The processed expression matrix is prepared for
#     downstream analysis in DESeq2.

# set working dir
setwd("/Users/richardzhang/Desktop/Northwestern/Y1/BME_311/Final Project/data/")

# read in coldata for PC9
samples_pc9 <- read.csv('colData.csv', nrows = 9)
samples_skmel28 <- read.csv('colData.csv', nrows = 9)

### 1. run this line when using GEO downloaded Supplementary file
# set file path
dir <- "/Users/richardzhang/Desktop/Northwestern/Y1/BME_311/Final Project/data/GSE162045_RAW/"
filename <- paste0(samples_pc9$running_id, "_PC9_", seq(1, 9), ".txt.gz")
# get file paths
files_pc9 <- file.path(dir, filename)
names(files_pc9) <- c(
  paste0("untreated_", 1:3),       # First 3 names: untreated1, untreated2, untreated3
  paste0("gefitinib-treated_day 3_", 1:3),          # Next 3 names: day 3_1, day 3_2, day 3_3
  paste0("gefitinib-treated_day 9_", 1:3)           # Last 3 names: day 9_1, day 9_2, day 9_3
)

# Initialize an empty list to store data frames
data_list <- lapply(files_pc9, function(file) {
  # Read each file
  data <- read.table(file, header = TRUE, stringsAsFactors = FALSE)
  # Rename the second column to match the sample name (use file name without extension)
  colnames(data)[2] <- gsub("\\.txt\\.gz$", "", basename(file))
  return(data)
})

# Merge all data frames by "ID"
merged_data <- Reduce(function(x, y) merge(x, y, by = "feature", all = TRUE), data_list)

# Set rownames to "ID" and drop the column
rownames(merged_data) <- merged_data$feature
count_matrix <- merged_data[, -1]

samples_pc9$file_id <- paste0(samples_pc9$running_id, "_PC9_", seq(1, 9))



### 2. run this line when using RSEM output expression data
merged_data <- read.table("PC9/pc9_expression_matrix.txt", sep = "\t", header = T)
# Set rownames to "ID" and drop the column
rownames(merged_data) <- merged_data$X
count_matrix <- merged_data[, -1]

# Create the colname vector
colname <- c(
  "PC9_Untreated_Sample_1", "PC9_Untreated_Sample_2", "PC9_Untreated_Sample_3",
  "PC9_Day_3_Sample_1", "PC9_Day_3_Sample_2", "PC9_Day_3_Sample_3",
  "PC9_Day_9_Sample_1", "PC9_Day_9_Sample_2", "PC9_Day_9_Sample_3"
)

colnames(count_matrix) <- colname

# add column IDs into coldata
samples_pc9$file_id <- colnames(count_matrix)



### 3. save the generated expression and coldata
write.csv(samples_pc9, "PC9/colData_PC9.csv", row.names = F)
# save expression matrix
write.csv(count_matrix, "PC9/PC9_expression_matrix.csv", row.names = T)

