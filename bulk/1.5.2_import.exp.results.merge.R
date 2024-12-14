"
Project: Final Project for BMD_ENG_311
Author: Group 7
Date: 2024-11-21
"

# set working dir
setwd("/projects/e30836/project/group7/bulk/expmtx/")

dir <- "/projects/e30836/project/group7/bulk/expression/"

# collect annotation dataframe manually, named as 'metadata_pc9.csv'
# read in
samples_pc9 <- read.table('metadata_pc9.csv')

# set file path
files_pc9 <- file.path(dir, samples_pc9$run, paste0(samples_pc9$run, ".genes.results.gz"))
names(files_pc9) <- c(
  paste0("untreated", 1:3),       # First 3 names: untreated1, untreated2, untreated3
  paste0("day 3_", 1:3),          # Next 3 names: day 3_1, day 3_2, day 3_3
  paste0("day 9_", 1:3)           # Last 3 names: day 9_1, day 9_2, day 9_3
)
txi.rsem <- tximport(files_pc9, type = "rsem", txIn = FALSE, txOut = FALSE)
head(txi.rsem$counts)
