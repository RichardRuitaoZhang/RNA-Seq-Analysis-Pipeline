# Title: DESeq2 Analysis for PC9 Cell Line Groups
# Author: Group 7
# Date: 2024-11-25
# Project: Final Project for BMD_ENG_311
# Description:
#     This script performs differential expression analysis on RNA-Seq data
#     for the PC9 cell line groups using the DESeq2 package. It focuses only
#     on samples 1–9, which include untreated, 3-day treated, and 9-day treated
#     groups, each in triplicate.
#     
#     The analysis includes:
#         - Reading the expression matrix for samples 1–9.
#         - Generating a colData metadata file specific to these samples.
#         - Setting up a DESeqDataSet object with the count matrix and colData.
#         - Performing differential expression analysis between conditions.
#         - Visualizing results with MA plots and heatmaps.


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

### 1.Differential expression analysis
# generate metadata for multi-comparison
metadata <- data.frame(Untreated = coldata$file_id[1:3],
                       Day3 = coldata$file_id[4:6],
                       Day9 = coldata$file_id[7:9])

# run multi-comparison
# load multi-group function
load("/Users/richardzhang/Desktop/Northwestern/Y1/BME_311/Final Project/data/DESeq_multi_group.RData")

multi <- DESeq_multi_group(exprSet = cts, meta = metadata, coldata = coldata)
dds_list <- multi$dds_list
results_list <- multi$results_list
resLFC_list <- multi$resLFC_list


### 2. Exploring and exporting results
# Names of the datasets for saving the plots
names_list <- c("PC9: Gefitinib-Treated Day 3 vs Untreated", "PC9: Gefitinib-Treated Day 9 vs Untreated", "PC9: Gefitinib-Treated Day 3 vs PC9: Gefitinib-Treated Day 9") # Update these names accordingly

load("/Users/richardzhang/Desktop/Northwestern/Y1/BME_311/Final Project/data/MA.RData")

resLFC_list <- MA(dds_list, names_list = names_list)

### 3.
load("/Users/richardzhang/Desktop/Northwestern/Y1/BME_311/Final Project/data/trans_vis.RData")

  
for (i in seq_along(results_list)) {
  trans_vis(names(results_list)[i], dds_list[[i]], coldata = coldata, output_dir="Visualization")
}


### 4. volcano plot
load("/Users/richardzhang/Desktop/Northwestern/Y1/BME_311/Final Project/data/volcano_plot.RData")
# plot
gene_panel <- c(
  "STAT3", "JAK1", "BCL6",       # STAT3 pathway
  "IGF1R", "IGFBP3", "IGFBP5",   # IGF pathway
  "PLCG1", "GNAS", "AXL",        # PLC/PKC and resistance-related genes
  "GSTM1", "GSTM3", "GPX4",      # Glutathione-related and oxidative stress response
  "ULK1", "GABARAPL1", "LC3A",   # Autophagy-related genes
  "CDKN1A", "POLQ", "RAD18",     # DNA repair and senescence
  "YAP", "BCL2"                  # YAP and anti-apoptotic regulation
)

gene_panel_whole <- c(
  # STAT3 Pathway
  "STAT3", "JAK1", "BCL6",
  
  # IGF Pathway
  "IGF1R", "IGFBP3", "IGFBP5", "IGFBP7",
  
  # PLC/PKC Pathway
  "ADCY6", "GNAL", "GNAS", "PDE1C", "PLCG1",
  
  # YAP Pathway
  "YAP", "PAI-1",
  
  # mTOR Pathway
  "RPS6", "EIF4G1", "EEF2K",
  
  # Resistance-Related Genes
  "AXL", "GAS6", "GPX4",
  
  # Glutathione-Related Genes
  "GSTA4", "GSTM1", "GSTM2", "GSTM3", "GSTM4",
  
  # Autophagy-Related Genes
  "ULK1", "GABARAPL1", "LC3A", "LC3B",
  
  # Senescence-Related Genes
  "CDKN1A", "CDKN2B",
  
  # DNA Repair Genes and Low-Fidelity Polymerases
  "POLH", "POLI", "POLK", "POLQ", "RAD18",
  
  # PI3K-AKT-mTOR Pathway
  "PIK3CA", "PIK3CB", "PIK3CD", "PIK3R1", "AKT1", "AKT2", "AKT3",
  "MTOR", "RPTOR", "RICTOR", "TSC1", "TSC2", "PTEN", "FOXO1", "FOXO3", "FOXO4"
)


for (i in 1:length(results_list)) {
  # Generate the volcano plot
  vcn <- volcano_plot(results_list[[i]], 
                      title = names_list[i], 
                      oncogenes = gene_panel_whole, 
                      padj_threshold = 0.05, 
                      lfc_threshold = 1)
  
  # Save each plot as a PDF
  pdf_file <- file.path("volcano_results", paste0(names_list[i], "_volcano_plot.pdf"))
  ggsave(filename = pdf_file, plot = vcn, width = 8, height = 8, device = cairo_pdf)
}

# 5. get DEGs
load("/Users/richardzhang/Desktop/Northwestern/Y1/BME_311/Final Project/data/get_exp.RData")
# result
exp_matrix <- list()

# Iterate through the list
for (dataset_name in names(results_list)) {
  # Get the DESeq2 result for the current dataset
  res <- results_list[[dataset_name]]
  
  # Generate a custom prefix for the output
  output_prefix <- paste0(dataset_name, "_processed")
  
  # Process the current result using your function and save it in the list
  exp_matrix[[dataset_name]] <- get_exp(res, output_prefix = output_prefix)
}

load("/Users/richardzhang/Desktop/Northwestern/Y1/BME_311/Final Project/data/get_DEGs.RData")
# Apply the save_DEGs function to each dataset in the list
for (sample in names(exp_matrix)) {
  res <- exp_matrix[[sample]]
  
  # Generate custom output prefix based on the dataset name
  output_prefix <- paste0(sample, "_deg_05")
  
  # Run the save_DEGs function
  get_DEGs(res, output_dir = output_prefix)
}
