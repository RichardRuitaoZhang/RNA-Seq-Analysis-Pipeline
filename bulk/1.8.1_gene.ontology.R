# Load libraries
library(clusterProfiler)
library(org.Hs.eg.db)  # Use org.Mm.eg.db for mouse, etc.
library(enrichplot)
library(DOSE)
library(GO.db)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)  # For combining multiple plots
library("RColorBrewer") # provides several color palette

#### Example: A list of gene symbols
gene_symbols <- c("TP53", "BRCA1", "EGFR", "MYC", "CDK2")

# Convert gene symbols to Entrez IDs
gene_list <- bitr(gene_symbols, fromType = "SYMBOL", 
                  toType = "ENTREZID", 
                  OrgDb = org.Hs.eg.db)

# For ranked gene lists (e.g., logFC values), you can prepare:
ranked_gene_list <- c(TP53 = 2.5, BRCA1 = 1.8, EGFR = -1.2, MYC = 0.5, CDK2 = -2.0)
ranked_gene_list <- sort(ranked_gene_list, decreasing = TRUE)

go_results <- enrichGO(gene          = gene_list$ENTREZID,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = "ENTREZID",
                       ont           = "BP",  # Biological Process; can be "MF" or "CC"
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.2)

# View results
head(go_results)

# Visualize results
dotplot(go_results, showCategory = 20)  # Show top 20 categories

# set working dir
setwd("/Users/richardzhang/Desktop/Northwestern/Y1/BME_311/Final Project/data/PC9/")

# read in DE data
d3_u_up <- read.csv("Day3_vs_Untreated_deg_05/deg_05_up.csv")
d3_u_down <- read.csv("Day3_vs_Untreated_deg_05/deg_05_down.csv")
d9_u_up <- read.csv("Day9_vs_Untreated_deg_05/deg_05_up.csv")
d9_u_down <- read.csv("Day9_vs_Untreated_deg_05/deg_05_down.csv")
d9_d3_up <- read.csv("Day9_vs_Day3_deg_05/deg_05_up.csv")
d9_d3_down <- read.csv("Day9_vs_Day3_deg_05/deg_05_down.csv")

load("/Users/richardzhang/Desktop/Northwestern/Y1/BME_311/Final Project/data/GO_analysis.RData")
d3_u_up_GO <- GO_analysis(d3_u_up, output_dir = "d3_u_up")
d3_u_down_GO <- GO_analysis(d3_u_down, output_dir = "d3_u_down")
d9_u_up_GO <- GO_analysis(d9_u_up, output_dir = "d9_u_up")
d9_u_down_GO <- GO_analysis(d9_u_down, output_dir = "d9_u_down")
d9_d3_up_GO <- GO_analysis(d9_d3_up, output_dir = "d9_d3_up")
d9_d3_down_GO <- GO_analysis(d9_d3_down, output_dir = "d9_d3_down")

load("/Users/richardzhang/Desktop/Northwestern/Y1/BME_311/Final Project/data/combined_bubble_plot2.RData")
# plot the bubble
bubble(as.data.frame(d3_u_up_GO$BP), title = "Day 3 vs Untreated", top_n = 10)
combined_bubble_plot(d3_u_up_GO)

# combine bubble plots
combined_bubble_plot2(d3_u_up_GO, title = "Day 3 vs Untreated Up-Regulated")
combined_bubble_plot2(d3_u_down_GO, title = "Day 3 vs Untreated Down-Regulated")
combined_bubble_plot2(d9_u_up_GO, title = "Day 9 vs Untreated Up-Regulated")
combined_bubble_plot2(d9_u_down_GO, title = "Day 9 vs Untreated Down-Regulated")
combined_bubble_plot2(d9_d3_up_GO, title = "Day 9 vs Day 3 Up-Regulated")
combined_bubble_plot2(d9_d3_down_GO, title = "Day 9 vs Day 3 Down-Regulated")

data <- data %>%
  mutate(GeneRatio = as.numeric(sub("/.*", "", GeneRatio)) / as.numeric(sub(".*/", "", GeneRatio)))

### SKMEL28
# set working dir
setwd("/Users/richardzhang/Desktop/Northwestern/Y1/BME_311/Final Project/data/SKMEL28/")

# read in DE data
d3_u_up <- read.csv("Day3_vs_Untreated_deg_05/deg_05_up.csv")
d3_u_down <- read.csv("Day3_vs_Untreated_deg_05/deg_05_down.csv")
d9_u_up <- read.csv("Day9_vs_Untreated_deg_05/deg_05_up.csv")
d9_u_down <- read.csv("Day9_vs_Untreated_deg_05/deg_05_down.csv")
d9_d3_up <- read.csv("Day9_vs_Day3_deg_05/deg_05_up.csv")
d9_d3_down <- read.csv("Day9_vs_Day3_deg_05/deg_05_down.csv")

load("/Users/richardzhang/Desktop/Northwestern/Y1/BME_311/Final Project/data/GO_analysis.RData")
d3_u_up_GO <- GO_analysis(d3_u_up, output_dir = "d3_u_up")
d3_u_down_GO <- GO_analysis(d3_u_down, output_dir = "d3_u_down")
d9_u_up_GO <- GO_analysis(d9_u_up, output_dir = "d9_u_up")
d9_u_down_GO <- GO_analysis(d9_u_down, output_dir = "d9_u_down")
d9_d3_up_GO <- GO_analysis(d9_d3_up, output_dir = "d9_d3_up")
d9_d3_down_GO <- GO_analysis(d9_d3_down, output_dir = "d9_d3_down")

load("/Users/richardzhang/Desktop/Northwestern/Y1/BME_311/Final Project/data/combined_bubble_plot2.RData")
# plot the bubble
bubble(as.data.frame(d3_u_up_GO$BP), title = "Day 3 vs Untreated", top_n = 10)
combined_bubble_plot(d3_u_up_GO)

# combine bubble plots
combined_bubble_plot2(d3_u_up_GO, title = "Day 3 vs Untreated Up-Regulated")
combined_bubble_plot2(d3_u_down_GO, title = "Day 3 vs Untreated Down-Regulated")
combined_bubble_plot2(d9_u_up_GO, title = "Day 9 vs Untreated Up-Regulated")
combined_bubble_plot2(d9_u_down_GO, title = "Day 9 vs Untreated Down-Regulated")
combined_bubble_plot2(d9_d3_up_GO, title = "Day 9 vs Day 3 Up-Regulated")
combined_bubble_plot2(d9_d3_down_GO, title = "Day 9 vs Day 3 Down-Regulated")
