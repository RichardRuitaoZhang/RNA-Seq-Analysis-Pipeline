# Load libraries
library(clusterProfiler)
library(org.Hs.eg.db)  # Use org.Mm.eg.db for mouse, etc.
library(enrichplot)
library(DOSE)
library(KEGG.db)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)  # For combining multiple plots
library(RColorBrewer) # provides several color palette
# load functions
load("/Users/richardzhang/Desktop/Northwestern/Y1/BME_311/Final Project/data/KEGG_analysis.RData")
load("/Users/richardzhang/Desktop/Northwestern/Y1/BME_311/Final Project/data/bar.RData")

# set working dir of PC9
setwd("/Users/richardzhang/Desktop/Northwestern/Y1/BME_311/Final Project/data/PC9/")

# read in DE data
d3_u_up <- read.csv("Day3_vs_Untreated_deg_05/deg_05_up.csv")
d3_u_down <- read.csv("Day3_vs_Untreated_deg_05/deg_05_down.csv")
d9_u_up <- read.csv("Day9_vs_Untreated_deg_05/deg_05_up.csv")
d9_u_down <- read.csv("Day9_vs_Untreated_deg_05/deg_05_down.csv")
d9_d3_up <- read.csv("Day9_vs_Day3_deg_05/deg_05_up.csv")
d9_d3_down <- read.csv("Day9_vs_Day3_deg_05/deg_05_down.csv")

d3_u_up_KEGG <- KEGG_analysis(d3_u_up, output_dir = "d3_u_up")
d3_u_down_KEGG <- KEGG_analysis(d3_u_down, output_dir = "d3_u_down")
d9_u_up_KEGG <- KEGG_analysis(d9_u_up, output_dir = "d9_u_up")
d9_u_down_KEGG <- KEGG_analysis(d9_u_down, output_dir = "d9_u_down")
d9_d3_up_KEGG <- KEGG_analysis(d9_d3_up, output_dir = "d9_d3_up")
d9_d3_down_KEGG <- KEGG_analysis(d9_d3_down, output_dir = "d9_d3_down")

# bar plots
bar(d3_u_up_KEGG, output_file = "Day 3 vs Untreated Up-Regulated", output_dir = "KEGG")
bar(d3_u_down_KEGG, output_file = "Day 3 vs Untreated Down-Regulated", output_dir = "KEGG")
bar(d9_u_up_KEGG, output_file = "Day 9 vs Untreated Up-Regulated", output_dir = "KEGG")
bar(d9_u_down_KEGG, output_file = "Day 9 vs Untreated Down-Regulated", output_dir = "KEGG")
bar(d9_d3_up_KEGG, output_file = "Day 9 vs Day 3 Up-Regulated", output_dir = "KEGG")
bar(d9_d3_down_KEGG, output_file = "Day 9 vs Day 3 Down-Regulated", output_dir = "KEGG")


# set working dir of skmel28
setwd("/Users/richardzhang/Desktop/Northwestern/Y1/BME_311/Final Project/data/SKMEL28/")

# read in DE data
d3_u_up <- read.csv("Day3_vs_Untreated_deg_05/deg_05_up.csv")
d3_u_down <- read.csv("Day3_vs_Untreated_deg_05/deg_05_down.csv")
d9_u_up <- read.csv("Day9_vs_Untreated_deg_05/deg_05_up.csv")
d9_u_down <- read.csv("Day9_vs_Untreated_deg_05/deg_05_down.csv")
d9_d3_up <- read.csv("Day9_vs_Day3_deg_05/deg_05_up.csv")
d9_d3_down <- read.csv("Day9_vs_Day3_deg_05/deg_05_down.csv")

d3_u_up_KEGG <- KEGG_analysis(d3_u_up, output_dir = "d3_u_up")
d3_u_down_KEGG <- KEGG_analysis(d3_u_down, output_dir = "d3_u_down")
d9_u_up_KEGG <- KEGG_analysis(d9_u_up, output_dir = "d9_u_up")
d9_u_down_KEGG <- KEGG_analysis(d9_u_down, output_dir = "d9_u_down")
d9_d3_up_KEGG <- KEGG_analysis(d9_d3_up, output_dir = "d9_d3_up")
d9_d3_down_KEGG <- KEGG_analysis(d9_d3_down, output_dir = "d9_d3_down")

# bar plots
bar(d3_u_up_KEGG, output_file = "Day 3 vs Untreated Up-Regulated", output_dir = "KEGG")
bar(d3_u_down_KEGG, output_file = "Day 3 vs Untreated Down-Regulated", output_dir = "KEGG")
bar(d9_u_up_KEGG, output_file = "Day 9 vs Untreated Up-Regulated", output_dir = "KEGG")
bar(d9_u_down_KEGG, output_file = "Day 9 vs Untreated Down-Regulated", output_dir = "KEGG")
bar(d9_d3_up_KEGG, output_file = "Day 9 vs Day 3 Up-Regulated", output_dir = "KEGG")
bar(d9_d3_down_KEGG, output_file = "Day 9 vs Day 3 Down-Regulated", output_dir = "KEGG")
