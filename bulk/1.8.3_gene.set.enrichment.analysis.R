# Load necessary libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(ggplot2)
library(enrichplot)
library(viridis)
library(RColorBrewer) # provides several color palette
library(msigdbr)
library(fgsea)
load("/Users/richardzhang/Desktop/Northwestern/Y1/BME_311/Final Project/data/GSEA_analysis.RData")

# set working dir
setwd("/Users/richardzhang/Desktop/Northwestern/Y1/BME_311/Final Project/data/")

# read in DE data
d3_u_up <- read.csv("Day3_vs_Untreated_deg_05/deg_05_up.csv")
d3_u_down <- read.csv("Day3_vs_Untreated_deg_05/deg_05_down.csv")
d9_u_up <- read.csv("Day9_vs_Untreated_deg_05/deg_05_up.csv")
d9_u_down <- read.csv("Day9_vs_Untreated_deg_05/deg_05_down.csv")
d9_d3_up <- read.csv("Day9_vs_Day3_deg_05/deg_05_up.csv")
d9_d3_down <- read.csv("Day9_vs_Day3_deg_05/deg_05_down.csv")

d3_u_up_GSEA <- GSEA_analysis(d3_u_up, output_dir = "d3_u_up")
d3_u_down_GSEA <- GSEA_analysis(d3_u_down, output_dir = "d3_u_down")
d9_u_up_GSEA <- GSEA_analysis(d9_u_up, output_dir = "d9_u_up")
d9_u_down_GSEA <- GSEA_analysis(d9_u_down, output_dir = "d9_u_down")
d9_d3_up_GSEA <- GSEA_analysis(d9_d3_up, output_dir = "d9_d3_up")
d9_d3_down_GSEA <- GSEA_analysis(d9_d3_down, output_dir = "d9_d3_down")


# Function to draw enrichment plots for multiple gene sets
draw_multiple_enrichment_plots <- function(gsea_result, gene_sets, output_dir) {
  # Ensure output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Loop through each gene set and generate a plot
  for (gene_set in gene_sets) {
    # Check if the gene set exists in the results
    if (gene_set %in% gsea_result@result$Description) {
      # Generate enrichment plot
      enrichment_plot <- gseaplot2(
        gsea_result,
        geneSetID = gene_set,
        title = paste("Enrichment Plot for", gene_set),
        color = "blue"
      )
      
      # Save the plot
      output_file <- file.path(output_dir, paste0("enrichment_plot_", gene_set, ".pdf"))
      ggsave(output_file, plot = enrichment_plot, width = 8, height = 6)
      message("Enrichment plot saved for ", gene_set, " at: ", output_file)
    } else {
      warning("Gene set '", gene_set, "' not found in GSEA results.")
    }
  }
}


# Define your gene sets of interest
gene_sets <- c("IGFR1", "PLC-PKC", "YAP", "STAT3", "mTOR")

# Specify the output directory
output_dir <- "enrichment_plots"

# Call the function with your GSEA result and gene sets
draw_multiple_enrichment_plots(gsea_result = d9_u_down_GSEA, 
                               gene_sets = gene_sets, 
                               output_dir = output_dir)

gseaplot2(
  d3_u_up_GSEA,
  geneSetID = "hsa04150",
  title = "mTOR signaling pathway"
)

gseaplot2(
  d3_u_up_GSEA,
  geneSetID = "hsa04150",
  title = "mTOR signaling pathway"
)

# fgesa
# Example: Get human (Homo sapiens) Hallmark pathways
pathways_h_df <- msigdbr(species = "Homo sapiens", category = "H")  # 'H' is for Hallmark gene sets

# Convert to a list format for fgsea
pathways_h <- split(pathways_h_df$gene_symbol, pathways_h_df$gs_name)

pathways_c2_df <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
pathways_c2 <- split(pathways_c2_df$gene_symbol, pathways_c2_df$gs_name)

# Retrieve Reactome pathways for humans
reactome_pathways <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
# Convert to a list of pathways
reactome_list <- split(reactome_pathways$gene_symbol, reactome_pathways$gs_name)

fgsea_d9_u_up <- fgsea_analysis(
  deg_df = d9_u_up,
  output_dir = "fgsea_results",
  pathways = reactome_list,
  minSize = 5,  # Adjust minimum gene set size as needed
  maxSize = 500,  # Adjust maximum gene set size as needed
  nperm = 1000  # Number of permutations
)

pathway_genes <- fgsea_d9_u_up[fgsea_d9_u_up$pathway == "REACTOME_APOPTOSIS", "leadingEdge"]
# Extract the first (and only) element from the list
pathway_genes <- unlist(pathway_genes)

# Check the content of pathway_genes
print(pathway_genes)

# Plot enrichment for a specific pathway
plotEnrichment(
  pathway = pathway_genes,
  stats = ranked_list
) + ggtitle("Enrichment Plot for REACTOME_APOPTOSIS
")

# Extract gene list and log2FoldChange values
ranked_list <- d9_u_up$log2FoldChange
names(ranked_list) <- d9_u_up$SYMBOL
ranked_list <- sort(ranked_list, decreasing = TRUE)
ranked_list <- jitter(ranked_list)

pathway <- unname(pathway)
overlap <- intersect(pathway, names(ranked_list))
length(overlap)  # Should be greater than 0

# Create a custom enrichment plot
pathway <- pathway_genes
pathway <- unname(pathway)
enrichment_scores <- calcGseaStat(
  stats = ranked_list,
  selectedStats = pathway
)
gene_ranks <- seq_along(ranked_list)

ggplot() +
  geom_line(aes(x = ranks, y = es), color = "green") +
  geom_vline(xintercept = which(names(ranks) %in% pathway), color = "blue", alpha = 0.5) +
  labs(
    title = "Custom Enrichment Plot",
    x = "Rank in Ordered Dataset",
    y = "Enrichment Score"
  ) +
  theme_minimal()


gseaplot2(
  stats = ranked_list,
  pathway = pathways[[top_pathway]],
  title = top_pathway,
  NES = fgsea_results$NES[1]
)