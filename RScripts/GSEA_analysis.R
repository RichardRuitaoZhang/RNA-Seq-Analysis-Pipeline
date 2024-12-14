# GSEA analysis
GSEA_analysis <- function(deg_df, output_dir, org_db = org.Hs.eg.db, species = "hsa") {
  # Ensure the output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Ensure 'deg_df' contains the required columns: ENSEMBL and logFC
  if (!all(c("ENSEMBL", "log2FoldChange") %in% colnames(deg_df))) {
    stop("Input dataframe must contain 'ENSEMBL' and 'log2FoldChange' columns.")
  }
  
  # Extract gene list and logFC values
  gene_list <- deg_df$ENSEMBL
  log_fc <- deg_df$log2FoldChange
  
  # Convert ENSEMBL IDs to ENTREZ IDs
  gene_entrez <- bitr(gene_list, fromType = "ENSEMBL", 
                      toType = "ENTREZID", 
                      OrgDb = org_db)
  
  # Merge logFC values with ENTREZ IDs
  gene_df <- merge(deg_df, gene_entrez, by.x = "ENSEMBL", by.y = "ENSEMBL")
  
  # Ensure conversion worked
  if (nrow(gene_entrez) == 0) {
    stop("No valid ENTREZ IDs were found for the provided gene list.")
  }
  
  # Prepare ranked list for GSEA
  ranked_list <- gene_df$log2FoldChange
  names(ranked_list) <- gene_df$ENTREZID
  ranked_list <- sort(ranked_list, decreasing = TRUE)
  
  # Perform GSEA
  message("Performing GSEA analysis...")
  
  gsea_result <- gseKEGG(geneList     = ranked_list,
                         organism     = species,
                         pAdjustMethod = "BH",
                         minGSSize    = 10,
                         maxGSSize    = 500,
                         pvalueCutoff = 0.05)
  
  # Check if there are any significant pathways
  if (is.null(gsea_result) || nrow(gsea_result@result) == 0) {
    warning("No significant pathways were identified by GSEA.")
  } else {
    # Save the GSEA results to a CSV file
    output_file <- file.path(output_dir, "GSEA_results.csv")
    write.csv(as.data.frame(gsea_result), output_file, row.names = FALSE)
    message("GSEA analysis completed. Results saved in: ", output_file)
  }
  
  # Plot enrichment results
  if (!is.null(gsea_result) && nrow(gsea_result@result) > 0) {
    gsea_plot <- dotplot(gsea_result, showCategory = 20) +
      scale_color_distiller(palette = "RdBu", name = "Adjusted p-value") +
      theme_bw()
    plot_file <- file.path(output_dir, "GSEA_dotplot.pdf")
    ggsave(plot_file, plot = gsea_plot)
    message("GSEA dotplot saved in: ", plot_file)
  }
  
  return(gsea_result)
}

# Save the function to a file
save(GSEA_analysis, file = "GSEA_analysis.RData")
