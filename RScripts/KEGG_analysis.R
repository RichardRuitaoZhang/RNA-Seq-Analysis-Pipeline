KEGG_analysis <- function(deg_df, output_dir, org_db = org.Hs.eg.db, species = "hsa") {
  # Ensure the output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Extract the gene list in ENSEMBL ID format
  gene_list <- deg_df$ENSEMBL
  
  # Convert ENSEMBL IDs to ENTREZ IDs
  gene_entrez <- bitr(gene_list, fromType = "ENSEMBL", 
                      toType = "ENTREZID", 
                      OrgDb = org_db)
  
  # Ensure conversion worked
  if (nrow(gene_entrez) == 0) {
    stop("No valid ENTREZ IDs were found for the provided gene list.")
  }
  
  # Perform KEGG pathway analysis
  message("Performing KEGG pathway enrichment analysis...")
  
  kegg_result <- enrichKEGG(gene          = gene_entrez$ENTREZID,
                            organism      = species,
                            keyType       = "kegg",
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.05,
                            qvalueCutoff  = 0.2)
  
  # Check if there are any significant pathways
  if (is.null(kegg_result) || nrow(kegg_result@result) == 0) {
    warning("No significant KEGG pathways were found.")
  } else {
    # Save the KEGG results to a CSV file
    output_file <- file.path(output_dir, "KEGG_results.csv")
    write.csv(as.data.frame(kegg_result), output_file, row.names = FALSE)
    message("KEGG pathway analysis completed. Results saved in: ", output_file)
  }
  
  return(kegg_result)
}
save(KEGG_analysis, file = "KEGG_analysis.RData")
