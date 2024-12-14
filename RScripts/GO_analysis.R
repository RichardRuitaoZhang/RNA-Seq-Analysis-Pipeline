
GO_analysis <- function(deg_df, output_dir, org_db = org.Hs.eg.db) {
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
  
  # Initialize a list to store GO results for BP, CC, MF
  go_results_list <- list()
  
  # Perform GO analysis for BP, CC, MF
  for (ontology in c("BP", "CC", "MF")) {
    message("Performing GO analysis for ontology: ", ontology)
    
    go_result <- enrichGO(gene          = gene_entrez$ENTREZID,
                          OrgDb         = org_db,
                          keyType       = "ENTREZID",
                          ont           = ontology,
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.05,
                          qvalueCutoff  = 0.2)
    
    # Add the 'type' column dynamically when converted to a DataFrame
    go_result@result$type <- switch(ontology,
                                    "BP" = "Biological Process",
                                    "CC" = "Cellular Component",
                                    "MF" = "Molecular Function")
    
    # Save result to the list
    go_results_list[[ontology]] <- go_result
    
    # Save the results as a CSV file
    output_file <- file.path(output_dir, paste0("GO_", ontology, "_results.csv"))
    write.csv(as.data.frame(go_result), output_file, row.names = FALSE)
  }
  
  message("GO analysis completed. Results saved in: ", output_dir)
  
  return(go_results_list)
}
save(GO_analysis, file = "GO_analysis.RData")