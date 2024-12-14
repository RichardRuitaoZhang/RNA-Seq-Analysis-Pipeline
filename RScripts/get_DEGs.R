# function for get the differential expression genes
get_DEGs <- function(res_df, output_dir = "output_folder") {
  # Ensure the output directory exists, create it if not
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Get DEGs based on adjusted p-value (adj.p < 0.05) and log2 Fold Change (|logFC| > 1)
  degenes <- res_df %>% filter(padj < 0.05 & abs(log2FoldChange) > 1)
  
  # Mark if a gene is up- or down-regulated
  degenes$regulation <- ifelse(degenes$log2FoldChange > 1, 'up', 'down')
  
  # Print summary of up- and down-regulated genes
  cat(paste("Up-regulated genes:", nrow(subset(degenes, regulation == "up")), "\n"))
  cat(paste("Down-regulated genes:", nrow(subset(degenes, regulation == "down")), "\n"))
  
  # Define file paths
  all_file <- file.path(output_dir, "deg_05.csv")
  up_file <- file.path(output_dir, "deg_05_up.csv")
  down_file <- file.path(output_dir, "deg_05_down.csv")
  
  # Save DEG list to the directory
  write.csv(degenes, all_file, row.names = FALSE)
  write.csv(subset(degenes, regulation == "up"), up_file, row.names = FALSE)
  write.csv(subset(degenes, regulation == "down"), down_file, row.names = FALSE)
  
  # Return the DEG data frame for further use
  return(degenes)
}
save(get_DEGs, file = "get_DEGs.RData")