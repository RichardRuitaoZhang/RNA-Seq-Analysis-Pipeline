# function of get the result list into readable dataframe with gene id

get_exp <- function(res, output_prefix = "processed_deseq") {
  library(org.Hs.eg.db) # Ensure this library is installed
  library(dplyr)
  
  # Convert DESeq2 results into a dataframe
  res_df <- as.data.frame(res)
  
  # Remove version numbers from Ensembl IDs (everything after the dot)
  res_df$ENSEMBL <- sub("\\.\\d+$", "", rownames(res_df))
  
  # Map ENSEMBL IDs to Gene Symbols
  gene_annotations <- select(
    org.Hs.eg.db,
    keys = res_df$ENSEMBL,
    columns = c("SYMBOL"),      # Retrieve gene symbols
    keytype = "ENSEMBL"         # Input IDs are ENSEMBL
  )
  
  # Merge DESeq2 results with gene annotations
  res_df <- merge(res_df, gene_annotations, by = "ENSEMBL", all.x = TRUE)
  
  # Save the merged results to a CSV file
  output_file <- paste0(output_prefix, ".csv")
  write.csv(res_df, output_file, row.names = FALSE)
  
  # Return the merged data frame for further use
  return(res_df)
}
save(get_exp, file = "get_exp.RData")