# function for heatmap

heatmap_count <- function(data, type, top_n = 20) {
  # load
  library(DESeq2)
  library(pheatmap)
  # Ensure the type is valid
  if (!type %in% c("ntd", "vsd", "rld")) {
    stop("Invalid type specified. Use 'ntd', 'vsd', or 'rld'.")
  }
  
  # Select top N genes based on normalized counts
  select <- order(rowMeans(counts(data, normalized = TRUE)), 
                  decreasing = TRUE)[1:top_n]
  
  # Create annotation dataframe from colData
  df <- as.data.frame(colData(data)[, c("condition", "read_type")])
  colnames(df) <- c("Condition", "Read_Type")
  
  # Apply the appropriate transformation based on the type
  transformed_data <- switch(type,
                             ntd = assay(data),
                             vsd = assay(data),
                             rld = assay(data))
  
  # Rename columns to include condition in sample names
  conditions <- colData(data)$condition
  sample_names <- paste0("Sample_", seq_along(conditions), "_", conditions)
  colnames(transformed_data) <- sample_names
  
  # Generate the heatmap
  pheatmap(transformed_data[select, ], 
           cluster_rows = FALSE, 
           show_rownames = FALSE, 
           cluster_cols = FALSE, 
           annotation_col = df, 
           main = paste(type, "Transformation Heatmap"))
}

save(heatmap_count, file = "heatmap_count.RData")