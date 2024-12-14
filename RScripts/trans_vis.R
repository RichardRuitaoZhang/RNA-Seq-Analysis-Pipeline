library(DESeq2)
library(vsn)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)

trans_vis <- function(result, dds, coldata, output_dir="output") {
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  vsd <- vst(dds, blind = FALSE)
  rld <- rlog(dds, blind = FALSE)
  ntd <- normTransform(dds)
  
  # Save meanSd plots using png() with higher resolution
  meanSdPlot(assay(ntd))
  png(filename = file.path(output_dir, paste0("meanSd_ntd_", result, ".png")), 
      width = 1600, height = 1200, res = 300)
  print(meanSdPlot(assay(ntd))$gg)
  dev.off()
  
  meanSdPlot(assay(vsd))
  png(filename = file.path(output_dir, paste0("meanSd_vsd_", result, ".png")), 
      width = 1600, height = 1200, res = 300)
  print(meanSdPlot(assay(vsd))$gg)
  dev.off()
  
  meanSdPlot(assay(rld))
  png(filename = file.path(output_dir, paste0("meanSd_rld_", result, ".png")), 
      width = 1600, height = 1200, res = 300)
  print(meanSdPlot(assay(rld))$gg)
  dev.off()
  
  # Heatmaps
  select <- order(rowMeans(counts(dds, normalized = TRUE)), decreasing = TRUE)[1:20]
  df <- as.data.frame(coldata[, c("condition", "read_type")])
  
  png(filename = file.path(output_dir, paste0("heatmap_ntd_", result, ".png")), 
      width = 1600, height = 1200, res = 300)
  pheatmap(assay(ntd)[select, ], cluster_rows = FALSE, show_rownames = FALSE, 
           cluster_cols = FALSE, annotation_col = df, 
           main = paste("Heatmap (Normal Transformation) -", result))
  dev.off()
  
  png(filename = file.path(output_dir, paste0("heatmap_vsd_", result, ".png")), 
      width = 1600, height = 1200, res = 300)
  pheatmap(assay(vsd)[select, ], cluster_rows = FALSE, show_rownames = FALSE, 
           cluster_cols = FALSE, annotation_col = df, 
           main = paste("Heatmap (VST) -", result))
  dev.off()
  
  png(filename = file.path(output_dir, paste0("heatmap_rld_", result, ".png")), 
      width = 1600, height = 1200, res = 300)
  pheatmap(assay(rld)[select, ], cluster_rows = FALSE, show_rownames = FALSE, 
           cluster_cols = FALSE, annotation_col = df, 
           main = paste("Heatmap (RLD) -", result))
  dev.off()
  
  # Sample-to-sample distances heatmap
  sampleDists <- dist(t(assay(vsd)))
  sampleDistMatrix <- as.matrix(sampleDists)
  
  # Ensure the number of rows in coldata matches the number of samples in the matrix
  coldata <- coldata[rownames(sampleDistMatrix), , drop = FALSE]
  rownames(sampleDistMatrix) <- paste(coldata$condition, coldata$read_type, sep = "-")
  
  colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
  
  png(filename = file.path(output_dir, paste0("sample_distance_heatmap_", result, ".png")), 
      width = 1600, height = 1200, res = 300)
  pheatmap(sampleDistMatrix, 
           clustering_distance_rows = sampleDists, 
           clustering_distance_cols = sampleDists, 
           col = colors, 
           main = paste("Sample Distance Heatmap -", result))
  dev.off()
  
  # PCA plot
  pcaData <- plotPCA(vsd, intgroup = c("condition", "read_type"), returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  png(filename = file.path(output_dir, paste0("pca_plot_", result, ".png")), 
      width = 1600, height = 1200, res = 300)
  p <- ggplot(pcaData, aes(PC1, PC2, color = condition, shape = read_type)) +
    geom_point(size = 3) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    coord_fixed() +
    ggtitle(paste("PCA Plot -", result))
  print(p)
  dev.off()
  
  return(list(vsd = vsd, rld = rld, ntd = ntd))
}


heatmap_count <- function(data, coldata, type, top_n = 20) {
  library(DESeq2)
  library(pheatmap)
  
  # Ensure the type is valid
  if (!type %in% c("ntd", "vsd", "rld")) {
    stop("Invalid type specified. Use 'ntd', 'vsd', or 'rld'.")
  }
  
  # Select top N genes based on normalized counts
  select <- order(rowMeans(counts(data, normalized = TRUE)), 
                  decreasing = TRUE)[1:top_n]
  
  # Use provided coldata for annotations
  df <- as.data.frame(coldata[, c("condition", "read_type")])
  colnames(df) <- c("Condition", "Read_Type")
  
  # Apply the appropriate transformation based on the type
  transformed_data <- switch(type,
                             ntd = assay(data),
                             vsd = assay(data),
                             rld = assay(data))
  
  # Rename columns to include condition in sample names
  conditions <- coldata$condition
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

# Example of calling the function in a loop
# for (i in seq_along(results_list)) {
#   trans_vis(names(results_list)[i], dds_list[[i]], coldata_list[[i]], output_dir="results")
# }

save(trans_vis, file = "trans_vis.RData")

