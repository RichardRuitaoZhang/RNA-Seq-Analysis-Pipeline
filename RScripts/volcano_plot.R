# function for volcano plot
volcano_plot <- function(volcano, title = NULL, label = TRUE, label.size = 2, 
                         oncogenes = NULL, padj_threshold = 0.05, lfc_threshold = 1, 
                         out_dir = "volcano_results") {
  library(ggplot2)
  library(ggrepel)
  library(colorspace)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(AnnotationDbi)
  
  # Ensure output directory exists
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  # Get the dataset
  volcano <- as.data.frame(volcano)
  volcano$ENSEMBL <- sub("\\.\\d+$", "", rownames(volcano)) # Remove everything after the dot
  
  # Map ENSEMBL IDs to Gene Symbols
  gene_annotations <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = volcano$ENSEMBL,
    columns = c("SYMBOL"),      # Retrieve gene symbols
    keytype = "ENSEMBL"         # Input IDs are ENSEMBL
  )
  
  volcano <- merge(volcano, gene_annotations, by = "ENSEMBL")
  
  # Categorize genes as upregulated, downregulated, or not significant
  volcano$significance <- "NS"
  volcano$significance[volcano$padj < padj_threshold & volcano$log2FoldChange > lfc_threshold] <- "Upregulated"
  volcano$significance[volcano$padj < padj_threshold & volcano$log2FoldChange < -lfc_threshold] <- "Downregulated"
  
  # Filter the dataset for oncogenes if provided
  if (!is.null(oncogenes)) {
    volcano$oncogene_label <- ifelse(volcano$SYMBOL %in% oncogenes, volcano$SYMBOL, NA)
  } else {
    volcano$oncogene_label <- NA
  }
  
  # Save the filtered oncogenes data
  if (!is.null(oncogenes)) {
    oncogene_data <- volcano[volcano$SYMBOL %in% oncogenes, ]
    write.csv(oncogene_data, file = file.path(out_dir, paste(title, "oncogenes.csv")), row.names = FALSE)
  }
  
  # Create volcano plot
  p <- ggplot(volcano, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = significance), alpha = 0.5, size = 0.7) +
    scale_color_manual(values = c("Upregulated" = "#009E8E", "Downregulated" = "#AC7F21", "NS" = "grey")) +
    theme_bw() +
    geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", color = "black", size = 0.5) +
    geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), linetype = "dashed", color = "black", size = 0.5) +
    labs(
      title = title,
      x = "Log2 Fold Change",
      y = "-Log10 Adjusted P-Value",
      color = "Significance"
    ) +
    theme(plot.title = element_text(hjust = 0.5, size = 14)) +
    geom_text_repel(
      data = volcano[!is.na(volcano$oncogene_label), ],
      aes(x = log2FoldChange, y = -log10(padj), label = oncogene_label),
      size = label.size, box.padding = unit(0.5, "lines"),
      point.padding = unit(0.8, "lines"), 
      segment.color = "black", 
      show.legend = FALSE,
      max.overlaps = Inf,
      segment.alpha = 0.5
    )
  
  # Save the plot
  ggsave(filename = file.path(out_dir, paste0(title, "_volcano_plot.png")), plot = p, width = 8, height = 6, dpi = 300)
  
  return(p)
}

save(volcano_plot, file = "volcano_plot.RData")
