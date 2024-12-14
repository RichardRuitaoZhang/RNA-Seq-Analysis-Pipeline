# bubble plot
bubble <- function(data, category, title = "Bubble Plot of GO Terms", output_file = NULL) {
  library(ggplot2)
  library(tidyr)
  library(RColorBrewer)
  
  # Preprocess the data
  data$logP <- -log10(data$PValue)  # Compute -log10(p-value)
  data <- data %>%
    separate(Term, into = c("GO_ID", "Term"), sep = "~")  # Split GO ID and Term
  data$Category <- rep(category, nrow(data))  # Add category label
  data <- data[data$PValue < 0.05, ]  # Filter for significant terms
  
  # Create the bubble plot
  plot <- ggplot(data, aes(x = Fold.Enrichment, y = reorder(Term, Fold.Enrichment), size = Count, color = logP)) +
    geom_point(alpha = 0.7) +  # Add transparency
    scale_size(range = c(3, 10)) +  # Adjust bubble sizes
    scale_color_gradientn(colors = brewer.pal(9, "RdYlBu")) +  # Color scale
    theme_bw() +
    labs(
      title = title,
      x = "Fold Enrichment",
      y = "GO Term",
      color = "-log10(p-value)",
      size = "Gene Count"
    ) +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 10))  # Adjust axis text size for clarity
  
  # Save the plot if output_file is specified
  if (!is.null(output_file)) {
    ggsave(output_file, plot, width = 10, height = 8)
    message("Bubble plot saved to: ", output_file)
  }
  
  return(plot)
}
