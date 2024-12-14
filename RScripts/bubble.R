bubble <- function(data, category = NULL, title = "Bubble Plot of Terms", top_n = 20, output_file = NULL) {
  library(ggplot2)
  library(tidyr)
  library(RColorBrewer)
  library(dplyr)
  
  # Preprocess the data
  data <- data %>%
    mutate(
      GeneRatio = as.numeric(sub("/.*", "", GeneRatio)) / as.numeric(sub(".*/", "", GeneRatio)),  # Convert GeneRatio to numeric
      logP = -log10(pvalue),  # Compute -log10(p-value)
      Category = category     # Add category if provided
    ) %>%
    filter(pvalue < 0.05) %>%  # Filter for significant terms
    arrange(pvalue) %>%        # Sort by p-value
    head(top_n)                # Select the top_n terms
  
  # Create the bubble plot
  plot <- ggplot(data, aes(x = GeneRatio, y = reorder(Description, GeneRatio), size = Count, color = logP)) +
    geom_point(alpha = 0.7) +  # Add transparency
    scale_size(range = c(3, 10)) +  # Adjust bubble sizes
    scale_color_gradientn(colors = brewer.pal(9, "RdYlBu")) +  # Color scale
    theme_bw() +
    labs(
      title = title,
      x = "Gene Ratio",
      y = NULL,
      color = "-log10(p-value)",
      size = "Gene Count"
    ) +
    theme(axis.text.y = element_text(size = 10), plot.title = element_text(hjust = 0.5))  # Adjust axis text size and title alignment
  
  # Save the plot if output_file is specified
  if (!is.null(output_file)) {
    ggsave(output_file, plot, width = 10, height = 8)
    message("Bubble plot saved to: ", output_file)
  }
  
  return(plot)
}

save(bubble, file = "bubble.RData")
