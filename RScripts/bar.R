bar <- function(data, title = "Bar Plot of Terms", top_n = 20, category = NULL, output_dir = "bar_plots", output_file = NULL) {
  library(ggplot2)
  library(tidyr)
  library(RColorBrewer)
  library(dplyr)
  
  # Ensure the output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
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
  
  # Create the bar plot
  plot <- ggplot(data, aes(x = GeneRatio, y = reorder(Description, GeneRatio), fill = logP)) +
    geom_bar(stat = "identity", color = "black", alpha = 0.8) +
    scale_fill_gradientn(colors = brewer.pal(9, "RdYlBu"), name = "-log10(p-value)") +
    theme_bw() +
    labs(
      title = title,
      x = "Gene Ratio",
      y = NULL
    ) +
    theme(
      axis.text.y = element_text(size = 10),
      plot.title = element_text(hjust = 0.5)
    )
  
  # Set up the output file path
  if (!is.null(output_file)) {
    # Ensure the output file has an extension
    if (!grepl("\\.[a-zA-Z]+$", output_file)) {
      output_file <- paste0(output_file, ".png")
    }
    # Sanitize the filename to replace spaces and special characters
    sanitized_file <- gsub("[^a-zA-Z0-9._-]", "_", output_file)
    output_path <- file.path(output_dir, sanitized_file)
    ggsave(output_path, plot, width = 10, height = 8)
    message("Bar plot saved to: ", output_path)
  }
  
  return(plot)
}

save(bar, file = "bar.RData")
