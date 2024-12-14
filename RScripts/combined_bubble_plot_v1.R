# combine bubble plot
combined_bubble_plot <- function(go_results_list, top_n = 10, title = "Combined GO Bubble Plot", output_dir = "bubble_plots") {
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(patchwork)  # For combining multiple plots
  
  # Ensure the output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Helper function to process GO result for bubble plotting
  prepare_data <- function(go_result, category, top_n) {
    data <- as.data.frame(go_result)
    data$logP <- -log10(data$pvalue)  # Add -log10(p-value)
    data <- data %>%
      separate(Description, into = c("GO_ID", "Term"), sep = "~", fill = "right")  # Separate GO ID and Term
    data <- data[data$pvalue < 0.05, ]  # Filter significant terms
    data$Category <- category  # Add category
    data <- data %>% arrange(desc(logP)) %>% head(top_n)  # Select top_n terms
    return(data)
  }
  
  # Prepare data for each ontology
  bp_data <- prepare_data(go_results_list[["BP"]], "Biological Process", top_n)
  cc_data <- prepare_data(go_results_list[["CC"]], "Cellular Component", top_n)
  mf_data <- prepare_data(go_results_list[["MF"]], "Molecular Function", top_n)
  
  # Create individual bubble plots
  plot_bp <- bubble(bp_data, category = "Biological Process", title = "Biological Process")
  plot_cc <- bubble(cc_data, category = "Cellular Component", title = "Cellular Component")
  plot_mf <- bubble(mf_data, category = "Molecular Function", title = "Molecular Function")
  
  # Save individual plots
  ggsave(file.path(output_dir, "bubble_plot_BP.png"), plot_bp, width = 8, height = 6)
  ggsave(file.path(output_dir, "bubble_plot_CC.png"), plot_cc, width = 8, height = 6)
  ggsave(file.path(output_dir, "bubble_plot_MF.png"), plot_mf, width = 8, height = 6)
  
  # Combine plots vertically
  combined_plot <- plot_bp / plot_cc / plot_mf + plot_annotation(title = title)
  
  # Save the combined plot
  combined_output_file <- file.path(output_dir, "combined_bubble_plot.png")
  ggsave(combined_output_file, combined_plot, width = 10, height = 12)
  
  message("All plots saved to: ", output_dir)
  
  return(combined_plot)
}
