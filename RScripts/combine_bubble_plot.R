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
  
  # Prepare data for each ontology
  bp_data <- as.data.frame(go_results_list[["BP"]])
  cc_data <- as.data.frame(go_results_list[["CC"]])
  mf_data <- as.data.frame(go_results_list[["MF"]])
  
  bp_data <- bp_data %>%
    mutate(GeneRatio = as.numeric(sub("/.*", "", GeneRatio)) / as.numeric(sub(".*/", "", GeneRatio)))
  cc_data <- cc_data %>%
    mutate(GeneRatio = as.numeric(sub("/.*", "", GeneRatio)) / as.numeric(sub(".*/", "", GeneRatio)))
  mf_data <- mf_data %>%
    mutate(GeneRatio = as.numeric(sub("/.*", "", GeneRatio)) / as.numeric(sub(".*/", "", GeneRatio)))
  
  # Create individual bubble plots
  plot_bp <- bubble(bp_data, category = "Biological Process", title = "Biological Process", top_n)
  plot_cc <- bubble(cc_data, category = "Cellular Component", title = "Cellular Component", top_n)
  plot_mf <- bubble(mf_data, category = "Molecular Function", title = "Molecular Function", top_n)
  
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
save(combined_bubble_plot, file = "combined_bubble_plot.RData")