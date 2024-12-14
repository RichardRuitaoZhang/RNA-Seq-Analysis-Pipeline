combined_bubble_plot2 <- function(go_results_list, top_n = 10, title = "Combined GO Bubble Plot", output_dir = "bubble_plots") {
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(patchwork)  # For combining multiple plots
  
  # Ensure the output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Helper function to prepare data for plotting
  prepare_data <- function(data, category) {
    prepared_data <- data %>%
      filter(p.adjust < 0.05) %>%
      top_n(-top_n, p.adjust) %>%
      mutate(
        GeneRatio = as.numeric(sub("/.*", "", GeneRatio)) / as.numeric(sub(".*/", "", GeneRatio)),
        Category = category
      ) %>%
      arrange(p.adjust)
    
    return(prepared_data)
  }
  
  # Combine data for all categories, skipping if no significant terms
  combined_data <- list()
  
  if (nrow(as.data.frame(go_results_list[["BP"]]) %>% filter(p.adjust < 0.05)) > 0) {
    bp_data <- prepare_data(as.data.frame(go_results_list[["BP"]]), "Biological Process")
    combined_data <- append(combined_data, list(bp_data))
  } else {
    message("Biological Process doesn't contain any significant terms.")
  }
  
  if (nrow(as.data.frame(go_results_list[["CC"]]) %>% filter(p.adjust < 0.05)) > 0) {
    cc_data <- prepare_data(as.data.frame(go_results_list[["CC"]]), "Cellular Component")
    combined_data <- append(combined_data, list(cc_data))
  } else {
    message("Cellular Component doesn't contain any significant terms.")
  }
  
  if (nrow(as.data.frame(go_results_list[["MF"]]) %>% filter(p.adjust < 0.05)) > 0) {
    mf_data <- prepare_data(as.data.frame(go_results_list[["MF"]]), "Molecular Function")
    combined_data <- append(combined_data, list(mf_data))
  } else {
    message("Molecular Function doesn't contain any significant terms.")
  }
  
  # Check if there is any data to plot
  if (length(combined_data) == 0) {
    message("No significant terms found in any category. Skipping plot generation.")
    return(NULL)
  }
  
  combined_data <- bind_rows(combined_data)
  
  # Create the combined bubble plot
  combined_plot <- ggplot(combined_data, aes(x = GeneRatio, y = Description, size = Count, color = p.adjust)) +
    geom_point(alpha = 0.8) +
    facet_wrap(~Category, scales = "free_y", ncol = 1, strip.position = "right") +  # Arrange vertically with labels on the right
    scale_color_gradient(low = "blue", high = "red", name = "p.adjust") +
    scale_size(range = c(3, 10), name = "Count") +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_blank(),
      strip.text.y = element_text(size = 12, face = "bold"),  # Bold category labels on the right
      strip.background = element_rect(fill = "grey80", color = "black"),  # Grey background for category labels
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      panel.background = element_rect(fill = "white", color = NA),  # White plot background
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    ) +
    labs(x = "Gene Ratio", y = NULL, title = title)
  
  # Save the combined plot
  combined_output_file <- file.path(output_dir, paste(title, "combined_bubble_plot.png"))
  ggsave(combined_output_file, combined_plot, width = 10, height = 12)
  
  message("Combined plot saved to: ", combined_output_file)
  
  return(combined_plot)
}
save(combined_bubble_plot2, file = "combined_bubble_plot2.RData")