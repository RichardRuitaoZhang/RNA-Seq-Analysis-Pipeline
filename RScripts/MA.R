# function for generate lfcshrink and plot MA
# because we are interested in treated vs untreated, we set 'coef=2'

MA <- function(dds_list, coef = 2, output_dir = "MAplots", shr_type = "ALL", plot = TRUE, names_list = NULL) {
  # Validate shr_type
  valid_types <- c("normal", "apeglm", "ashr", "ALL")
  if (!shr_type %in% valid_types) {
    stop("Invalid `shr_type`. Choose from 'normal', 'apeglm', 'ashr', or 'ALL'.")
  }
  
  # Validate plot
  if (!is.logical(plot)) {
    stop("`plot` must be TRUE or FALSE.")
  }
  
  # Initialize result lists
  resNorm_list <- list()
  resApm_list <- list()
  resAsh_list <- list()
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Iterate through the dds_list
  for (i in seq_along(dds_list)) {
    dataset_name <- if (!is.null(names(dds_list))) names(dds_list)[i] else paste0("dataset_", i)
    plot_title <- if (!is.null(names_list) && i <= length(names_list)) names_list[i] else dataset_name
    
    # Apply LFC shrinkage based on shr_type
    resNorm <- resApm <- resAsh <- NULL
    if (shr_type == "normal" || shr_type == "ALL") {
      resNorm <- tryCatch(lfcShrink(dds_list[[i]], coef = coef, type = "normal"),
                          error = function(e) { warning(paste("Error in normal for", dataset_name)); NULL })
      if (!is.null(resNorm)) resNorm_list[[dataset_name]] <- resNorm
    }
    if (shr_type == "apeglm" || shr_type == "ALL") {
      resApm <- tryCatch(lfcShrink(dds_list[[i]], coef = coef, type = "apeglm"),
                         error = function(e) { warning(paste("Error in apeglm for", dataset_name)); NULL })
      if (!is.null(resApm)) resApm_list[[dataset_name]] <- resApm
    }
    if (shr_type == "ashr" || shr_type == "ALL") {
      resAsh <- tryCatch(lfcShrink(dds_list[[i]], coef = coef, type = "ashr"),
                         error = function(e) { warning(paste("Error in ashr for", dataset_name)); NULL })
      if (!is.null(resAsh)) resAsh_list[[dataset_name]] <- resAsh
    }
    
    # Skip plotting if plot = FALSE
    if (!plot) next
    
    # Generate MA plots
    pdf_file <- file.path(output_dir, paste0(dataset_name, "_MA_plots.pdf"))
    pdf(pdf_file, width = 10, height = 4)
    
    # Ensure the device is closed after this iteration
    on.exit(dev.off(), add = TRUE)
    
    # Set plotting layout based on shr_type
    if (shr_type == "ALL") {
      par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))
    } else {
      par(mfrow = c(1, 1), mar = c(4, 4, 2, 1))
    }
    
    # Dynamic axis limits
    xlim <- range(c(resNorm$baseMean, resApm$baseMean, resAsh$baseMean), na.rm = TRUE)
    ylim <- range(c(resNorm$log2FoldChange, resApm$log2FoldChange, resAsh$log2FoldChange), na.rm = TRUE)
    
    # Plot based on shr_type
    if (shr_type == "normal") plotMA(resNorm, xlim = xlim, ylim = ylim, main = paste("Normal -", plot_title))
    if (shr_type == "apeglm") plotMA(resApm, xlim = xlim, ylim = ylim, main = paste("Apeglm -", plot_title))
    if (shr_type == "ashr") plotMA(resAsh, xlim = xlim, ylim = ylim, main = paste("Ashr -", plot_title))
    if (shr_type == "ALL") {
      plotMA(resApm, xlim = xlim, ylim = ylim, main = paste("Apeglm -", plot_title))
      plotMA(resNorm, xlim = xlim, ylim = ylim, main = paste("Normal -", plot_title))
      plotMA(resAsh, xlim = xlim, ylim = ylim, main = paste("Ashr -", plot_title))
    }
    
    #dev.off() # Close the PDF device
  }
  
  # Return the lists as a named list
  return(list(resNorm_list = resNorm_list, resApm_list = resApm_list, resAsh_list = resAsh_list))
}

save(MA, file = "MA.RData")
