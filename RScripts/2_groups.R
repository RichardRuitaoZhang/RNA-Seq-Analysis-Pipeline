


# function (2-group comparison, latest version)
DESeq_multi_group <- function(exprSet, meta, coldata, repNum1 = NULL, repNum2 = NULL, 
                              test = NULL, control = NULL, 
                              min.size = 3, Sum = 10, ref = NULL, lfc = "normal") {
  # Load required libraries
  library(DESeq2)
  library(dplyr)
  
  # Helper function: Generate condition and subset count data
  get_condition <- function(exprSet, meta, coldata, group1, group2) {
    # Extract sample names for the groups
    group1_samples <- unlist(meta[group1], use.names = FALSE)
    group2_samples <- unlist(meta[group2], use.names = FALSE)
    
    # Subset expression set for the two groups
    samples <- c(group1_samples, group2_samples)
    count_data <- exprSet[, samples]
    
    # Prepare colData
    colData <- coldata[samples,]
    
    return(list(count_data = count_data, colData = colData))
  }
  
  # If test and control are not specified, perform all pairwise comparisons
  # Identify unique group names in metadata
  groups <- colnames(meta)
  
  # Initialize list to store results
  dds_list <- list()
  results_list <- list()
  resLFC_list <- list()
  
  # Perform pairwise comparisons
  for (i in 1:(length(groups) - 1)) {
    for (j in (i + 1):length(groups)) {
      group1 <- groups[i]
      group2 <- groups[j]
      
      # Subset data and create condition
      data <- get_condition(exprSet, meta, coldata, group1, group2)
      
      # Run DESeq2 analysis
      dds <- DESeqDataSetFromMatrix(countData = data$count_data, colData = data$colData, design = ~ condition)
      
      # Pre-filter low-count genes
      keep <- rowSums(counts(dds) >= Sum) >= min.size
      dds <- dds[keep, ]
      
      # Relevel the reference level if specified
      if (is.null(ref)) {
        dds$condition <- relevel(dds$condition, ref = group1)  # Default: relevel to group1
      }
      
      # Run DESeq2 analysis
      dds <- DESeq(dds)
      res <- results(dds, name=paste("condition", group2, "vs", group1, sep = "_"))
      res <- results(dds, contrast = c("condition", group2, group1))
      resLFC <- lfcShrink(dds, coef=paste("condition", group2, "vs", group1, sep = "_"), type = lfc)
      register(MulticoreParam(4))
      resOrdered <- res[order(res$pvalue),] # sort res by the smallest p-value
      
      # Store results in the list
      dds_list[[paste(group2, "vs", group1, sep = "_")]] <- dds
      results_list[[paste(group2, "vs", group1, sep = "_")]] <- res
      resLFC_list[[paste(group2, "vs", group1, sep = "_")]] <- resLFC
    }
  }
  
  # Return the list of results
  return(list(results_list = results_list, dds_list = dds_list, resLFC_list = resLFC_list))
}


# If test and control are specified, perform two-group comparison
if (!is.null(test) & !is.null(control)) {
  # Validate replicates
  if (!is.null(repNum1) & !is.null(repNum2)) {
    if (ncol(exprSet[, unlist(meta[test])]) != repNum1 |
        ncol(exprSet[, unlist(meta[control])]) != repNum2) {
      stop("Mismatch between provided replicate numbers and sample counts.")
    }
  }
  
  
  # Subset data and create condition
  data <- get_condition(exprSet, meta, test, control)
  
  # Run DESeq2 analysis
  dds <- DESeqDataSetFromMatrix(countData = data$count_data, colData = data$colData, design = ~ condition)
  
  # Pre-filter low-count genes
  keep <- rowSums(counts(dds) >= Sum) >= min.size
  dds <- dds[keep, ]
  
  # Relevel the reference level if specified
  if (!is.null(ref)) {
    dds$condition <- relevel(dds$condition, ref = ref)
  }
  
  # Run DESeq2 analysis
  dds <- DESeq(dds)
  res <- results(dds, name=paste("condition", test, "vs", control, sep = "_"))
  res <- results(dds, contrast = c("condition", test, control))
  
  # Return results as a data frame
  return(list(results_list = list(single_result = as.data.frame(res)), 
              dds_list = list(single_dds = dds)))
}
