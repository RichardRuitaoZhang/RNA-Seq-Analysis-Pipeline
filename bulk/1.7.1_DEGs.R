# set working dir
setwd("/Users/richardzhang/Desktop/Northwestern/Y1/BME_311/Final Project/data/PC9/")

load("/Users/richardzhang/Desktop/Northwestern/Y1/BME_311/Final Project/data/get_exp.RData")
# result
exp_matrix <- list()

# Iterate through the list
for (dataset_name in names(results_list)) {
  # Get the DESeq2 result for the current dataset
  res <- results_list[[dataset_name]]
  
  # Generate a custom prefix for the output
  output_prefix <- paste0(dataset_name, "_processed")
  
  # Process the current result using your function and save it in the list
  exp_matrix[[dataset_name]] <- get_exp(res, output_prefix = output_prefix)
}

load("/Users/richardzhang/Desktop/Northwestern/Y1/BME_311/Final Project/data/get_DEGs.RData")
# Apply the save_DEGs function to each dataset in the list
for (sample in names(exp_matrix)) {
  res <- exp_matrix[[sample]]
  
  # Generate custom output prefix based on the dataset name
  output_prefix <- paste0(sample, "_deg_05")
  
  # Run the save_DEGs function
  get_DEGs(res, output_dir = output_prefix)
}
