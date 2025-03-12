#### RSCRIPT COX ANALYSIS - ALL METRICS 
#### HIGOR ALMEIDA, EMANNUEL, ENRIQUE MEDINA-ACOSTA 
#### LAST VERSION - 13/06/2024

# Directory
setwd("C:/Users/Emamnuell/Desktop/panoptosis/genes/genes_24_optosis/codigos/3°_cox")

# Load required packages
library(dplyr)
library(rio)
library(UCSCXenaTools)
library(UCSCXenaShiny)

# Dataframe
genes <- import("big_table_contacetado_FxG_miRNA.xlsx")

# Function to perform Cox analysis for a specific metric
perform_cox_analysis <- function(genes, measure) {
  results <- data.frame()  # Create an empty data frame to store results
  
  # Use unique to avoid unnecessary repetitions
  unique_genes <- unique(genes$genes)
  
  for (gene in unique_genes) {
    cox_result <- vis_unicox_tree(
      Gene = gene,
      measure = measure,
      data_type = "miRNA",
      threshold = 0.5
    )
    
    # Extract values of p.value, Type, and cancer
    p_value <- cox_result[["data"]][["p.value"]]
    Type <- cox_result[["data"]][["Type"]]
    cancer <- cox_result[["data"]][["cancer"]]
    
    # Create a temporary data frame with the current results
    temp_df <- data.frame(
      Gene = gene,
      p.value = p_value,
      Type = Type,
      cancer = cancer
    )
    
    # Append the temporary data frame to the main results
    results <- rbind(results, temp_df)
  }
  
  # Rename columns in results to avoid conflicts when merging
  colnames(results) <- c("Gene", paste0("p.value_Cox_", measure), paste0("Type_Cox_", measure), "cancer")
  
  # Merge the genes and results data frames by Gene and cancer_types
  merged_genes <- merge(genes, results[, c("Gene", paste0("p.value_Cox_", measure), paste0("Type_Cox_", measure), "cancer")], 
                        by.x = c("genes", "cancer_types"),
                        by.y = c("Gene", "cancer"),
                        all.x = TRUE)
  
  return(merged_genes)
}

# Apply the function for each metric and store results in the genes data frame
measures <- c("OS", "DSS", "DFI", "PFI")
genes_merged <- genes

# List to store completion times
end_times <- list()

for (measure in measures) {
  cat("Analyzing metric:", measure, "\n")  # Indicate the current metric
  genes_merged <- perform_cox_analysis(genes_merged, measure)
  end_times[[measure]] <- Sys.time()  # Store completion time
}

# Save the final result to a CSV file
export(genes_merged, "final_cox_results_miRNA.xlsx", row.names = FALSE)

# Print completion times
print(end_times)
