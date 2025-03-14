# Load necessary packages
library(dplyr)
library(purrr)
library(fs)
library(readr)
library(tidyr)
library(rio)

# Set the working directory
setwd("C:/Users/Emamnuell/Desktop/panoptosis/genes/genes_24_optosis/codigos/4°_signature/miRNA")

# Import the final Cox results table
Tabela_Cancer_Genes_Pathways <- import("final_cox_results_miRNA_com_vias.xlsx")

# Function to extract the first optosis category
extract_first_optosis <- function(Pathway_Categorization) {
  sapply(strsplit(Pathway_Categorization, ","), function(x) trimws(x[1]))
}

# Function to format the gene signature according to requirements
format_gene_signature <- function(genes) {
  formatted_genes <- sapply(genes, function(gene) {
    if (grepl("-", gene)) {
      paste0("`", gene, "`")
    } else {
      gene
    }
  })
  signature <- paste(formatted_genes, collapse = " + ")
  paste0("(", signature, ")")
}

# Function to aggregate miRNAs by cancer type, optosis, phenotypic variable, and genotypic variable, separating by correlation
aggregate_genes <- function(data) {
  # Extract the first optosis category
  data <- data %>%
    mutate(Pathway_Categorization = extract_first_optosis(Pathway_Categorization))
  
  # Separate data with positive and negative correlation
  positive <- data %>% filter(correlation > 0)
  negative <- data %>% filter(correlation < 0)
  
  # Internal aggregation function with custom formatting
  aggregate_and_format <- function(df) {
    df %>%
      group_by(cancer_types, Pathway_Categorization, var_fenotipica, var_genotipica, Expression, Type_Cox_OS, Type_Cox_DSS, Type_Cox_DFI, Type_Cox_PFI) %>%
      summarise(miRNA_mature = format_gene_signature(miRNA_mature), .groups = 'drop')
  }
  
  # Aggregate and format for positive and negative correlations
  positive_aggregated <- aggregate_and_format(positive)
  negative_aggregated <- aggregate_and_format(negative)
  
  # Add a column to identify the correlation type
  positive_aggregated <- positive_aggregated %>%
    mutate(correlation_type = "positive")
  negative_aggregated <- negative_aggregated %>%
    mutate(correlation_type = "negative")
  
  # Combine both aggregated dataframes
  aggregated_data <- bind_rows(positive_aggregated, negative_aggregated)
  return(aggregated_data)
}

# Apply the function to the example dataframe
aggregated_df <- aggregate_genes(Tabela_Cancer_Genes_Pathways)

# Add a column that counts the number of miRNAs in each signature
aggregated_df <- aggregated_df %>%
  mutate(number_of_genes = sapply(strsplit(miRNA_mature, " \\+ "), length))

# Rename the DataFrame columns
aggregated_df <- aggregated_df %>%
  rename(
    Cancer_type = cancer_types,
    miRNA_signature = miRNA_mature,
    RCD = Pathway_Categorization,
    Correlation_sign = correlation_type,
    Phenotype = var_fenotipica,
    Genotype = var_genotipica,
    Members = number_of_genes
  )

# Sort by 'Cancer_type' and 'Members' in descending order within each cancer type
aggregated_df <- aggregated_df %>%
  arrange(Cancer_type, desc(Members))

# Create a 'signature_position' column that enumerates entries for each cancer type
aggregated_df <- aggregated_df %>%
  group_by(Cancer_type) %>%
  mutate(Signature_Position = row_number()) %>%
  ungroup()

# Reorder columns to place 'signature_position' as the first column
aggregated_df <- aggregated_df %>%
  select(Signature_Position, everything())

# View final results
print(aggregated_df)
View(aggregated_df)

# Export results to an Excel file
export(aggregated_df, "signature_miRNA.xlsx")
