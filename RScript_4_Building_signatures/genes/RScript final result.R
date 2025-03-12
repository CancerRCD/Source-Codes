# Load necessary packages
library(dplyr)
library(purrr)
library(fs)
library(readr)
library(tidyr)
library(rio)

#### PART A ####
# Define the base path
base_path <- "C:/Users/Emamnuell/Desktop/panoptosis/PANoptosis_radar_analysis"

# Identify the phenotypic variables and their respective paths
phenotypes <- c("MSI", "Stemness", "TMB")  # Phenotypic prefixes
phenotypic_suffix <- "_PanCan"  # Common suffix for phenotypes
genotypes <- c("CNV", "Methylation", "mRNA", "Mutation", "Protein")

# Function to find the .tsv file with the highest suffix in a folder
get_max_suffix_tsv <- function(path) {
  files <- dir_ls(path, regexp = "gene_cancer_correlation_network_.*\\.tsv$")
  if (length(files) == 0) return(NULL)
  max_file <- files[which.max(as.numeric(stringr::str_extract(basename(files), "\\d\\.\\d")))]
  return(max_file)
}

# Function to import and annotate the data
import_and_label_data <- function(phenotype, genotype) {
  phenotype_dir = paste(phenotype, phenotypic_suffix, sep="")
  folder_path <- file.path(base_path, phenotype_dir, paste(genotype, "vs", phenotype, sep = "_"))
  result_path <- file.path(folder_path, "results")
  
  if (!dir.exists(result_path)) {
    return(NULL)  # Returns NULL if the directory does not exist
  }
  
  file_path <- get_max_suffix_tsv(result_path)
  
  if (is.null(file_path)) {
    return(NULL)  # Returns NULL if the .tsv file is not found
  }
  
  data <- read_tsv(file_path, col_types = cols())
  mutate(data, var_fenotipica = phenotype, var_genotipica = genotype)
}

# Process all files, combining them into a single dataframe
all_data <- map(phenotypes, function(phenotype) {
  map_df(genotypes, ~ import_and_label_data(phenotype, .x))
}) %>% bind_rows() %>% filter(!is.null(.))

# View the general table
print(all_data)



#### PART B ####
##### MERGING OPTOSIS CATEGORIES WITH GENES FOUND IN SIGNATURES
setwd("C:/Users/Emamnuell/Desktop/panoptosis/genes")

# Import Data
Pathways <- import("ALL_genes_five_optosis_final_com_vias.xlsx") # Provided by oct4h
Cor <- all_data

# Merge 'Pathways' and 'Cor' based on the gene
combined_table <- merge(Cor, Pathways, by.x = "genes", by.y = "name", all.x = TRUE)

# Aggregate all optosis categories associated with each gene into a single cell
final_table <- aggregate(Pathway_Categorization ~ genes, combined_table, FUN = function(x) paste(unique(x), collapse = ", "))

# Aggregate all optosis categories associated with each gene into a single cell
final_table <- aggregate(cbind(Pathway_Categorization) ~ genes, combined_table, FUN = function(x) paste(unique(x), collapse = ", "))

# Merge aggregated data with the original table
final_table <- merge(Cor, final_table, by = "genes", all.x = TRUE, all.y = FALSE)

# Reorder columns
New_column_order <- c("cancer_types", "genes", "Pathway_Categorization", "correlation", "var_fenotipica", "var_genotipica")

final_table <- final_table[New_column_order]

# Sort observations in the first column alphabetically
Tabela_Cancer_Genes_Pathways <- final_table[order(final_table$cancer_types), ]



##### PART C ####
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

# Function to aggregate genes by cancer type, optosis, phenotypic variable, and genotypic variable, separating by correlation
aggregate_genes <- function(data) {
  # Extract the first optosis
  data <- data %>%
    mutate(Pathway_Categorization = extract_first_optosis(Pathway_Categorization))
  
  # Separate data with positive and negative correlation
  positive <- data %>% filter(correlation > 0)
  negative <- data %>% filter(correlation < 0)
  
  # Internal aggregation function with custom formatting
  aggregate_and_format <- function(df) {
    df %>%
      group_by(cancer_types, Pathway_Categorization, var_fenotipica, var_genotipica) %>%
      summarise(genes = format_gene_signature(genes), .groups = 'drop')
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
# Assume Tabela_Cancer_Genes_Pathways is already loaded in the R environment
aggregated_df <- aggregate_genes(Tabela_Cancer_Genes_Pathways)

# Add a column that counts the number of genes in each signature
aggregated_df <- aggregated_df %>%
  mutate(number_of_genes = sapply(strsplit(genes, " \\+ "), length))

# Rename the DataFrame columns
aggregated_df <- aggregated_df %>%
  rename(
    Cancer_type = cancer_types,
    Gene_signature = genes,
    Optosis_Category = Pathway_Categorization,
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



##### PART D ####
# Re-add the 'signature_position' column and format it (if necessary)
aggregated_df <- aggregated_df %>%
  group_by(Cancer_type) %>%
  mutate(signature_position = sprintf("%02d", row_number())) %>%
  ungroup()

# Correcting the creation of the 'Name' column
aggregated_df_names <- aggregated_df %>%
  mutate(
    # Creating the new 'Name' column according to the specified format
    Name = paste(
      Cancer_type, 
      paste(
        signature_position,
        case_when(
          Phenotype == "MSI" ~ "1",
          Phenotype == "Stemness" ~ "2",
          Phenotype == "TMB" ~ "3",
          TRUE ~ NA_character_  # Keeps NA if none of the cases match
        ),
        case_when(
          Genotype == "CNV" ~ "1",
          Genotype == "Methylation" ~ "2",
          Genotype == "mRNA" ~ "3",
          Genotype == "Mutation" ~ "4",
          Genotype == "Protein" ~ "5",
          TRUE ~ NA_character_  # Keeps NA if none of the cases match
        ),
        case_when(
          Correlation_sign == "positive" ~ "P",
          Correlation_sign == "negative" ~ "N",
          TRUE ~ NA_character_  # Keeps NA if none of the cases match
        ),
        sep = "."
      ),
      sep = "-"
    )
  )

# Checking the first rows of the modified dataframe to confirm the changes
print(head(aggregated_df_names))

# Export results to a TSV file
export(aggregated_df_names, "aggregated_df_names.tsv")
