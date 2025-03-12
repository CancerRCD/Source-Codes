# Load necessary packages
library(rio)
library(dplyr)
library(stringr)

# Set the working directory
setwd("C:/Users/Emamnuell/Desktop/panoptosis/genes/genes_24_optosis/codigos/4°_signature/all")

# Import data frames
signatures <- import("signature_protein.xlsx")
main_concatened <- import("big_table_contacetado_signatures_protein.tsv")

# Function to process compound and simple gene names
process_gene_name <- function(gene) {
  if (grepl("\\+", gene)) {
    return(gene)
  } else {
    return(gsub("`|\\(|\\)", "", gene))
  }
}

# Apply the processing function to the signature column
signatures$signature <- sapply(signatures$Protein_signature, process_gene_name)

# Transform the word "transcript" into "Transcript" in the Genotype column
signatures$Genotype <- str_replace_all(signatures$Genotype, "\\btranscript\\b", "Transcript")

# Filter main_concatened based on values present in signatures
filtered_main_concatened <- main_concatened %>%
  semi_join(signatures, by = c("cancer_types" = "Cancer_type", 
                               "var_fenotipica" = "Phenotype", 
                               "var_genotipica" = "Genotype", 
                               "genes" = "signature"))

# Identify rows in signatures that are not in main_concatened
not_in_main <- anti_join(signatures, main_concatened, 
                         by = c("Cancer_type" = "cancer_types", 
                                "Phenotype" = "var_fenotipica", 
                                "Genotype" = "var_genotipica", 
                                "signature" = "genes"))

# Function to dynamically create data frames
create_data_frames <- function(data) {
  unique_combinations <- unique(data[c("Phenotype", "Genotype")])
  
  for (i in 1:nrow(unique_combinations)) {
    phenotype <- unique_combinations$Phenotype[i]
    genotype <- unique_combinations$Genotype[i]
    
    subset_data <- data %>%
      filter(Phenotype == phenotype & Genotype == genotype)
    
    df_name <- paste0("Phenotype_", phenotype, "_Genotype_", genotype)
    assign(df_name, subset_data, envir = .GlobalEnv)
  }
}

# Dynamically create data frames from not_in_main
create_data_frames(not_in_main)
export(filtered_unique, "all_signature_expresson_filtred_protein.tsv")
# export(Phenotype_MSI_Genotype_CNV, "Target_genes_MSI_CNV.txt")
# export(Phenotype_Stemness_Genotype_Methylation, "Target_genes_Stemness_Methylation.txt")
# export(Phenotype_Stemness_Genotype_mRNA, "Target_genes_Stemness_mRNA.txt")
