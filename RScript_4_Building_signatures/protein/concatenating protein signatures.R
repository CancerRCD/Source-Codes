# Load necessary packages
library(dplyr)
library(purrr)
library(fs)
library(readr)
library(tidyr)
library(rio)

#### PART A ####
# Define the base path
base_path <- "C:/Users/Emamnuell/Desktop/panoptosis/genes/genes_24_optosis/codigos/5°_Radar_siganature"

# Identify the phenotypic variables and their respective paths
phenotypes <- c("MSI", "Stemness", "TMB")  # Phenotype prefixes
phenotypic_suffix <- "_PanCan"  # Common suffix for phenotypes
genotypes <- c("Protein")

# Function to find the .tsv file with the specific "0.1" suffix in a folder
get_specific_suffix_tsv <- function(path) {
  files <- dir_ls(path, regexp = "gene_cancer_correlation_network_.*\\.tsv$")
  if (length(files) == 0) return(NULL)
  specific_file <- files[stringr::str_detect(basename(files), "0\\.1\\.tsv$")]
  if (length(specific_file) == 0) return(NULL)
  return(specific_file)
}

# Function to import and annotate the data
import_and_label_data <- function(phenotype, genotype) {
  phenotype_dir = paste(phenotype, phenotypic_suffix, sep="")
  folder_path <- file.path(base_path, phenotype_dir, paste(genotype, "vs", phenotype, sep = "_"))
  result_path <- file.path(folder_path, "results")
  
  if (!dir.exists(result_path)) {
    return(NULL)  # Returns NULL if the directory does not exist
  }
  
  file_path <- get_specific_suffix_tsv(result_path)
  
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
setwd("C:/Users/Emamnuell/Desktop/panoptosis/genes/genes_24_optosis/codigos/4°_signature/all")
export(all_data, "big_table_contacetado_signatures_protein.tsv")
