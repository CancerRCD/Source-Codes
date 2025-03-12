# Setting the working directory
setwd("")

# Load necessary packages
library(rio)
library(dplyr)

# Import the data
genes <- import("all_genes_duplicatas_excel.xlsx")

# Check the structure of the data frame to identify duplicate columns
str(genes)

# Group by gene name and concatenate pathways
genes_simplificado <- genes %>%
  group_by(name) %>%
  summarise(
    uid = first(uid),
    description = first(description),
    chromosome = first(chromosome),
    summary = first(summary),
    Category_Search = first(Category_Search),  # Using the only occurrence after manual removal
    Pathway_Categorization = paste(unique(Pathway_Categorization), collapse = "/")
  )

# Check the result
print(genes_simplificado)

# Export the simplified dataset
export(genes_simplificado, "genes_simplificado_model_optosis.tsv")
