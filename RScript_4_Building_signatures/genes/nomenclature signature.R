# Load necessary packages
library(dplyr)
library(purrr)
library(fs)
library(readr)
library(tidyr)
library(rio)

setwd("C:/Users/Emamnuell/Desktop/panoptosis/genes/genes_24_optosis/codigos/4°_signature/genes")

aggregated_df <- import("signature_genes.xlsx")

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
