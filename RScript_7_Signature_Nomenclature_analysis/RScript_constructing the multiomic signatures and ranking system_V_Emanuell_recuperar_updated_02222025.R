## RScript_constructing the multiomic signatures and ranking system.
## Validation rerun and amended for validation on 05/03/2025 to 07/03/2025 
## bind or merge columns of different dataframes by common column values
## Note: if there are duplicated values if the common columns, do not used merge because the way merge works scrambles the rows
## Emanuell Rodrigues de Souza, Higor Almeida Cordeiro Nogueira, Ana Beatriz Garcia and Enrique Medina-Acostai

# Define the list of required packages
# remove.packages("rio")
# install.packages("rio", dependencies = TRUE)

library(rio)
packageVersion("rio")

required_packages <- c(
  "dplyr", "stringr", "tidyr", "writexl", "rio", "ggplot2", "boot", 
  "openxlsx", "gganimate", "png", "gifski",
  "geomtextpath", "ggdist", "gghighlight", "ggiraph", "ggpubr", 
  "ggrepel", "ggstatsplot", "ggtext", "patchwork", "ggforce","lemon", "png",
  "gganimate", "plotly", "patchwork", "correctR", "openxlsx", "ggpattern", "forcats"
)

# Install ggpattern if not already installed
# if (!requireNamespace("ggpattern", quietly = TRUE)) {
#   install.packages("ggpattern", repos = "https://cinc.rud.is")

# Function to check, install, and load missing packages
install_and_load_packages <- function(packages) {
  # Find packages that are not installed
  missing_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
  
  # Install any missing packages
  if (length(missing_packages) > 0) {
    message("Installing missing packages: ", paste(missing_packages, collapse = ", "))
    install.packages(missing_packages)
  }
  
  # Load all packages
  sapply(packages, require, character.only = TRUE)
}

# Call the function to install and load packages
install_and_load_packages(required_packages)

# For any GitHub-hosted packages (if CRAN installation fails)
if (!"geomtextpath" %in% installed.packages()[, "Package"]) {
  remotes::install_github("clauswilke/geomtextpath")
}

###### cbinding HRC dataframes #######
setwd("D:/Pré-artigo 5-optosis model/Dataframes all metrics")

# Get the current working directory and print it
current_wd <- getwd()
print(current_wd)

# Set the working directory back to the original (if needed)
setwd(current_wd)

# List all the TSV files in the working directory
tsv_files <- list.files(pattern = "\\.tsv$")

# Read the first TSV file with the column names
df1 <- read.delim(tsv_files[1])
df2 <- read.delim(tsv_files[2])
df3 <- read.delim(tsv_files[3])
df4 <- read.delim(tsv_files[4])

# Order each dataframe alphabetically by the "Gene", "tumor" and "gene_var" columns
df1 <- df1[order(df1$Gene, df1$tumor, df1$gene_var), ]
df2 <- df2[order(df2$Gene, df2$tumor, df2$gene_var), ]
df3 <- df3[order(df3$Gene, df3$tumor, df3$gene_var), ]
df4 <- df4[order(df4$Gene, df4$tumor, df4$gene_var), ]

# Ensure the first three columns are the same across all dataframes

# Compare the first three columns of df1 with df2
identical_1_2 <- all(df1[, 1:3] == df2[, 1:3])

# Compare the first three columns of df1 with df3
identical_1_3 <- all(df1[, 1:3] == df3[, 1:3])

# Compare the first three columns of df1 with df4
identical_1_4 <- all(df1[, 1:3] == df4[, 1:3])

# Check if all comparisons are true
if(identical_1_2 && identical_1_3 && identical_1_4) {
  print("The first three columns are identical across all four dataframes.")
} else {
  print("The first three columns are NOT identical across all four dataframes.")
}

# Detailed comparison (optional, in case of differences)
if(!identical_1_2) {
  print("Differences found between df1 and df2 in the first three columns.")
}
if(!identical_1_3) {
  print("Differences found between df1 and df3 in the first three columns.")
}
if(!identical_1_4) {
  print("Differences found between df1 and df4 in the first three columns.")
}

# List of dataframes
dfs <- list(df1, df2, df3, df4)

# Apply the filtering logic to each dataframe in the list
filtered_dfs_NS_1 <- lapply(dfs, function(df) {
  df %>%
    filter(p_val > 0.05, worst_prognosis_group != "NS")
})

# Assign the filtered dataframes to individual variables
filtered_df1_NS_1 <- filtered_dfs_NS_1[[1]]
filtered_df2_NS_1 <- filtered_dfs_NS_1[[2]]
filtered_df3_NS_1 <- filtered_dfs_NS_1[[3]]
filtered_df4_NS_1 <- filtered_dfs_NS_1[[4]]

# Testing for possible inconsistencies - Apply the filtering logic to each dataframe in the list
filtered_dfs_NS_2 <- lapply(dfs, function(df) {
  df %>%
    filter(p_val < 0.05, worst_prognosis_group == "NS")
})

# Assign the filtered dataframes to individual variables
filtered_df1_NS_2 <- filtered_dfs_NS_2[[1]]
filtered_df2_NS_2 <- filtered_dfs_NS_2[[2]]
filtered_df3_NS_2 <- filtered_dfs_NS_2[[3]]
filtered_df4_NS_2 <- filtered_dfs_NS_2[[4]]

DFI <- df1
DSS <- df2
OS <- df3
PFI <- df4

DFI <- DFI %>% rename(status_DFI = status, p_val_DFI = p_val, log_rank_DFI = log_rank, 
                      worst_prognosis_group_DFI = worst_prognosis_group)

DSS <- DSS %>% rename(status_DSS = status, p_val_DSS = p_val, log_rank_DSS = log_rank, 
                      worst_prognosis_group_DSS = worst_prognosis_group)

OS <- OS %>% rename(status_OS = status, p_val_OS = p_val, log_rank_OS = log_rank, 
                    worst_prognosis_group_OS = worst_prognosis_group)

PFI <- PFI %>% rename(status_PFI = status, p_val_PFI = p_val, log_rank_PFI = log_rank, 
                      worst_prognosis_group_PFI = worst_prognosis_group)

# Unir as colunas de DSS na tabela OS
resultado <- OS %>%
  bind_cols(
    DSS %>%
      select(status_DSS, p_val_DSS, log_rank_DSS, worst_prognosis_group_DSS)
  )

# Unir as colunas de DFI na tabela resultado
resultado <- resultado %>%
  bind_cols(
    DFI %>%
      select(status_DFI, p_val_DFI, log_rank_DFI, worst_prognosis_group_DFI)
  )

# Unir as colunas de PFI na tabela resultado
resultado <- resultado %>%
  bind_cols(
    PFI %>%
      select(status_PFI, p_val_PFI, log_rank_PFI, worst_prognosis_group_PFI)
  )

# Remove parentheses only for simple solo values in the "Gene" column of df1
resultado$Gene <- resultado$Gene %>%
  gsub("^\\(([a-zA-Z0-9_-]+)\\)$", "\\1", .)

# Apply trimws() to remove leading and trailing spaces from all columns in 'resultado'
resultado <- resultado %>%
  mutate(across(everything(), ~trimws(.)))  # Applies trimws to all columns

# Para o dataframe 'unicos' by "Gene"
unicos_gene <- data.frame(Gene = unique(resultado$Gene)) %>%
  arrange(Gene)  # Ordena os valores alfabeticamente

# Para o dataframe 'unicos', incluindo apenas as variáveis "Gene" e "tumor"
unicos_gene_tumor <- resultado %>%
  select(Gene, tumor) %>%  # Seleciona as colunas Gene e tumor
  distinct() %>%  # Mantém apenas as combinações únicas
  arrange(Gene, tumor)  # Ordena os valores alfabeticamente

# Para o dataframe 'unicos', incluindo apenas as variáveis "Gene", "tumor" e "gene_var"
unicos_gene_tumor_genvar <- resultado %>%
  select(Gene, tumor, gene_var) %>%  # Seleciona as colunas Gene e tumor
  distinct() %>%  # Mantém apenas as combinações únicas
  arrange(Gene, tumor, gene_var)  # Ordena os valores alfabeticamente

# Identify and create a dataframe 'non_unique_gene_tumor_genvar' with the non-unique rows
non_unique_gene_tumor_genvar <- resultado %>%
  select(Gene, tumor, gene_var) %>%  # Select the columns Gene, tumor, and gene_var
  group_by(Gene, tumor, gene_var) %>%  # Group by the selected columns
  filter(n() > 1) %>%  # Filter groups with more than 1 instance (non-unique)
  arrange(Gene, tumor, gene_var)  # Sort alphabetically by Gene, tumor, and gene_var

# Apply trimws() to remove leading and trailing spaces from all columns in 'resultado'
resultado_trimmed <- resultado %>%
  mutate(across(everything(), ~trimws(.)))  # Applies trimws to all columns

resultado <- resultado_trimmed

# Identify and exclude rows that are identical across all columns, and save those rows to a new dataframe
duplicate_rows <- resultado_trimmed %>%
  group_by(across(everything())) %>%  # Group by all columns
  filter(n() > 1) %>%  # Keep only the rows that appear more than once (duplicates)
  ungroup()  # Ungroup the dataframe after filtering

# Create a new dataframe excluding the duplicate rows
resultado_unique <- resultado_trimmed %>%
  distinct()  # Keep only unique rows across all columns

setwd("D:/Pré-artigo 5-optosis model/Dataframes all metrics")

rio::export(duplicate_rows, "duplicated_rows.tsv")

rio::export(resultado, "HRC_survival_analysis_all_signature.tsv")

###### Merging microenvinronment and TIL dataframes ###############
all_signature <- import("HRC_survival_analysis_all_signature.tsv")

til_mir <- import("TIL/TIL_results_miRNA.tsv")
til_mrna <- import("TIL/TIL_results_mRNA.tsv")
til_transcript <- import("TIL/TIL_results_transcript.tsv")
til_protein <- import("TIL/TIL_results_protein.tsv")

til_mir <- til_mir %>% rename(tumor = Cancer_type, Gene = miRNA_signature, immune_classification = classification,
                              immune_score_details = score_details, microenvironment_classification = Classification,
                              microenviroment_score_details = Score_Details, )


til_mrna <- til_mrna %>% rename(tumor = Cancer_type, Gene = Gene_signature, immune_classification = classification, 
                                immune_score_details = score_details, microenvironment_classification = Classification, 
                              microenviroment_score_details = Score_Details, )

til_transcript <- til_transcript %>% rename(tumor = Cancer_type, Gene = Transcript, immune_classification = classification, 
                                            immune_score_details = score_details, microenvironment_classification = Classification, 
                              microenviroment_score_details = Score_Details, )

til_protein <- til_protein %>% rename(tumor = Cancer_type, Gene = Protein_signature, immune_classification = classification, 
                                            immune_score_details = score_details, microenvironment_classification = Classification, 
                                            microenviroment_score_details = Score_Details, )

# Perform pairwise comparison of column names and their order for til_mir and til_mrna
if (identical(colnames(til_mir), colnames(til_mrna))) {
  print("The column names and their order are identical between til_mir and til_mrna.")
} else {
  print("The column names or their order are different between til_mir and til_mrna.")
  # Print the differences
  print("Columns in til_mir but not in til_mrna:")
  print(setdiff(colnames(til_mir), colnames(til_mrna)))
  
  print("Columns in til_mrna but not in til_mir:")
  print(setdiff(colnames(til_mrna), colnames(til_mir)))
}

# Perform pairwise comparison of column names and their order for til_mir and til_transcript
if (identical(colnames(til_mir), colnames(til_transcript))) {
  print("The column names and their order are identical between til_mir and til_transcript.")
} else {
  print("The column names or their order are different between til_mir and til_transcript.")
  # Print the differences
  print("Columns in til_mir but not in til_transcript:")
  print(setdiff(colnames(til_mir), colnames(til_transcript)))
  
  print("Columns in til_transcript but not in til_mir:")
  print(setdiff(colnames(til_transcript), colnames(til_mir)))
}

# Perform pairwise comparison of column names and their order for til_mir and til_protein
if (identical(colnames(til_mir), colnames(til_protein))) {
  print("The column names and their order are identical between til_mir and til_protein.")
} else {
  print("The column names or their order are different between til_mir and til_protein.")
  # Print the differences
  print("Columns in til_mir but not in til_protein:")
  print(setdiff(colnames(til_mir), colnames(til_protein)))
  
  print("Columns in til_protein but not in til_mir:")
  print(setdiff(colnames(til_protein), colnames(til_mir)))
}

# Perform pairwise comparison of column names and their order for til_mrna and til_transcript
if (identical(colnames(til_mrna), colnames(til_transcript))) {
  print("The column names and their order are identical between til_mrna and til_transcript.")
} else {
  print("The column names or their order are different between til_mrna and til_transcript.")
  # Print the differences
  print("Columns in til_mrna but not in til_transcript:")
  print(setdiff(colnames(til_mrna), colnames(til_transcript)))
  
  print("Columns in til_transcript but not in til_mrna:")
  print(setdiff(colnames(til_transcript), colnames(til_mrna)))
}

# Perform pairwise comparison of column names and their order for til_mrna and til_protein
if (identical(colnames(til_mrna), colnames(til_protein))) {
  print("The column names and their order are identical between til_mrna and til_protein.")
} else {
  print("The column names or their order are different between til_mrna and til_protein.")
  # Print the differences
  print("Columns in til_mrna but not in til_protein:")
  print(setdiff(colnames(til_mrna), colnames(til_protein)))
  
  print("Columns in til_protein but not in til_mrna:")
  print(setdiff(colnames(til_protein), colnames(til_mrna)))
}

# Perform pairwise comparison of column names and their order for til_transcript and til_protein
if (identical(colnames(til_transcript), colnames(til_protein))) {
  print("The column names and their order are identical between til_transcript and til_protein.")
} else {
  print("The column names or their order are different between til_transcript and til_protein.")
  # Print the differences
  print("Columns in til_transcript but not in til_protein:")
  print(setdiff(colnames(til_transcript), colnames(til_protein)))
  
  print("Columns in til_protein but not in til_transcript:")
  print(setdiff(colnames(til_protein), colnames(til_transcript)))
}

#### Deleting the old "protein" symbols signatures (n=97)
df80 <- resultado
sum(df80$gene_var == "protein") # old 97 rows; new 288
df81 <- til_mrna

# Assuming df80 and df81 are your dataframes
df82 <- merge(df80, df81, by = c("Gene", "tumor"), all = FALSE)

# Assuming df82 has a column named "gene_var"
df83 <- df82[df82$gene_var == "protein", ]

# Assuming df81 and df83 are your dataframes
df81_filtered <- df81[!paste(df81$Gene, df81$tumor) %in% paste(df83$Gene, df83$tumor), ]

til_mrna <- df81_filtered

##########################
##########################
##########################
all_signature_til <- rbind(til_mir, til_mrna, til_transcript,til_protein)

# Apply trimws() to remove leading and trailing spaces from all columns in 'all_signature_til'
all_signature_til_trimmed <- all_signature_til %>%
  mutate(across(everything(), ~trimws(.)))  # Applies trimws to all columns

all_signature_til <- all_signature_til_trimmed

# Identify and exclude rows that are identical across all columns, and save those rows to a new dataframe
duplicate_rows_til <- all_signature_til_trimmed %>%
  group_by(across(everything())) %>%  # Group by all columns
  filter(n() > 1) %>%  # Keep only the rows that appear more than once (duplicates)
  ungroup()  # Ungroup the dataframe after filtering

# Create a new dataframe excluding the duplicate rows
all_signature_til_unique <- all_signature_til_trimmed %>%
  distinct()  # Keep only unique rows across all columns

rio::export(duplicate_rows_til, "duplicated_rows_TIL.tsv")

# Para o dataframe 'unicos2', incluindo apenas as variáveis "Gene" e "tumor"
unicos_gene_2 <- all_signature_til %>%
  select(Gene) %>%  # Seleciona as colunas Gene e tumor
  distinct() %>%  # Mantém apenas as combinações únicas
  arrange(Gene)  # Ordena os valores alfabeticamente

unicos_gene_tumor_2 <- all_signature_til %>%
  select(Gene, tumor) %>%  # Seleciona as colunas Gene e tumor
  distinct() %>%  # Mantém apenas as combinações únicas
  arrange(Gene, tumor)  # Ordena os valores alfabeticamente

# Order alphabetically by Gene in unique_gene_1
unicos_gene <- unicos_gene %>%
  arrange(Gene)

# Order alphabetically by Gene in unique_gene_2
unicos_gene_2 <- unicos_gene_2 %>%
  arrange(Gene)

# Verificar se todos os valores na coluna 'Gene' são idênticos em ambos os dataframes
identical_genes <- identical(unicos_gene$Gene, unicos_gene_2$Gene)

# Verificar se todos os valores nas colunas 'Gene' e 'tumor' são idênticos em ambos os dataframes
identical_genes_tumores <- identical(unicos_gene_tumor[, c("Gene", "tumor")], unicos_gene_tumor_2[, c("Gene", "tumor")])

# Exibir o resultado
if (identical_genes) {
  print("Todos os valores na coluna 'Gene' são idênticos em ambos os dataframes.")
} else {
  print("Os valores na coluna 'Gene' NÃO são idênticos em ambos os dataframes.")
}

# Exibir o resultado
if (identical_genes_tumores) {
  print("Os valores de 'Gene' e 'tumor' são idênticos em ambos os dataframes.")
} else {
  print("Os valores de 'Gene' e 'tumor' são diferentes entre os dataframes.")
}

# Para verificar quais valores são diferentes (se houver)
differences <- setdiff(unicos_gene$Gene, unicos_gene_2$Gene)
if (length(differences) > 0) {
  print("Os seguintes valores estão presentes em 'unicos', mas não em 'unicos2':")
  print(differences)
} else {
  print("Não há diferenças entre os valores de 'Gene' nos dois dataframes but there may differences in row number.")
}

### Note: RICTOR is present in "all_signature" associated with mutation, 
### but not in "all_signature_til"
# Check for duplicates in unicos_gene
duplicates_unicos_gene <- unicos_gene %>%
  group_by(Gene) %>%
  filter(n() > 1)

# Check for duplicates in unicos_gene_2
duplicates_unicos_gene_2 <- unicos_gene_2 %>%
  group_by(Gene) %>%
  filter(n() > 1)

# Find genes in unicos_gene that are not in unicos_gene_2
genes_only_in_unicos_gene <- setdiff(unicos_gene$Gene, unicos_gene_2$Gene)

# Find genes in unicos_gene_2 that are not in unicos_gene
genes_only_in_unicos_gene_2 <- setdiff(unicos_gene_2$Gene, unicos_gene$Gene)

# Display the results
print("Genes in unicos_gene but not in unicos_gene_2:")
genes_only_in_unicos_gene

print("Genes in unicos_gene_2 but not in unicos_gene:")
genes_only_in_unicos_gene_2

# Merge selected columns from all_signatures_til into all_signatures
merged_data <- all_signature %>%
  left_join(all_signature_til %>% 
              select(Gene, tumor, immune_classification, immune_score_details, 
                     microenvironment_classification, microenviroment_score_details),
            by = c("Gene", "tumor"))

# Rename columns that are misspelled
colnames(merged_data) <- gsub("microenviroment", "microenvironment", colnames(merged_data))

# Replace empty values in the "immune_classification" column with "NS"
merged_data$immune_classification[merged_data$immune_classification == ""] <- "No data"

# Replace "pró-tumoral" values with "pro-tumoral" in the "microenvironment_classification" column
merged_data$microenvironment_classification[merged_data$microenvironment_classification == "pró-tumoral"] <- "pro-tumoral"

# Replace "dualista" values with "dual" in the "microenvironment_classification" column
merged_data$microenvironment_classification[merged_data$microenvironment_classification == "dualista"] <- "dual"

# Replace multiple values in the "immune_classification" column
merged_data$immune_classification <- gsub("Frio", "Cold", merged_data$immune_classification)
merged_data$immune_classification <- gsub("Quente", "Hot", merged_data$immune_classification)
merged_data$immune_classification <- gsub("Variável", "Variable", merged_data$immune_classification)

# View the updated column names
colnames(merged_data)

rio::export(merged_data, "til_survival_analysis_all_signature.tsv")

TIL_signature <- import("til_survival_analysis_all_signature.tsv")

###### TIL dataframe treatment ###################
# Filter the rows in TIL_signature that are duplicated by "Gene" and "tumor" variable values
duplicated_signatures <- TIL_signature %>%
  group_by(Gene, tumor) %>%
  filter(n() > 1) %>%
  ungroup()  # To remove the grouping structure

# Remove only one duplicated row for each duplicated combination of "Gene" and "tumor"
reduced_duplicates <- duplicated_signatures %>%
  group_by(Gene, tumor) %>%
  slice(-1) %>%  # Removes the first duplicate, keeping the rest
  ungroup()

# Filter rows that are duplicates across all columns in TIL_signature
duplicated_signatures_accros <- TIL_signature %>%
  filter(duplicated(.))  # Filter rows that are duplicates across all columns

###### Merging COX_signature and treatment ########
COX_signature <- import("COX df/result_cox_combined.tsv")

# Relocate the "Genotype" column after the "Cancer_type" column in the "COX_signature" data frame
COX_signature <- COX_signature %>%
  relocate(Genotype, .after = Cancer_type)

# Rename the column "signature" to "Signature" in the COX_signature data frame
COX_signature <- COX_signature %>%
  rename(Signature = signature)

# Exclude rows where the value in the "Genotype" column is "protein"
COX_signature_filtered <- COX_signature %>%
  filter(Genotype != "Protein")

COX_signature <- COX_signature_filtered

COX_signature_protein <- import("COX df/result_cox_combined_protein.xlsx")

# Rename the column "signature" to "Signature" in the COX_signature_protein data frame
COX_signature_protein <- COX_signature_protein %>%
  rename(Signature = Protein_signature)

COX_signature_protein$Signature <- COX_signature_protein$Signature %>%
  gsub("^\\(([a-zA-Z0-9_-]+)\\)$", "\\1", .)

# Relocate the "Genotype" column after the "Cancer_type" column in the "COX_signature_protein" data frame
COX_signature_protein <- COX_signature_protein %>%
  relocate(Genotype, .after = Cancer_type)

# Exclude specific columns from COX_signature_protein
COX_signature_protein_filtered <- COX_signature_protein %>%
  select(-Signature_Position, -RCD, -Phenotype, -Expression, -Correlation_sign,-Members)

COX_signature_protein <- COX_signature_protein_filtered

# Verify if the column names in COX_signature and COX_signature_protein are identical and in the same order
columns_identical <- identical(names(COX_signature), names(COX_signature_protein))

# Print the result
if (columns_identical) {
  print("The column names are identical and in the same order.")
} else {
  print("The column names are NOT identical or are in a different order.")
}

# Combine COX_signature and COX_signature_protein using rbind
COX_signature_combined <- rbind(COX_signature, COX_signature_protein)

COX_signature <-  COX_signature_combined 

COX_signature$Signature <- COX_signature$Signature %>%
  gsub("^\\(([a-zA-Z0-9_-]+)\\)$", "\\1", .)

# Rename multiple columns in TIL_signature
TIL_signature <- TIL_signature %>%
  rename(Signature = Gene, Cancer_type = tumor, Genotype = gene_var)

# Remove the specified columns from the TIL_signature data frame
TIL_signature <- TIL_signature %>%
  select(-status_OS, -status_DSS, -status_DFI, -status_PFI)

# Rename multiple columns in TIL_signature
TIL_signature <- TIL_signature %>%
  rename(Type_log_rank_OS = worst_prognosis_group_OS, Type_log_rank_DSS = worst_prognosis_group_DSS, 
         Type_log_rank_DFI = worst_prognosis_group_DFI, 
         Type_log_rank_PFI = worst_prognosis_group_PFI)

# Rename multiple columns in TIL_signature
TIL_signature <- TIL_signature %>%
  rename(p.value_log_rank_OS = p_val_OS, p.value_log_rank_DSS = p_val_DSS, 
         p.value_log_rank_DFI = p_val_DFI, 
         p.value_log_rank_PFI = p_val_PFI)

# Rename values under the "Genotype" column
TIL_signature <- TIL_signature %>%
  mutate(Genotype = recode(Genotype, 
                           "mRNA" = "mRNA", 
                           "methylation" = "Methylation", 
                           "mutation" = "Mutation", 
                           "miRNA" = "miRNA", 
                           "cnv" = "CNV", 
                           "transcript" = "Transcript", 
                           "protein" = "Protein"))

# Rename only the value "transcript" to "Transcript" in the "Genotype" column
COX_signature <- COX_signature %>%
  mutate(Genotype = recode(Genotype, "transcript" = "Transcript"))

# Use gsub to remove the parentheses and backticks only for values with a single "hsa-miR-..." or "hsa-let-..." gene
COX_signature$Signature <- gsub("^\\(`(hsa-(miR|let)-[^`]+)`\\)$", "\\1", COX_signature$Signature)

# Use gsub to remove the parentheses around "Signature" values that start with "ENS" and are enclosed in parentheses
COX_signature$Signature <- gsub("^\\((ENST\\w+)\\)$", "\\1", COX_signature$Signature)

# Ensure the first three columns are the same across all dataframes

# Sort both data frames by "Signature", "Cancer_type", and "Genotype" as before
TIL_signature <- TIL_signature %>%
  arrange(Signature, Cancer_type, Genotype)

COX_signature <- COX_signature %>%
  arrange(Signature, Cancer_type, Genotype)

unique_COX_Sig <- unique(COX_signature$Signature)
unique_TIL_Sig <- unique(TIL_signature$Signature)

unique_COX_Ct <- unique(COX_signature$Cancer_type)
unique_TIL_Ct <- unique(TIL_signature$Cancer_type)

unique_COX_Gt<- unique(COX_signature$Genotype)
unique_TIL_Gt <- unique(TIL_signature$Genotype)

# Check if any values in the "Signature" column of COX_signature contain "/" or "\"
contains_special_chars_COX <- grepl("[/\\\\]", COX_signature$Signature)

# Find rows where "/" or "\" is found in the "Signature" column of COX_signature
rows_with_special_chars_COX <- which(contains_special_chars_COX)

# Display the rows in COX_signature where "/" or "\" is found in the "Signature" column
if (length(rows_with_special_chars_COX) > 0) {
  print("Values containing '/' or '\\' found in the 'Signature' column of COX_signature:")
  print(COX_signature[rows_with_special_chars_COX, "Signature"])
} else {
  print("No values containing '/' or '\\' found in the 'Signature' column of COX_signature.")
}

# Check if any values in the "Signature" column of TIL_signature contain "/" or "\"
contains_special_chars_TIL <- grepl("[/\\\\]", TIL_signature$Signature)

# Find rows where "/" or "\" is found in the "Signature" column of TIL_signature
rows_with_special_chars_TIL <- which(contains_special_chars_TIL)

# Display the rows in TIL_signature where "/" or "\" is found in the "Signature" column
if (length(rows_with_special_chars_TIL) > 0) {
  print("Values containing '/' or '\\' found in the 'Signature' column of TIL_signature:")
  print(TIL_signature[rows_with_special_chars_TIL, "Signature"])
} else {
  print("No values containing '/' or '\\' found in the 'Signature' column of TIL_signature.")
}

# Use gsub to rename the values in the "Signature" column of COX_signature
COX_signature$Signature <- gsub("ago/02", "AGO2", COX_signature$Signature)
COX_signature$Signature <- gsub("ago/04", "AGO4", COX_signature$Signature)

# Ensure both data frames are ordered by "Signature", "Cancer_type", and "Genotype"
TIL_signature <- TIL_signature[order(TIL_signature$Signature, TIL_signature$Cancer_type, TIL_signature$Genotype), ]
COX_signature <- COX_signature[order(COX_signature$Signature, COX_signature$Cancer_type, COX_signature$Genotype), ]

# Apply trimws() to remove leading and trailing spaces from all columns 
TIL_signature <- TIL_signature %>%
  mutate(across(everything(), ~trimws(.)))  # Applies trimws to all columns

# Apply trimws() to remove leading and trailing spaces from all columns 
COX_signature <- COX_signature %>%
  mutate(across(everything(), ~trimws(.)))  # Applies trimws to all columns

# Determine the minimum number of rows to compare
min_rows <- min(nrow(COX_signature), nrow(TIL_signature))

# Compare the first three columns, up to the minimum number of rows
identical_COX_TIL <- all(COX_signature[1:min_rows, 1:3] == TIL_signature[1:min_rows, 1:3])

# Print result
if (identical_COX_TIL) {
  print("The first three columns are identical for the shared rows.")
} else {
  print("The first three columns differ for the shared rows.")
}

# Check for duplicates in the COX_signature data frame based on the three common columns
duplicates_COX <- COX_signature[duplicated(COX_signature[, c("Signature", "Cancer_type", "Genotype")]), ]
if (nrow(duplicates_COX) > 0) {
  print("Duplicates found in COX_signature!")
  } else {
  print("No duplicates in COX_signature.")
}

# Check for duplicates in the TIL_signature data frame based on the three common columns
duplicates_TIL <- TIL_signature[duplicated(TIL_signature[, c("Signature", "Cancer_type", "Genotype")]), ]
if (nrow(duplicates_TIL) > 0) {
  print("Duplicates found in TIL_signature!")
  } else {
  print("No duplicates in TIL_signature.")
}

# Even if the dimensions are different, join by the common columns
# Perform a full join to merge COX_signature and TIL_signature
merged_data_COX_TIL <- full_join(
  TIL_signature,
  COX_signature,
  by = c("Signature", "Cancer_type", "Genotype")
)

HRC_TMC_TIC_COX_signatures <- merged_data_COX_TIL

# Relocate the specified columns after "Type_Cox_PFI"
HRC_TMC_TIC_COX_signatures <- HRC_TMC_TIC_COX_signatures %>%
  relocate(microenvironment_classification, microenvironment_score_details,
           immune_classification, immune_score_details, .after = Type_Cox_PFI)

rio::export(HRC_TMC_TIC_COX_signatures, "HRC_TMC_TIC_COX_all_signatures.tsv")

##### Adding corresponding codes to values in variables ###########
##### Genotype to GFC
HRC_TMC_TIC_COX_all_signatures <-  import("HRC_TMC_TIC_COX_all_signatures.tsv")

# Define the correspondence table where the new column is named "GFC"
genotype_to_gfc <- data.frame(
  Genotype = c("Protein", "Mutation", "CNV", "miRNA", "Transcript", "mRNA", "Methylation"),
  GFC = c(1, 2, 3, 4, 5, 6, 7)  # The new column is defined as "GFC" here
)

# Add the new "GFC" column to HRC_TMC_TIC_COX_all_signatures
HRC_TMC_TIC_COX_all_signatures <- HRC_TMC_TIC_COX_all_signatures %>%
  left_join(genotype_to_gfc, by = "Genotype") %>%  # "GFC" gets added based on the "Genotype"
  relocate(GFC, .after = Genotype)  # Ensure "GFC" is placed after "Genotype"

##### Replacing Type_Cox_X to A, B or C, where A signifies no effect or no data, B indicates risk, and C denotes protection

# Define the columns to replace empty values with "NS"
columns_to_replace <- c("Type_Cox_OS", "Type_Cox_DSS", "Type_Cox_DFI", "Type_Cox_PFI")

# Replace empty "" values in the specified columns with "NS"
HRC_TMC_TIC_COX_all_signatures <- HRC_TMC_TIC_COX_all_signatures %>%
  mutate_at(vars(one_of(columns_to_replace)), ~ ifelse(. == "", "NS", .))

# Define the correspondence table where the new column is named "CoxI_OS"
Type_Cox_OS_to_CoxI_OS <- data.frame(
  Type_Cox_OS = c("NS", "Risky", "Protective"),  # Fix the column name and values
  CoxI_OS = c("A", "B", "C")  # Values for the new "CoxI_OS" column (using strings here)
)

# Add the new "CoxI_OS" column to HRC_TMC_TIC_COX_all_signatures
HRC_TMC_TIC_COX_all_signatures <- HRC_TMC_TIC_COX_all_signatures %>%
  left_join(Type_Cox_OS_to_CoxI_OS, by = "Type_Cox_OS") %>%  # Join based on "Type_Cox_OS"
  relocate(CoxI_OS, .after = Type_Cox_OS)  # Move "CoxI_OS" after "Type_Cox_OS"

# Define the correspondence table where the new column is named "CoxI_DSS"
Type_Cox_DSS_to_CoxI_DSS <- data.frame(
  Type_Cox_DSS = c("NS", "Risky", "Protective"),  # Fix the column name and values
  CoxI_DSS = c("A", "B", "C")  # Values for the new "CoxI_DSS" column (using strings here)
)

# Add the new "CoxI_DSS" column to HRC_TMC_TIC_COX_all_signatures
HRC_TMC_TIC_COX_all_signatures <- HRC_TMC_TIC_COX_all_signatures %>%
  left_join(Type_Cox_DSS_to_CoxI_DSS, by = "Type_Cox_DSS") %>%  # Join based on "Type_Cox_DSS"
  relocate(CoxI_DSS, .after = Type_Cox_DSS)  # Move "CoxI_DSS" after "Type_Cox_DSS"

# Define the correspondence table where the new column is named "CoxI_DFI"
Type_Cox_DFI_to_CoxI_DFI <- data.frame(
  Type_Cox_DFI = c("NS", "Risky", "Protective"),  # Fix the column name and values
  CoxI_DFI = c("A", "B", "C")  # Values for the new "CoxI_DFI" column (using strings here)
)

# Add the new "CoxI_DFI" column to HRC_TMC_TIC_COX_all_signatures
HRC_TMC_TIC_COX_all_signatures <- HRC_TMC_TIC_COX_all_signatures %>%
  left_join(Type_Cox_DFI_to_CoxI_DFI, by = "Type_Cox_DFI") %>%  # Join based on "Type_Cox_DFI"
  relocate(CoxI_DFI, .after = Type_Cox_DFI)  # Move "CoxI_DFI" after "Type_Cox_DFI"

# Define the correspondence table where the new column is named "CoxI_PFI"
Type_Cox_PFI_to_CoxI_PFI <- data.frame(
  Type_Cox_PFI = c("NS", "Risky", "Protective"),  # Fix the column name and values
  CoxI_PFI = c("A", "B", "C")  # Values for the new "CoxI_PFI" column (using strings here)
)

# Add the new "CoxI_PFI" column to HRC_TMC_TIC_COX_all_signatures
HRC_TMC_TIC_COX_all_signatures <- HRC_TMC_TIC_COX_all_signatures %>%
  left_join(Type_Cox_PFI_to_CoxI_PFI, by = "Type_Cox_PFI") %>%  # Join based on "Type_Cox_PFI"
  relocate(CoxI_PFI, .after = Type_Cox_PFI)  # Move "CoxI_PFI" after "Type_Cox_PFI"

########
# Create and populate new columns, then relocate them according to the template
HRC_TMC_TIC_COX_all_signatures <- HRC_TMC_TIC_COX_all_signatures %>%
  
  # Create and populate new columns
  mutate(Cox_DSS_id = 1,
         Cox_DFI_id = 2,
         Cox_PFI_id = 3,
         Cox_OS_id = 4,
         logrank_DSS_id = 1,
         logrank_DFI_id = 2,
         logrank_PFI_id = 3,
         logrank_OS_id = 4) %>%
  
  # Relocate new columns according to the specified template
  relocate(Cox_DSS_id, .after = p.value_Cox_DSS) %>%
  relocate(Cox_DFI_id, .after = p.value_Cox_DFI) %>%
  relocate(Cox_PFI_id, .after = p.value_Cox_PFI) %>%
  relocate(Cox_OS_id, .after = p.value_Cox_OS) %>%
  relocate(logrank_DSS_id, .after = p.value_log_rank_DSS) %>%
  relocate(logrank_DFI_id, .after = p.value_log_rank_DFI) %>%
  relocate(logrank_PFI_id, .after = p.value_log_rank_PFI) %>%
  relocate(logrank_OS_id, .after = p.value_log_rank_OS)

# Define the columns to replace "No data" or "" values with "NS"
columns_to_replace <- c("immune_classification", "microenvironment_classification", "immune_score_details", "microenvironment_score_details")

# Replace "No data" or "" values in the specified columns with "NS"
HRC_TMC_TIC_COX_all_signatures <- HRC_TMC_TIC_COX_all_signatures %>%
  mutate_at(vars(one_of(columns_to_replace)), ~ ifelse(. %in% c("No data", ""), "NS", .))

# Define the mappings for immune_classification to TIL_ID
immune_mapping <- data.frame(
  immune_classification = c("Hot", "Variable", "Cold", "NS"),
  TIL_ID = c(1, 2, 3, 4)
)

# Define the mappings for microenvironment_classification to TMC_ID
microenvironment_mapping <- data.frame(
  microenvironment_classification = c("anti-tumoral", "dual", "pro-tumoral", "NS"),
  TMC_ID = c(1, 2, 3, 4)
)

# Create and populate the new columns in the specified dataframe
HRC_TMC_TIC_COX_all_signatures <- HRC_TMC_TIC_COX_all_signatures %>%
  # Join the immune classification mapping
  left_join(immune_mapping, by = "immune_classification") %>%
  relocate(TIL_ID, .after = immune_classification) %>%  # Place TIL_ID after immune_classification
  
  # Join the microenvironment classification mapping
  left_join(microenvironment_mapping, by = "microenvironment_classification") %>%
  relocate(TMC_ID, .after = microenvironment_classification)  # Place TMC_ID after microenvironment_classification

rio::export(HRC_TMC_TIC_COX_all_signatures, "HRC_TMC_TIC_COX_all_signatures_with_IDs.tsv")

#####################

HRC_TMC_TIC_COX_all_signatures_with_IDs <-  import("HRC_TMC_TIC_COX_all_signatures_with_IDs.tsv")

# Filter rows where Genotype is 'CNV' and return unique values based on the 'Type_log_rank_OS' column
filtered_data_CNV <- HRC_TMC_TIC_COX_all_signatures_with_IDs[
  HRC_TMC_TIC_COX_all_signatures_with_IDs$Genotype == "CNV" &
    HRC_TMC_TIC_COX_all_signatures_with_IDs$Type_log_rank_OS %in% unique(HRC_TMC_TIC_COX_all_signatures_with_IDs$Type_log_rank_OS), 
]

# Select column names that match the structure "Type_XXX"
type_columns <- grep("^Type_\\w+", colnames(HRC_TMC_TIC_COX_all_signatures_with_IDs), value = TRUE)

# Print the selected column names in the console
cat("Columns with 'Type_XXX' structure:\n", paste(type_columns, collapse = "\n"))

# Modify the dataframe by creating new columns based on existing ones and mapping the values
HRC_TMC_TIC_COX_all_signatures_with_IDs <- HRC_TMC_TIC_COX_all_signatures_with_IDs %>%
  mutate(
    # Populate Type_log_rank_OS_I based on Genotype and Type_log_rank_OS
    Type_log_rank_OS_I = case_when(
      Genotype == "Protein" & Type_log_rank_OS == "NS" ~ "A",
      Genotype == "Protein" & Type_log_rank_OS == "High" ~ "B",
      Genotype == "Protein" & Type_log_rank_OS == "Low" ~ "C",
      
      Genotype == "Mutation" & Type_log_rank_OS == "NS" ~ "A",
      Genotype == "Mutation" & Type_log_rank_OS == "MT" ~ "B",
      Genotype == "Mutation" & Type_log_rank_OS == "WT" ~ "C",
      
      Genotype == "CNV" & Type_log_rank_OS == "NS" ~ "A",
      Genotype == "CNV" & Type_log_rank_OS == "Deleted" ~ "B",
      Genotype == "CNV" & Type_log_rank_OS == "Duplicated" ~ "C",
      Genotype == "CNV" & Type_log_rank_OS == "Deleted/Duplicated" ~ "D",
      
      Genotype == "miRNA" & Type_log_rank_OS == "NS" ~ "A",
      Genotype == "miRNA" & Type_log_rank_OS == "High" ~ "B",
      Genotype == "miRNA" & Type_log_rank_OS == "Low" ~ "C",
      
      Genotype == "Transcript" & Type_log_rank_OS == "NS" ~ "A",
      Genotype == "Transcript" & Type_log_rank_OS == "High" ~ "B",
      Genotype == "Transcript" & Type_log_rank_OS == "Low" ~ "C",
      
      Genotype == "mRNA" & Type_log_rank_OS == "NS" ~ "A",
      Genotype == "mRNA" & Type_log_rank_OS == "High" ~ "B",
      Genotype == "mRNA" & Type_log_rank_OS == "Low" ~ "C",
      
      Genotype == "Methylation" & Type_log_rank_OS == "NS" ~ "A",
      Genotype == "Methylation" & Type_log_rank_OS == "High" ~ "B",
      Genotype == "Methylation" & Type_log_rank_OS == "Low" ~ "C",
      TRUE ~ NA_character_
    ),
    
    # Populate Type_log_rank_DSS_I based on Genotype and Type_log_rank_DSS
    Type_log_rank_DSS_I = case_when(
      Genotype == "Protein" & Type_log_rank_DSS == "NS" ~ "A",
      Genotype == "Protein" & Type_log_rank_DSS == "High" ~ "B",
      Genotype == "Protein" & Type_log_rank_DSS == "Low" ~ "C",
      
      Genotype == "Mutation" & Type_log_rank_DSS == "NS" ~ "A",
      Genotype == "Mutation" & Type_log_rank_DSS == "MT" ~ "B",
      Genotype == "Mutation" & Type_log_rank_DSS == "WT" ~ "C",
      
      Genotype == "CNV" & Type_log_rank_DSS == "NS" ~ "A",
      Genotype == "CNV" & Type_log_rank_DSS == "Deleted" ~ "B",
      Genotype == "CNV" & Type_log_rank_DSS == "Duplicated" ~ "C",
      Genotype == "CNV" & Type_log_rank_DSS == "Deleted/Duplicated" ~ "D",
      
      Genotype == "miRNA" & Type_log_rank_DSS == "NS" ~ "A",
      Genotype == "miRNA" & Type_log_rank_DSS == "High" ~ "B",
      Genotype == "miRNA" & Type_log_rank_DSS == "Low" ~ "C",
      
      Genotype == "Transcript" & Type_log_rank_DSS == "NS" ~ "A",
      Genotype == "Transcript" & Type_log_rank_DSS == "High" ~ "B",
      Genotype == "Transcript" & Type_log_rank_DSS == "Low" ~ "C",
      
      Genotype == "mRNA" & Type_log_rank_DSS == "NS" ~ "A",
      Genotype == "mRNA" & Type_log_rank_DSS == "High" ~ "B",
      Genotype == "mRNA" & Type_log_rank_DSS == "Low" ~ "C",
      
      Genotype == "Methylation" & Type_log_rank_DSS == "NS" ~ "A",
      Genotype == "Methylation" & Type_log_rank_DSS == "High" ~ "B",
      Genotype == "Methylation" & Type_log_rank_DSS == "Low" ~ "C",
      TRUE ~ NA_character_
    ),
    
    # Populate Type_log_rank_DFI_I based on Genotype and Type_log_rank_DFI
    Type_log_rank_DFI_I = case_when(
      Genotype == "Protein" & Type_log_rank_DFI == "NS" ~ "A",
      Genotype == "Protein" & Type_log_rank_DFI == "High" ~ "B",
      Genotype == "Protein" & Type_log_rank_DFI == "Low" ~ "C",
      
      Genotype == "Mutation" & Type_log_rank_DFI == "NS" ~ "A",
      Genotype == "Mutation" & Type_log_rank_DFI == "MT" ~ "B",
      Genotype == "Mutation" & Type_log_rank_DFI == "WT" ~ "C",
      
      Genotype == "CNV" & Type_log_rank_DFI == "NS" ~ "A",
      Genotype == "CNV" & Type_log_rank_DFI == "Deleted" ~ "B",
      Genotype == "CNV" & Type_log_rank_DFI == "Duplicated" ~ "C",
      Genotype == "CNV" & Type_log_rank_DFI == "Deleted/Duplicated" ~ "D",
      
      Genotype == "miRNA" & Type_log_rank_DFI == "NS" ~ "A",
      Genotype == "miRNA" & Type_log_rank_DFI == "High" ~ "B",
      Genotype == "miRNA" & Type_log_rank_DFI == "Low" ~ "C",
      
      Genotype == "Transcript" & Type_log_rank_DFI == "NS" ~ "A",
      Genotype == "Transcript" & Type_log_rank_DFI == "High" ~ "B",
      Genotype == "Transcript" & Type_log_rank_DFI == "Low" ~ "C",
      
      Genotype == "mRNA" & Type_log_rank_DFI == "NS" ~ "A",
      Genotype == "mRNA" & Type_log_rank_DFI == "High" ~ "B",
      Genotype == "mRNA" & Type_log_rank_DFI == "Low" ~ "C",
      
      Genotype == "Methylation" & Type_log_rank_DFI == "NS" ~ "A",
      Genotype == "Methylation" & Type_log_rank_DFI == "High" ~ "B",
      Genotype == "Methylation" & Type_log_rank_DFI == "Low" ~ "C",
      TRUE ~ NA_character_
    ),
    
    # Populate Type_log_rank_PFI_I based on Genotype and Type_log_rank_PFI
    Type_log_rank_PFI_I = case_when(
      Genotype == "Protein" & Type_log_rank_PFI == "NS" ~ "A",
      Genotype == "Protein" & Type_log_rank_PFI == "High" ~ "B",
      Genotype == "Protein" & Type_log_rank_PFI == "Low" ~ "C",
      
      Genotype == "Mutation" & Type_log_rank_PFI == "NS" ~ "A",
      Genotype == "Mutation" & Type_log_rank_PFI == "MT" ~ "B",
      Genotype == "Mutation" & Type_log_rank_PFI == "WT" ~ "C",
      
      Genotype == "CNV" & Type_log_rank_PFI == "NS" ~ "A",
      Genotype == "CNV" & Type_log_rank_PFI == "Deleted" ~ "B",
      Genotype == "CNV" & Type_log_rank_PFI == "Duplicated" ~ "C",
      Genotype == "CNV" & Type_log_rank_PFI == "Deleted/Duplicated" ~ "D",
      
      Genotype == "miRNA" & Type_log_rank_PFI == "NS" ~ "A",
      Genotype == "miRNA" & Type_log_rank_PFI == "High" ~ "B",
      Genotype == "miRNA" & Type_log_rank_PFI == "Low" ~ "C",
      
      Genotype == "Transcript" & Type_log_rank_PFI == "NS" ~ "A",
      Genotype == "Transcript" & Type_log_rank_PFI == "High" ~ "B",
      Genotype == "Transcript" & Type_log_rank_PFI == "Low" ~ "C",
      
      Genotype == "mRNA" & Type_log_rank_PFI == "NS" ~ "A",
      Genotype == "mRNA" & Type_log_rank_PFI == "High" ~ "B",
      Genotype == "mRNA" & Type_log_rank_PFI == "Low" ~ "C",
      
      Genotype == "Methylation" & Type_log_rank_PFI == "NS" ~ "A",
      Genotype == "Methylation" & Type_log_rank_PFI == "High" ~ "B",
      Genotype == "Methylation" & Type_log_rank_PFI == "Low" ~ "C",
      TRUE ~ NA_character_
    )
  )

# Relocate the specified columns to be positioned after "Type_log_rank_PFI_I" in the correct order
HRC_TMC_TIC_COX_all_signatures_with_IDs <- HRC_TMC_TIC_COX_all_signatures_with_IDs %>%
  relocate(logrank_DSS_id, .after = immune_score_details) %>%
  relocate(Type_log_rank_DSS_I, .after = logrank_DSS_id) %>%
  relocate(logrank_DFI_id, .after = Type_log_rank_DSS_I) %>%
  relocate(Type_log_rank_DFI_I, .after = logrank_DFI_id) %>%
  relocate(logrank_PFI_id, .after = Type_log_rank_DFI_I) %>%
  relocate(Type_log_rank_PFI_I, .after = logrank_PFI_id) %>%
  relocate(logrank_OS_id, .after = Type_log_rank_PFI_I) %>%
  relocate(Type_log_rank_OS_I, .after = logrank_OS_id)

# Relocate the specified columns to be positioned after "Type_log_rank_PFI_I" in the correct order
HRC_TMC_TIC_COX_all_signatures_with_IDs <- HRC_TMC_TIC_COX_all_signatures_with_IDs %>%
  relocate(Cox_DSS_id, .after = Type_log_rank_OS_I) %>%
  relocate(CoxI_DSS, .after = Cox_DSS_id) %>%
  relocate(Cox_DFI_id, .after = CoxI_DSS) %>%
  relocate(CoxI_DFI, .after = Cox_DFI_id) %>%
  relocate(Cox_PFI_id, .after = CoxI_DFI) %>%
  relocate(CoxI_PFI, .after = Cox_PFI_id) %>%
  relocate(Cox_OS_id, .after = CoxI_PFI) %>%
  relocate(CoxI_OS, .after = Cox_OS_id)

# Count the number of empty "" or NA values in the entire dataframe
empty_na_count <- sum(is.na(HRC_TMC_TIC_COX_all_signatures_with_IDs) | HRC_TMC_TIC_COX_all_signatures_with_IDs == "")

print(paste("Number of empty or NA values: ", empty_na_count))

# Replace universaly "NA" and "" empty values with "NS" in the entire dataframe
HRC_TMC_TIC_COX_all_signatures_with_IDs <- HRC_TMC_TIC_COX_all_signatures_with_IDs %>%
  mutate_all(~ ifelse(is.na(.) | . == "", "NS", .))

# Replace empty "" values with "NS" in all columns of the dataframe
HRC_TMC_TIC_COX_all_signatures_with_IDs <- HRC_TMC_TIC_COX_all_signatures_with_IDs %>%
  mutate_all(~ ifelse(. == "", "NS", .))

rio::export(HRC_TMC_TIC_COX_all_signatures_with_IDs, "HRC_TMC_TIC_COX_all_signatures_with_IDs.tsv")

HRC_TMC_TIC_COX_all_signatures_with_IDs <-  import("HRC_TMC_TIC_COX_all_signatures_with_IDs.tsv")

##### Concatenate values from specific columns into a new column variable

# Select column names that match the structure "Type_XXX"
type_columns <- grep("^Type_log_rank_\\w+", colnames(HRC_TMC_TIC_COX_all_signatures_with_IDs), value = TRUE)

type_columns

# Select column names that match the structure "Type_XXX"
type_columns <- grep("^logrank_\\w+", colnames(HRC_TMC_TIC_COX_all_signatures_with_IDs), value = TRUE)

type_columns

# 1. Convert the specified columns to character type
HRC_TMC_TIC_COX_all_signatures_with_IDs <- HRC_TMC_TIC_COX_all_signatures_with_IDs %>%
  mutate_at(vars(logrank_DSS_id, log_rank_DSS, logrank_DFI_id, log_rank_DFI, 
                 logrank_PFI_id, log_rank_PFI, logrank_OS_id, log_rank_OS), 
            as.character)

# 2. Create a new column variable "SMC_series" by concatenating the character values of the specified columns
HRC_TMC_TIC_COX_all_signatures_with_IDs <- HRC_TMC_TIC_COX_all_signatures_with_IDs %>%
  mutate(SMC_series = paste0(logrank_DSS_id, Type_log_rank_DSS_I, logrank_DFI_id, Type_log_rank_DFI_I, 
                             logrank_PFI_id, Type_log_rank_PFI_I, logrank_OS_id, Type_log_rank_OS_I))

# Select column names that match the structure "Type_XXX"
type_columns <- grep("^Cox_\\w+", colnames(HRC_TMC_TIC_COX_all_signatures_with_IDs), value = TRUE)

type_columns

# Select column names that match the structure "Type_XXX"
type_columns <- grep("^CoxI_\\w+", colnames(HRC_TMC_TIC_COX_all_signatures_with_IDs), value = TRUE)

type_columns

# 1. Convert the specified columns to character type
HRC_TMC_TIC_COX_all_signatures_with_IDs <- HRC_TMC_TIC_COX_all_signatures_with_IDs %>%
  mutate_at(vars(Cox_DSS_id, CoxI_DSS, Cox_DFI_id, CoxI_DFI, 
                 Cox_PFI_id, CoxI_PFI, Cox_OS_id, CoxI_OS), 
            as.character)

# 2. Create a new column variable "HRC_series" by concatenating the character values of the specified columns
HRC_TMC_TIC_COX_all_signatures_with_IDs <- HRC_TMC_TIC_COX_all_signatures_with_IDs %>%
  mutate(HRC_series = paste0(Cox_DSS_id, CoxI_DSS, Cox_DFI_id, CoxI_DFI, 
                             Cox_PFI_id, CoxI_PFI, Cox_OS_id, CoxI_OS))

# Creating the valid (useful) universe of array combination in HRC and SMC arrays for mapping purposes 

# Get the unique values from HRC_series and SMC_series
Unique_HRC_array <- unique(HRC_TMC_TIC_COX_all_signatures_with_IDs$HRC_series)
Unique_SMC_array <- unique(HRC_TMC_TIC_COX_all_signatures_with_IDs$SMC_series)

# Convert unique values to data frames
df_HRC <- data.frame(series = Unique_HRC_array)  # Rename the column as "series"
df_SMC <- data.frame(series = Unique_SMC_array)  # Rename the column as "series"

# Use rbind to combine the two data frames
combined_HRC_SMC_array <- rbind(df_HRC, df_SMC)

# Filter out duplicated values in the combined data frame
unique_combined_HRC_SMC_array_filtered <- combined_HRC_SMC_array %>%
  distinct()

# Order the dataframe alphabetically by the "series" column
unique_combined_HRC_SMC_array_filtered <- unique_combined_HRC_SMC_array_filtered %>%
  arrange(series)

# Add a column "array_series_number" and populate it with a series from 0 to n-1
unique_combined_HRC_SMC_array_filtered <- unique_combined_HRC_SMC_array_filtered %>%
  mutate(array_series_number = 0:(n() - 1)) %>%  # Generate the series from 0 to n-1
  relocate(array_series_number, .before = series)  # Place the new column before "series"

rio::export(unique_combined_HRC_SMC_array_filtered, "unique_combined_HRC_SMC_array_filtered.tsv")

# Count the number of occurrences of "1NS2NS3NS4NS" in df1
count_value <- sum(HRC_TMC_TIC_COX_all_signatures_with_IDs == "1NS2NS3NS4NS", na.rm = TRUE)

unique_combined_HRC_SMC_array_filtered <- import("unique_combined_HRC_SMC_array_filtered.tsv")

#### Mapping the arrays in HRC_map
# Assuming the df "unique_combined_HRC_SMC_array_filtered" has the mapping template with columns:
# "array_series_number" and "series"
# And the df "HRC_TMC_TIC_COX_all_signatures_with_IDs" has the "HRC_series" column

# Perform the left join to add the mapping to "HRC_TMC_TIC_COX_all_signatures_with_IDs"
HRC_TMC_TIC_COX_all_signatures_with_IDs <- HRC_TMC_TIC_COX_all_signatures_with_IDs %>%
  left_join(unique_combined_HRC_SMC_array_filtered, 
            by = c("HRC_series" = "series")) %>%  # Match "HRC_series" with "series" in the mapping template
  rename(HRC_map = array_series_number) %>%  # Rename the "array_series_number" column to "HRC_map"
  relocate(HRC_map, .after = HRC_series)  # Place the new "HRC_map" column after "HRC_series"

#### Mapping the arrays in SMC_map
# Assuming the df "unique_combined_HRC_SMC_array_filtered" has the mapping template with columns:
# "array_series_number" and "series"
# And the df "HRC_TMC_TIC_COX_all_signatures_with_IDs" has the "SMC_series" column

# Perform the left join to add the mapping to "HRC_TMC_TIC_COX_all_signatures_with_IDs"
HRC_TMC_TIC_COX_all_signatures_with_IDs <- HRC_TMC_TIC_COX_all_signatures_with_IDs %>%
  left_join(unique_combined_HRC_SMC_array_filtered, 
            by = c("SMC_series" = "series")) %>%  # Match "SMC_series" with "series" in the mapping template
  rename(SMC_map = array_series_number) %>%  # Rename the "array_series_number" column to "SMC_map"
  relocate(SMC_map, .after = SMC_series)  # Place the new "SMC_map" column after "SMC_series"

rio::export(HRC_TMC_TIC_COX_all_signatures_with_IDs, "HRC_TMC_TIC_COX_all_signatures_with_IDs_arrays.tsv")

HRC_TMC_TIC_COX_all_signatures_with_IDs_arrays <- import("HRC_TMC_TIC_COX_all_signatures_with_IDs_arrays.tsv")

# Replace column names
colnames(HRC_TMC_TIC_COX_all_signatures_with_IDs_arrays)[colnames(HRC_TMC_TIC_COX_all_signatures_with_IDs_arrays) == "Cancer_type"] <- "CTAB"
colnames(HRC_TMC_TIC_COX_all_signatures_with_IDs_arrays)[colnames(HRC_TMC_TIC_COX_all_signatures_with_IDs_arrays) == "HRC_map"] <- "HRC"
colnames(HRC_TMC_TIC_COX_all_signatures_with_IDs_arrays)[colnames(HRC_TMC_TIC_COX_all_signatures_with_IDs_arrays) == "SMC_map"] <- "SMC"
colnames(HRC_TMC_TIC_COX_all_signatures_with_IDs_arrays)[colnames(HRC_TMC_TIC_COX_all_signatures_with_IDs_arrays) == "TMC_ID"] <- "TMC"
colnames(HRC_TMC_TIC_COX_all_signatures_with_IDs_arrays)[colnames(HRC_TMC_TIC_COX_all_signatures_with_IDs_arrays) == "TIL_ID"] <- "TIC"

### Order the columns as for the signature nomenclature:
# Get the current order of columns
cols <- colnames(HRC_TMC_TIC_COX_all_signatures_with_IDs_arrays)

# Remove the columns that will be relocated
cols_to_move <- c("GFC", "HRC", "SMC", "TMC", "TIC")
remaining_cols <- setdiff(cols, cols_to_move)

# Find the index of "CTAB" to determine the correct location to insert
ctab_index <- which(remaining_cols == "CTAB")

# Rebuild the columns in the desired order
new_order <- c(
  remaining_cols[1:ctab_index],   # Columns up to and including "CTAB"
  "GFC",                          # Sequentially place "GFC" after "GSI"
  "HRC",                          # Sequentially place "HRC" after "GSI"
  "SMC",                          # Sequentially place "SMC" after "HRC"
  "TMC",                          # Sequentially place "TMC" after "HRC"
  "TIC",                          # Sequentially place "TIC" after "TMC"
  remaining_cols[(ctab_index+1):length(remaining_cols)]  # Rest of the columns
)

# Apply the new column order
HRC_TMC_TIC_COX_all_signatures_with_IDs_arrays <- HRC_TMC_TIC_COX_all_signatures_with_IDs_arrays[, new_order]

# Get column names from both dataframes
cols_array <- colnames(HRC_TMC_TIC_COX_all_signatures_with_IDs_arrays)
cols_ids <- colnames(HRC_TMC_TIC_COX_all_signatures_with_IDs)

# Find the differences in column names
diff_in_array <- setdiff(cols_array, cols_ids)  # Columns present in 'HRC_TMC_TIC_COX_all_signatures_with_IDs_arrays' but not in 'HRC_TMC_TIC_COX_all_signatures_with_IDs'
diff_in_ids <- setdiff(cols_ids, cols_array)    # Columns present in 'HRC_TMC_TIC_COX_all_signatures_with_IDs' but not in 'HRC_TMC_TIC_COX_all_signatures_with_IDs_arrays'

# Print the differences
cat("Columns in 'HRC_TMC_TIC_COX_all_signatures_with_IDs_arrays' but not in 'HRC_TMC_TIC_COX_all_signatures_with_IDs':\n")
print(diff_in_array)

cat("\nColumns in 'HRC_TMC_TIC_COX_all_signatures_with_IDs' but not in 'HRC_TMC_TIC_COX_all_signatures_with_IDs_arrays':\n")
print(diff_in_ids)

rio::export(HRC_TMC_TIC_COX_all_signatures_with_IDs_arrays, "HRC_TMC_TIC_COX_all_signatures_with_IDs_arrays.tsv")

######
######
######
###### Adjusting dataframes with pending variables to be merged #####
######
######
######
# Importar os data frames
df1 <- import("drafts/df_combinado_NEW.tsv")

# # Exclude rows where the value in the "Genotype" column is "Protein"
# df1 <- df1 %>%
#   filter(Genotype != "Protein")

# setwd("D:/Pré-artigo 5-optosis model/Dataframes all metrics")
df2 <- import("HRC_TMC_TIC_COX_all_signatures_with_IDs_arrays.tsv")

# Remove leading and trailing spaces globally across all columns in df2_ordered
df1 <- df1 %>%
  mutate_all(trimws)

df2 <- df2 %>%
  mutate_all(trimws)

# Substituir "ago/02" por "AGO2" e "ago/04" por "AGO4" na coluna Signature de df1
df1 <- df1 %>%
  mutate(Signature = ifelse(Signature == "ago/02", "AGO2",
                            ifelse(Signature == "ago/04", "AGO4", Signature)))

df1 <- df1 %>%
  mutate(Genotype = gsub("transcript", "Transcript", Genotype))

# Remove parentheses only for simple solo values in the "Signature" column 
df1$Signature <- df1$Signature %>%
  gsub("^\\(([a-zA-Z0-9_-]+)\\)$", "\\1", .)

# 1. Order both dataframes alphabetically by "Signature", "CTAB", and "Genotype"
df1_ordered <- df1[order(df1$Signature, df1$CTAB, df1$Genotype), ]
df2_ordered <- df2[order(df2$Signature, df2$CTAB, df2$Genotype), ]

# 2. Compare the values in the common columns directly after ordering
comparison_result <- all.equal(df1_ordered[, c("Signature", "CTAB", "Genotype")], 
                               df2_ordered[, c("Signature", "CTAB", "Genotype")], 
                               check.attributes = FALSE)

# Check the result and print appropriate message
if (isTRUE(comparison_result)) {
  cat("The common columns 'Signature', 'CTAB', and 'Genotype' are identical between df1 and df2.\n")
} else {
  cat("The common columns 'Signature', 'CTAB', and 'Genotype' are not identical between df1 and df2. Differences:\n")
  print(comparison_result)
}

#### Creating a Correspondence Table and mapping the corresponding identifiers
# Define the correspondence table where the new column is named "TNC_id"
TNC_to_TNC_id <- data.frame(
  TNC = c("No_data", "Unchanged", "Underexpression", "Overexpression"),  # Original TNC values
  TNC_id = c(0, 1, 2, 3)  # The corresponding IDs for each TNC value
)

# Add the new "TNC_id" column to df1_ordered
df1_ordered <- df1_ordered %>%
  left_join(TNC_to_TNC_id, by = "TNC")  # "TNC_id" gets added based on matching "TNC" values

# Define the correspondence table where the new column is named "PFC_id"
PFC_to_PFC_id <- data.frame(
  PFC = c("TMB", "MSI", "Stemness"),  # Original PFC values
  PFC_id = c(1, 2, 3)  # The corresponding IDs for each PFC value
)

# Add the new "PFC_id" column to df1_ordered
df1_ordered <- df1_ordered %>%
  left_join(PFC_to_PFC_id, by = "PFC")  # "PFC_id" gets added based on matching "PFC" values

# Define the correspondence table where the new column is named "SCS_id"
SCS_to_SCS_id <- data.frame(
  SCS = c("negative", "positive"),  # Original SCS values
  SCS_id = c("N", "P")  # The corresponding IDs for each SCS value
)

# Add the new "SCS_id" column to df1_ordered
df1_ordered <- df1_ordered %>%
  left_join(SCS_to_SCS_id, by = "SCS")  # "SCS_id" gets added based on matching "SCS" values

# Exclude the "GSI" column from df1_ordered
df1_ordered <- df1_ordered %>%
  select(-GSI)

# Reorder columns in df1_ordered
df1_ordered <- df1_ordered %>%
  select(Signature, CTAB, Genotype, PFC, TNC, SCS, PFC_id, SCS_id, TNC_id, RCD)

# Reorder columns in df1_ordered
df1_ordered <- df1_ordered %>%
  select(Signature, CTAB, Genotype, PFC_id, SCS_id, TNC_id, RCD)

# Rename columns in df1_ordered
df1_ordered <- df1_ordered %>%
  rename(PFC = PFC_id, SCS = SCS_id, TNC = TNC_id)

# Order df1_ordered alphabetically by "Signature", "CTAB", and "Genotype"
df1_ordered <- df1_ordered %>%
  arrange(Signature, CTAB, Genotype)

# Order df1_ordered alphabetically by "Signature", "CTAB", and "Genotype"
df2_ordered <- df2_ordered %>%
  arrange(Signature, CTAB, Genotype)

# Even if the dimensions are different, join by the common columns
# Perform a full join to merge df1_ordered  and df2_ordered
merged_data_df1_df2_ordered <- full_join(
  df2_ordered,
  df1_ordered,
  by = c("Signature", "CTAB", "Genotype")
)

# Relocate columns 
merged_data_df1_df2_ordered <- merged_data_df1_df2_ordered %>%
  relocate(PFC, .after = GFC) %>%
  relocate(SCS, .after = PFC) %>%
  relocate(TNC, .after = SCS) %>%
  relocate(RCD, .after = TIC)

df2_ordered <- merged_data_df1_df2_ordered

# Count total number of NA values in df2_ordered
na_count <- sum(is.na(df2_ordered))

# Print the result
print(na_count)

###############################
###############################
###############################
###############################
# Create a dataframe with unique values of the "RCD" column, name it "unique_RCD", and sort it alphabetically
df_RCD <- as.data.frame(sort(unique(df1_ordered$RCD)))
colnames(df_RCD) <- "unique_RCD"

# Identify rows with empty or missing values in the "RCD" column
empty_rcd_rows <- which(is.na(df1_ordered$RCD))

df1_empty_rcd <- df1_ordered[df1_ordered$RCD == "", ]

# Print the rows that have missing values in "RCD"
if(length(empty_rcd_rows) > 0) {
  cat("Rows with missing values in the 'RCD' column:\n")
  print(empty_rcd_rows)
} else {
  cat("No missing values found in the 'RCD' column.\n")
}

# Select rows with empty values ("") in the "TNC" column
df1_empty_tnc <- df1_ordered[df1_ordered$TNC == "", ]

########################
########################
########################
# Count the number of values separated by "/" in each row and create a new column
df_RCD <- df_RCD %>%
  mutate(count_RCD_types = sapply(strsplit(as.character(unique_RCD), "/"), length))

# Estimate the global frequency of every regulated cell death type
df_RCD_frequencies <- df_RCD %>%
  # Split the RCD values by "/" and unnest them into individual rows
  mutate(RCD_split = strsplit(as.character(unique_RCD), "/")) %>%
  unnest(RCD_split) %>%
  group_by(RCD_split) %>%
  summarise(frequency = n()) %>%
  arrange(desc(frequency))

####################
# Ensure that df_RCD has columns "unique_RCD" and "count_RCD_types"
head(df_RCD)

# Perform the mapping and create the "RCD_type_number" column in df2_ordered
df2_ordered <- df2_ordered %>%
  # Use left_join to map based on the "RCD" and "unique_RCD" columns
  left_join(df_RCD, by = c("RCD" = "unique_RCD")) %>%
  # Relocate the "count_RCD_types" column after the "RCD" column and rename it to "RCD_type_number"
  relocate(count_RCD_types, .after = RCD) %>%
  rename(RCD_type_number = count_RCD_types)

# Rename columns in df2_ordered
df2_ordered <- df2_ordered %>%
  rename(RCD_types = RCD,  # Rename "RCD" to "RCD_types"
         RCD = RCD_type_number) %>%  # Rename "RCD_type_number" to "RCD"
  relocate(RCD, .before = RCD_types)  # Relocate the newly named "RCD" after "RCD_types"

# Remove leading and trailing spaces globally across all columns in df2_ordered
df2_ordered <- df2_ordered %>%
  mutate_all(trimws)

# Signature components order: CTAB-(GSI).GFC.PFC.SCS.TNC.HRC.SMC.TMC.TIC.RCD 
first_11_columns <- colnames(df2_ordered[1:11])
first_11_columns

# Remove subsequent duplicate ROWS ACROSS ALL VARIABLES, keeping only the first occurrence
df2_ordered_unique <- df2_ordered[!duplicated(df2_ordered), ]

df2_GSI <-df2_ordered_unique
dim(df2_GSI)
#
# Create the GSI (Genomic Series Identifier) DISTRIBUTION variable based on CTAB and Signature
df2_GSI <- df2_GSI %>%
  group_by(CTAB) %>%  # Group by CTAB
  mutate(GSI = row_number()) %>%  # Create the series number within each CTAB group
  ungroup() %>%  # Ungroup to finish the transformation
  relocate(GSI, .after = CTAB)  # Move GSI column to be after the CTAB column

# Distribution of signatures per cancer type
# Summarize the DISTRIBUTION number of Signatures per CTAB and add total row
summary_signatures_GSI <- df2_GSI %>%
  group_by(CTAB) %>%
  summarise(num_signatures = n()) %>%  # Count the number of rows (Signatures) for each CTAB
  ungroup() %>%
  # Add a row for the total number of CTABs and Signatures
  bind_rows(
    summarise(., CTAB = "Total", num_signatures = sum(num_signatures))
  )

rio::export(df2_GSI, "df2_ordered_unique.tsv")
rio::export(summary_signatures_GSI, "Summary_distribution_of_signatures_per_CTAB.tsv")

# Create the bar plot for "num_signatures" per "CTAB"
# Exclude the row where CTAB is "Total"
filtered_data <- summary_signatures_GSI %>%
  filter(CTAB != "Total")

# Create the bar plot for "num_signatures" per "CTAB", excluding "Total"
ggplot(filtered_data, aes(x = CTAB, y = num_signatures, fill = CTAB)) +
  geom_bar(stat = "identity") +
  labs(title = "Number of Signatures per Cancer Type (Excluding 'Total')",
       x = "Cancer Type",
       y = "Number of Signatures") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for readability
  scale_fill_discrete(name = "Cancer Type")

dev.off()

df2_ordered_old <- df2_ordered

df2_ordered <-  df2_GSI

# Apply trimws() to remove leading and trailing spaces from all columns in 'resultado'
df2_ordered <- df2_ordered %>%
  mutate(across(everything(), ~trimws(.)))  # Applies trimws to all columns

########################
########################
########################
# Create the "Nomenclature" column and place it after the "Signature" column
df2_ordered <- df2_ordered %>%
  mutate(
    Nomenclature = apply(df2_ordered[, 2:13], 1, function(row) {
      # Concatenate the first and second columns with an "_" and the rest with "."
      paste(c(paste(row[1], row[2], sep = "-"), row[3:11]), collapse = ".")
    })
  ) %>%
  relocate(Nomenclature, .after = Signature)  # Place the new column after "Signature"

unique_signatures <-  as.data.frame(unique(df2_ordered$Nomenclature))

# Rename the column "unique(df2_ordered$Nomenclature)" back to "Nomenclature" in unique_signatures
unique_signatures <- unique_signatures %>%
  rename(Nomenclature = `unique(df2_ordered$Nomenclature)`)

# Count the number of unique values in the "Signature" column
number_unique_nomenclatures <- length(unique(unique_signatures$Nomenclature))

rio::export(unique_signatures, "Nomenclature.tsv")

# View the updated dataframe
print(colnames(unique_signatures))

final_signatures <-  df2_ordered
rio::export(final_signatures, "Final_Signatures.tsv")

# Count the number of unique values in the "Signature" column
number_unique_signatures <- length(unique(final_signatures$Signature))

# Check if there are any duplicate rows across all columns in final_signatures
has_duplicates <- any(duplicated(final_signatures))

# Print the result
if (has_duplicates) {
  print("There are duplicate rows in final_signatures.")
} else {
  print("There are no duplicate rows in final_signatures.")
}

# Order the dataframe df_RCD by the "count_RCD_types" column in ascending order
df_RCD <- df_RCD %>%
  arrange(count_RCD_types)

rio::export(df_RCD, "df_unique_RCD_types.tsv")

# None_effect value_signatures: Select all rows where the value in "HRC_series" is "1A2A3A4A"
filtered_rows_HRC <- final_signatures %>%
  filter(HRC_series == "1A2A3A4A")

filtered_rows_SMC <- final_signatures %>%
  filter(SMC_series == "1A2A3A4A")

#### No effect value signatures:
filtered_rows_HRC_SMC <- final_signatures %>%
  filter(HRC_series == "1A2A3A4A" | SMC_series == "1A2A3A4A")

#### Select all rows that are meaningful by "HRC_series" and "SMC_series"
# Select all rows where the value in "HRC_series" an  ""SMC_series) is different that "1A2A3A4A" 
filtered_rows_meaningful <- final_signatures %>%
  filter(HRC_series != "1A2A3A4A", SMC_series != "1A2A3A4A")

# Select all rows that are meaningful by the value in "HRC_series", "SMC_series", "immune_classification", and "microenvironment_classification") is different that "1A2A3A4A"
meaningful_filtered_signatures <- final_signatures %>%
  filter(HRC_series != "1A2A3A4A", SMC_series != "1A2A3A4A", immune_classification != "NS", microenvironment_classification != "NS")

rio::export(meaningful_filtered_signatures, "Meaningful_filtered_signatures.tsv") 

rio::export(final_signatures, "Final_Signatures.tsv") 

final_signatures <- import("Final_Signatures.tsv")

########################
########################
########################
# Get the unique values from HRC_series and SMC_series
Unique_HRC_array <- unique(final_signatures$HRC_series)
Unique_SMC_array <- unique(final_signatures$SMC_series)

# Convert unique values to data frames
df_HRC <- data.frame(series = Unique_HRC_array)  # Rename the column as "series"

df_SMC <- data.frame(series = Unique_SMC_array)  # Rename the column as "series"

# Order the dataframe alphabetically by the "series" column
unique_df_HRC_filtered <- df_HRC %>%
  arrange(series)

unique_df_SMC_filtered <- df_SMC %>%
  arrange(series)

columns <- as.data.frame(colnames(final_signatures))

##########
########## Making the ranking templates for HRC and SMC components
##########
# Creating the mapping template (df12) as a named vector for easy lookup
df12 <- unique_combined_HRC_SMC_array_filtered

HRC_rank_template <- c('A' = 0, 'B' = 1, 'C' = 1, '1' = 2, '2' = 2, '3' = 2, '4' = 1)

# Define a function to rank each series
HRC_rank_series <- function(series_value) {
  # Split the series into individual parts, e.g., "1A2A3D4A" -> c("1A", "2A", "3D", "4A")
  split_series <- strsplit(series_value, "(?<=[A-D])", perl = TRUE)[[1]]
  
  total_rank <- 0
  
  # Loop through each part of the series (e.g., "1A", "2A", etc.)
  for (part in split_series) {
    # Extract the number and letter
    num <- substr(part, 1, 1)  # First character is the number
    letter <- substr(part, 2, 2)  # Second character is the letter
    
    # Check if letter is "A"; if so, skip adding the rank for this number-letter pair
    if (letter == "A") {
      next  # Skip to the next part
    }
    
    # Add the ranks from the template
    total_rank <- total_rank + HRC_rank_template[num] + HRC_rank_template[letter]
  }
  
  return(total_rank)
}

# Apply the ranking function to the "series" column in df12
df12 <- df12 %>%
  mutate(HRC_Rank = sapply(series, HRC_rank_series))

# Create the mapping template (df13) as a named vector for easy lookup
SMC_rank_template <- c('A' = 0, 'B' = 1, 'C' = 1,  'D' = 1, '1' = 2, '2' = 2, '3' = 2, '4' = 1)

# Define a function to rank each series
SMC_rank_series <- function(series_value) {
  # Split the series into individual parts, e.g., "1A2A3D4A" -> c("1A", "2A", "3D", "4A")
  split_series <- strsplit(series_value, "(?<=[A-D])", perl = TRUE)[[1]]
  
  total_rank <- 0
  
  # Loop through each part of the series (e.g., "1A", "2A", etc.)
  for (part in split_series) {
    # Extract the number and letter
    num <- substr(part, 1, 1)  # First character is the number
    letter <- substr(part, 2, 2)  # Second character is the letter
    
    # Check if letter is "A"; if so, skip adding the rank for this number-letter pair
    if (letter == "A") {
      next  # Skip to the next part
    }
    
    # Add the ranks from the template
    total_rank <- total_rank + SMC_rank_template[num] + SMC_rank_template[letter]
  }
  
  return(total_rank)
}

# Apply the ranking function to the "series" column in df12
df12 <- df12 %>%
  mutate(SMC_Rank = sapply(series, SMC_rank_series))

unique_combined_HRC_SMC_array_filtered_rank_template <- df12

df12 <- unique_combined_HRC_SMC_array_filtered_rank_template

# Renaming variables: "array_series_number", "series", "HRC_Rank", "SMC_Rank" 
# df12 <-  rename(df12,HRC_SMC_rank = SMC_Rank)
df13 <- rename(df12, HRC_SMC = array_series_number)
df12 <-  df13

####
####
# Duplicate HRC_SMC and rename variables
df12 <- df12 %>%
  mutate(SMC = HRC_SMC) %>%  # Duplicate HRC_SMC and name it SMC
  rename(
    HRC = HRC_SMC,           # Rename HRC_SMC to HRC
    `Array series class` = series  # Rename series to Array series class
  ) %>%
  select(HRC, SMC, `Array series class`, HRC_Rank, SMC_Rank)  # Reorder variables

# View the updated dataframe
head(df12)
####

# Load the dplyr package
library(dplyr)

# Replace empty strings and NA values with "N/A" in all columns
df12 <- df12 %>%
  mutate(across(everything(), ~ ifelse(. == "" | is.na(.), "N/A", as.character(.))))

HRC_SMC_map <- df12

rio::export(df12, "unique_combined_HRC_SMC_array_filtered_rank_template.xlsx")

rio::export(df12, "unique_combined_HRC_SMC_array_filtered_rank_template.tsv")
rio::export(columns, "columns.tsv")

df12 <-  import("unique_combined_HRC_SMC_array_filtered_rank_template.tsv")

#######
#######
#######
### Signature components rank maps!
#######
#######
# Genomic Feature Ranking
GFC_map <- data.frame(
  GFC = c("1", "2", "3", "4", "5", "6", "7"),
  GFC_Rank = c(7, 6, 5, 4, 3, 2, 1) # Protein Expression = 7,
  # Mutations = 6, Copy Number Variation = 5, miRNA Expression = 4, 
  # Transcript Expression = 3, mRNA Expression = 2, CpG Methylation =1
)

# Phenotypic Feature Ranking
PFC_map <- data.frame(
  PFC = c("1", "2", "3"),
  PFC_Rank = c(3, 2, 1) # TMB=3, MSI=2, TSM=1
)

# Spearman Correlation Sign (SCS) Feature Ranking
# Ranking SCS features based on their correlation signs.
SCS_map <- data.frame(
  SCS = c("P", "N"),
  SCS_Rank = c(2, 1) # P=Positive=2 , N=Negative=1
)

# TCGA vs. GTEx expression Code (TNC) Feature Ranking
# Ranking TNC features based on their expression codes.
# NOTE: MODIFIED revised on 01/12/2024
TNC_map <- data.frame(
  TNC = c("0", "1", "2", "3"),
  TNC_Rank = c(0, 1, 2, 2) # No_data; unchanged=1, underexpressed=3, overexpressed=3 NOTE: MODIFIED revised on 01/12/2024
)

# TMC and TIC Feature Ranking
# A combined ranking for TMC and TIC features.
TMC_TIC_map <- data.frame(
  TMC = c("1", "2", "3", "4"),
  TIC = c("1", "2", "3", "4"),
  TMC_Rank = c(7, 4, 1, 0), 
  TIC_Rank = c(7, 4, 1, 0) # "Hot", "Variable", "Cold", "Anti-tumoral, "Dual", "Pro-tumoral", "NS"
)

# RCD Feature Ranking
RCD_map <- data.frame(
  RCD = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "12"),
  RCD_Rank = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 12) # 1 for signatures whose gene members are associated with one RCD type and so on and so for
) 

df14 <-import("Final_signatures.tsv") 

df20 <-df12

HRC_SMC_map <- df20

rio::export(HRC_SMC_map, "HRC_SMC_map.tsv")

### Extracting Signature nomenclature codes for ranking calculation
### # Extract the "Nomenclature" column from df14 and create a new data frame
df14 <- final_signatures

df14 <- df14 %>%
  mutate_all(trimws)

# Get the rows that contain any NA values
rows_with_NA <- df14[!complete.cases(df14), ]

# Print the rows with NA values
print(rows_with_NA)
Nomenclature <- df14 %>%
  select(Nomenclature)

Nomenclature <- Nomenclature %>%
  mutate_all(trimws)

df15 <- Nomenclature

# Step 1: Duplicate the "Nomenclature" column and rename it to "Signature_nomenclature"
df15$Signature_nomenclature <- df15$Nomenclature

# Step 2: Relocate "Nomenclature" to come after "Signature_nomenclature"
df15 <- df15[, c(setdiff(names(df15), "Nomenclature"), "Nomenclature")]

# Split the "Nomenclature" column using "-" and "." as delimiters
df15_split <- df15 %>%
  separate(Nomenclature, into = c("CTAB", "GSI", "GFC", "PFC", "SCS", "TNC", "HRC", "SMC", "TMC", "TIC", "RCD"), sep = "[-\\.]")

df15_split <- df15_split %>%
  mutate_all(trimws)

#####
#####
###################### 
###################### 
###################### Crafted logic for ranking the signatures 
###################### 
###################### 
######
###### 

# Add mappings for all other feature maps as needed (already provided)

# Function to map ranks and add new ranking columns at the right end of the dataframe
assign_ranks <- function(df, col_name, map_df, map_col, rank_col) {
  # Match the values in the target column with the mapping column and extract the rank
  rank_values <- map_df[[rank_col]][match(df[[col_name]], map_df[[map_col]])]
  # Create the new column name
  rank_col_name <- paste0(col_name, "_ranking")
  # Add the new ranking column at the right end of the dataframe
  df[[rank_col_name]] <- rank_values
  return(df)
}

# Apply the function for each column with its respective map

# GFC ranking
df15_split <- assign_ranks(df15_split, "GFC", GFC_map, "GFC", "GFC_Rank")

# PFC ranking
df15_split <- assign_ranks(df15_split, "PFC", PFC_map, "PFC", "PFC_Rank")

# SCS ranking
df15_split <- assign_ranks(df15_split, "SCS", SCS_map, "SCS", "SCS_Rank")

# TNC ranking
df15_split <- assign_ranks(df15_split, "TNC", TNC_map, "TNC", "TNC_Rank")

# HRC ranking
df15_split <- assign_ranks(df15_split, "HRC", HRC_SMC_map, "HRC", "HRC_Rank")

# SMC ranking
df15_split <- assign_ranks(df15_split, "SMC", HRC_SMC_map, "SMC", "SMC_Rank")

# TMC ranking
df15_split <- assign_ranks(df15_split, "TMC", TMC_TIC_map, "TMC", "TMC_Rank")

# TIC ranking
df15_split <- assign_ranks(df15_split, "TIC", TMC_TIC_map, "TIC", "TIC_Rank")

# RCD ranking
df15_split <- assign_ranks(df15_split, "RCD", RCD_map, "RCD", "RCD_Rank")

#### mutate as.numeruc the relevant column variables
# Specify the columns to check
columns_to_check <- c(3:5, 7:21)

# Check for "NA" or "NS" values in the specified columns
na_ns_check <- df15_split %>%
  select(all_of(columns_to_check)) %>%
  summarise_all(~ any(. %in% c("NA", "NS", "")))

# Print the results
print(na_ns_check)

# Convert the specified columns to numeric
df15_split <- df15_split %>%
  mutate(across(c(3:5, 7:21), as.numeric))

# Preview the updated dataframe
str(head(df15_split))

#######
#######
#######
# Ensure the Final_Rank column is numeric
df15_split$Final_Rank <- rowSums(
  sapply(df15_split[, c("GFC_ranking", "PFC_ranking", "SCS_ranking", 
                        "TNC_ranking", "TMC_ranking", "TIC_ranking", 
                        "HRC_ranking", "SMC_ranking", "RCD_ranking")], as.numeric),
  na.rm = TRUE
)

# Preview the updated dataframe
print(head(df15_split$Final_Rank))

### Handling "NA" values in the HRC_ranking component because of the absence of CNV variable "D"
# Check for NA values in the "Signature_nomenclature" column
NA_problem <- is.na(df15_split$Signature_nomenclature)

# Check for NA values in the entire dataframe "df15_split"
NA_problem <- is.na(df15_split)

# Optionally, to count the total number of NA values in the entire dataframe:
total_NA_count <- sum(NA_problem)

# Display the number of NA values
total_NA_count

# Get the rows that contain any NA values
rows_with_NA <- df15_split[!complete.cases(df15_split), ] # rows in which the HRC_ranking is not possible because the classification is rhe series array has at least one  "D" (CNV)

# Display the number of TRUE values (i.e., how many NA values are present)
sum(NA_problem)

final_signature_ranking <- df15_split

columnFR <- as.data.frame(colnames(final_signature_ranking))

# Display the updated dataframe with the new ranking columns added to the right
rio::export(final_signature_ranking, "Final_signature_ranking.tsv")

##########
##########
########## Adding the final raking columns to the global signature df##########
##########
##########
# Ensure both columns are trimmed of white space and are in the same case
df70 <- final_signatures
df71 <- final_signature_ranking

# Apply trimws() to remove leading and trailing spaces from all columns 
df70 <- df70 %>%
  mutate(across(everything(), ~trimws(.)))  # Applies trimws to all columns

df71 <- df71 %>%
  mutate(across(everything(), ~trimws(.)))  # Applies trimws to all columns

# Step 1: Order df70 and df71 alphabetically
df70 <- df70[order(df70$Nomenclature), ]
df71 <- df71[order(df71$Signature_nomenclature), ]

# # Create the 'Ranking' column in 'df70' by matching with 'df71'
# df70$Ranking <- df71$Final_Rank[match(df70$Nomenclature, df71$Signature_nomenclature)]

# Create the 'Ranking' column in df70 using as.numeric to ensure numeric conversion
df70$Ranking <- as.numeric(df71$Final_Rank[match(df70$Nomenclature, df71$Signature_nomenclature)])

# Preview the updated dataframe
print(head(df70$Ranking))

# Reorder the columns to place 'Ranking' before the 'Signature' column
df70 <- df70[, c('Ranking', setdiff(names(df70), 'Ranking'))]

final_signatures_full_ranking <- df70

rio::export(final_signatures_full_ranking, "final_signatures_full_ranking.tsv")


df16 <-  final_signature_ranking

# Convert all columns in the "df16" dataframe to numeric, except "Signature_nomenclature", "CTAB", and "SCS"
df16[] <- lapply(df16, function(x) {
  if (any(names(df16) %in% c("Signature_nomenclature", "CTAB", "SCS"))) {
    x  # Retain original values for specified columns
  } else {
    as.numeric(as.character(x))  # Convert to numeric for other columns
  }
})

# Estimate the mean of the values under "Final_Rank" in the "df16" dataframe
mean_Final_Rank <- mean(df16$Final_Rank, na.rm = TRUE)

# Estimate the median of the values 
median_Final_Rank <- median(df16$Final_Rank, na.rm = TRUE)

# Estimate the quartiles (1st, 2nd/median, and 3rd) of the values 
quartiles_Final_Rank <- quantile(df16$Final_Rank, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)

print(quartiles_Final_Rank)

# Select rows where HRC_Rank and SMC_Rank are 11, and TMC_Rank and TIC_Rank are 7
top_ranked_signatures <- df16[df16$SMC_ranking == 11 &
                              df16$TMC_ranking == 7 &
                              df16$TIC_ranking == 7, ]

top_ranked_signatures_v2 <- df16[
  df16$SMC_ranking == 11 & 
    df16$TMC_ranking == 7 & 
    df16$TIC_ranking == 7 & 
    !(df16$TNC_ranking %in% c(0, 1)), 
]

# Select rows where HRC_Rank and SMC_Rank are 11, and TMC_Rank and TIC_Rank are 7
worst_ranked_signatures <- df16[df16$SMC_ranking == 11 & 
                                df16$TMC_ranking == 1 & 
                                df16$TIC_ranking == 1, ]

# Select rows where HRC_Rank and SMC_Rank are 11, and TMC_Rank and TIC_Rank are 7
worst_ranked_signatures_v2 <- df16[df16$SMC_ranking == 11 & 
                                  df16$TMC_ranking == 1 & 
                                  df16$TIC_ranking == 1 & 
                                  !(df16$TNC_ranking %in% c(0, 1)),
]

rio::export(top_ranked_signatures, "top_ranked_signatures.tsv")

rio::export(worst_ranked_signatures, "worst_ranked_signatures.tsv")

rio::export(top_ranked_signatures_v2, "top_ranked_signatures_v2.tsv")

rio::export(worst_ranked_signatures_v2, "worst_ranked_signatures_v2.tsv")

# Select rows from df14 where "Nomenclature" matches the "Signature_nomenclature" in "top_ranked_signatures"
complete_top_signatures <- df14 %>%
  filter(Nomenclature %in% top_ranked_signatures$Signature_nomenclature)

complete_top_signatures_v2 <- df14 %>%
  filter(Nomenclature %in% top_ranked_signatures_v2$Signature_nomenclature)

df141 <- complete_top_signatures

df141_mRNA <- df141 %>%
  filter(Genotype == "mRNA")

# Select rows from df14 where "Nomenclature" matches the "Signature_nomenclature" in "worst_ranked_signatures"
complete_worst_signatures <- df14 %>%
  filter(Nomenclature %in% worst_ranked_signatures$Signature_nomenclature)

# Select rows from df14 where "Nomenclature" matches the "Signature_nomenclature" in "worst_ranked_signatures"
complete_worst_signatures_v2 <- df14 %>%
  filter(Nomenclature %in% worst_ranked_signatures_v2$Signature_nomenclature)

rio::export(complete_top_signatures, "complete_top_signatures.tsv")
rio::export(complete_top_signatures_v2, "complete_top_signatures_v2.tsv")

rio::export(complete_worst_signatures, "complete_worst_signatures.tsv")
rio::export(complete_worst_signatures_v2, "complete_worst_signatures_v2.tsv")

# Calculate the mean, median, and quartiles of the "Final_Rank" column in df16
mean_rank <- mean(df16$Final_Rank, na.rm = TRUE)
median_rank <- median(df16$Final_Rank, na.rm = TRUE)
quartiles_rank <- quantile(df16$Final_Rank, probs = c(0.25, 0.50, 0.75), na.rm = TRUE)

# Create a histogram for the "Final_Rank" column in df16 with mean, median, and quartile lines
#distribution_plot <-# to save data
ggplot(df16, aes(x = Final_Rank)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") + 
  geom_vline(aes(xintercept = mean_rank), color = "red", linetype = "dashed", size = 1, show.legend = TRUE) + 
  geom_vline(aes(xintercept = median_rank), color = "blue", linetype = "solid", size = 1, show.legend = TRUE) +
  geom_vline(aes(xintercept = quartiles_rank[1]), color = "green", linetype = "dotted", size = 1) + 
  geom_vline(aes(xintercept = quartiles_rank[2]), color = "green", linetype = "dotted", size = 1) + 
  labs(title = "Histogram of Final Rank Distribution with Mean, Median, and Quartiles",
       x = "Ranking distribution",
       y = "Signatures per rank number") +
  theme_minimal() +
  annotate("text", x = mean_rank, y = Inf, label = "Mean", color = "red", vjust = 2, hjust = 1) +
  annotate("text", x = median_rank, y = Inf, label = "Median", color = "blue", vjust = 2, hjust = 1) +
  annotate("text", x = quartiles_rank[1], y = Inf, label = "Q1", color = "green", vjust = 2, hjust = 1) +
  annotate("text", x = quartiles_rank[2], y = Inf, label = "Q3", color = "green", vjust = 2, hjust = 1)

dev.off()

# Perform the Kolmogorov-Smirnov test for normality
ks_test <- ks.test(df16$Final_Rank, "pnorm", mean(df16$Final_Rank, na.rm = TRUE), sd(df16$Final_Rank, na.rm = TRUE))

# Print the results of the Kolmogorov-Smirnov test
print(ks_test)

# Interpretation:
# - If the p-value from the Kolmogorov-Smirnov test is greater than 0.05, the data is considered to be normally distributed.

# Create a Q-Q plot to visually assess normality
qqnorm(df16$Final_Rank)
qqline(df16$Final_Rank, col = "red", lwd = 2)

dev.off()

# Count the frequency of each unique value in Final_Rank
table(df16$Final_Rank)

# Check the number of unique values
length(unique(df16$Final_Rank))

# Add small random noise to break ties
# Jittering (Adding Noise to Break Ties): 
# to break ties while maintaining the overall distribution structure:
 
df16$Jittered_Final_Rank <- jitter(df16$Final_Rank, factor = 0.1)

# Plot the jittered data to check the distribution
qqnorm(df16$Jittered_Final_Rank)
qqline(df16$Jittered_Final_Rank, col = "red", lwd = 2)

dev.off()

# Group Final_Rank into quantiles (e.g., deciles)
df16$Rank_Quantiles <- cut(df16$Final_Rank, breaks = quantile(df16$Final_Rank, probs = seq(0, 1, 0.1), na.rm = TRUE), include.lowest = TRUE)

# Check the distribution of ranks within each quantile
quartile_ranking <- as.data.frame(table(df16$Rank_Quantiles)) 
quartile_ranking

# Add small random noise (jitter) to break ties in Final_Rank
df16$Jittered_Final_Rank <- jitter(df16$Final_Rank, factor = 0.1)

# Plot the jittered data to check the distribution
qqnorm(df16$Jittered_Final_Rank)
qqline(df16$Jittered_Final_Rank, col = "red", lwd = 2)

dev.off()

# Perform the Kolmogorov-Smirnov test on the jittered data
ks_test_jittered <- ks.test(df16$Jittered_Final_Rank, "pnorm", 
                            mean(df16$Jittered_Final_Rank, na.rm = TRUE), 
                            sd(df16$Jittered_Final_Rank, na.rm = TRUE))

# Print the results of the KS test on jittered data
print(ks_test_jittered)

# Group Final_Rank into deciles (10 groups)
df16$Rank_Quantiles <- cut(df16$Final_Rank, breaks = quantile(df16$Final_Rank, 
                                                              probs = seq(0, 1, 0.1), na.rm = TRUE), 
                           include.lowest = TRUE)

# Display the count of values in each quantile group
table(df16$Rank_Quantiles)

# Visualize the quantiles with a box plot
ggplot(df16, aes(x = Rank_Quantiles, y = Final_Rank)) +
  geom_boxplot() +
  labs(title = "Final Rank by Quantile",
       x = "Quantile Group",
       y = "Final Rank") +
  theme_minimal()

dev.off()

##### Conclusion although the p value is quite small (non normal distribution)
#####  the plots show that most of the data are normally distributed for most purposes)

### Wilcoxon Test (Non-parametric alternative to the t-test)
### # Example: Split Final_Rank into two groups based on a median split
median_split <- median(df16$Final_Rank, na.rm = TRUE)
df16$Group <- ifelse(df16$Final_Rank <= median_split, "Low", "High")

# Perform the Wilcoxon test to compare the two groups
wilcox_test <- wilcox.test(Final_Rank ~ Group, data = df16)

# Print the results
print(wilcox_test)

# Perform the Wilcoxon test again to obtain W and p-value
wilcox_test <- wilcox.test(Final_Rank ~ Group, data = df16)
w_value <- wilcox_test$statistic  # Extract the W statistic
p_value <- wilcox_test$p.value    # Extract the p-value

# Create a box plot and include the Wilcoxon test results (W and p-value)
ggplot(df16, aes(x = Group, y = Final_Rank, fill = Group)) +
  geom_boxplot(outlier.color = "red", outlier.shape = 16, outlier.size = 2) +
  labs(title = "Box Plot of Final Rank by Group with Wilcoxon Test Results",
       x = "Group (Low vs High)",
       y = "Final Rank") +
  theme_minimal() +
  scale_fill_manual(values = c("Low" = "skyblue", "High" = "orange")) +
  annotate("text", x = 1.5, y = max(df16$Final_Rank), label = paste0("W = ", round(w_value, 2)), size = 5, color = "black") +
  annotate("text", x = 1.5, y = max(df16$Final_Rank) - 1, label = paste0("p-value = ", format(p_value, scientific = TRUE)), size = 5, color = "black")

dev.off()

#### Bootstrap Method
# Define a function to calculate the mean for bootstrapping
boot_mean <- function(data, indices) {
  return(mean(data[indices], na.rm = TRUE))
}

# Perform bootstrapping on the Final_Rank column
set.seed(123)  # for reproducibility
bootstrap_results <- boot(data = df16$Final_Rank, statistic = boot_mean, R = 1000)

# Print the bootstrap estimate and confidence intervals
print(bootstrap_results)

# Calculate the 95% confidence intervals
boot_ci <- boot.ci(bootstrap_results, type = "perc")
print(boot_ci)

# # Plot the distribution of the bootstrapped means
# Extract the relevant bootstrap statistics
boot_mean <- bootstrap_results$t0  # The original mean
boot_bias <- mean(bootstrap_results$t) - boot_mean  # Bias from the bootstrap
boot_se <- sd(bootstrap_results$t)  # Standard error from the bootstrap

# Plot the bootstrapped distribution of means and include the statistics
boot_means <- bootstrap_results$t  # Bootstrapped means

ggplot() +
  geom_histogram(aes(x = boot_means), bins = 30, fill = "skyblue", color = "black") +
  labs(title = "Bootstrapped Distribution of Means with Statistics",
       x = "Mean of Final_Rank",
       y = "Frequency") +
  theme_minimal() +
  # Annotate with bootstrap statistics, shifted to the upper left of the plot
  annotate("text", x = min(boot_means), y = Inf, label = paste0("Mean = ", round(boot_mean, 4)), 
           vjust = 2, hjust = -0.1, size = 5, color = "black") +
  annotate("text", x = min(boot_means), y = Inf, label = paste0("Bias = ", round(boot_bias, 6)), 
           vjust = 4, hjust = -0.1, size = 5, color = "black") +
  annotate("text", x = min(boot_means), y = Inf, label = paste0("SE = ", round(boot_se, 6)), 
           vjust = 6, hjust = -0.1, size = 5, color = "black")

dev.off()

#######################
#######################
##### Analysis of top signatures ####
#######################
#######################
# setwd("D:/Pré-artigo 5-optosis model/Dataframes all metrics")

data <- import("Drafts_2/gene_info_trancritos_excel2.xlsx")

data <- data %>%
  mutate_all(trimws)

miRNA_database <-  import("Drafts_2/miRNA_database.xlsx")

miRNA_database <- miRNA_database %>%
  mutate_all(trimws)

df24 <- top_ranked_signatures

df24 <- df24 %>%
  mutate_all(trimws)

df14 <- df14 %>%
  mutate_all(trimws)

# Step 0: Filter df14 based on matching values in the "Nomenclature" column of df14 and "Signature_nomenclature" column in df24
df14_filtered <- df14 %>%
  filter(Nomenclature %in% df24$Signature_nomenclature)

# Step 1: Copy the filtered columns "Signature" and "Nomenclature" from df14_filtered to df25
df25 <- df14_filtered %>%
  select(Signature, Nomenclature)

# Step 2: Create a new column "Gene_conversion" after "Signature"
df25 <- df25 %>%
  mutate(Gene_conversion = NA)  # Initialize the column with NA values

# Step 3: Populate "Gene_conversion" using left_join for better matching reliability

# For values starting with "ENST", join with the "data" dataframe on "Transcript_ID"
df25 <- df25 %>%
  left_join(data %>% select(Transcript_ID, Display_Name), by = c("Signature" = "Transcript_ID")) %>%
  mutate(Gene_conversion = ifelse(!is.na(Display_Name), Display_Name, Gene_conversion)) %>%
  select(-Display_Name)  # Remove the temporary column

# For values starting with "hsa-miR", join with the "miRNA_database" on "Mature1" and "Mature2"
df25 <- df25 %>%
  left_join(miRNA_database %>% select(Mature1, Gene_Symbol), by = c("Signature" = "Mature1")) %>%
  mutate(Gene_conversion = ifelse(!is.na(Gene_Symbol), Gene_Symbol, Gene_conversion)) %>%
  select(-Gene_Symbol)  # Remove the temporary column

df25 <- df25 %>%
  left_join(miRNA_database %>% select(Mature2, Gene_Symbol), by = c("Signature" = "Mature2")) %>%
  mutate(Gene_conversion = ifelse(!is.na(Gene_Symbol), Gene_Symbol, Gene_conversion)) %>%
  select(-Gene_Symbol)  # Remove the temporary column

# Step 1: Extract gene symbols, transcript names, and miRNA from the "Signature" column
df30 <- df25 %>%
  # Use str_replace_all to remove parentheses and separate by "+" into individual elements
  mutate(Signature_cleaned = str_replace_all(Signature, "[()]", "")) %>%  # Remove parentheses
  separate_rows(Signature_cleaned, sep = " \\+ ") %>%  # Split by " + " into separate rows
  distinct() %>%  # Remove duplicate rows
  filter(str_detect(Signature_cleaned, "^ENST|^hsa-miR|^[A-Za-z]")) %>%  # Keep only gene symbols, ENST transcripts, or miRNA
  arrange(Signature_cleaned)  # Order alphabetically

# Step 2: Create the df30 with the unique and ordered entries
df30 <- df30 %>%
  select(Signature_cleaned) %>%
  rename(Extracted_Signature = Signature_cleaned)

# Step 1: Create a new column "Gene_Symbol" in df30 and populate it based on the mapping logic
df30 <- df30 %>%
  rowwise() %>%
  mutate(Gene_Symbol = case_when(
    
    # For values with prefix "ENST", look up in the "data" dataframe using "Transcript_ID"
    str_detect(Extracted_Signature, "^ENST") ~ 
      ifelse(Extracted_Signature %in% data$Transcript_ID, 
             data$Display_Name[which(data$Transcript_ID == Extracted_Signature)], 
             NA_character_),
    
    # For values with prefix "hsa-miR", look up in "miRNA_database" for Mature1
    str_detect(Extracted_Signature, "^hsa-miR") & Extracted_Signature %in% miRNA_database$Mature1 ~ 
      ifelse(Extracted_Signature %in% miRNA_database$Mature1, 
             miRNA_database$Gene_Symbol[which(miRNA_database$Mature1 == Extracted_Signature)], 
             NA_character_),
    
    # For values with prefix "hsa-miR", look up in "miRNA_database" for Mature2
    str_detect(Extracted_Signature, "^hsa-miR") & Extracted_Signature %in% miRNA_database$Mature2 ~ 
      ifelse(Extracted_Signature %in% miRNA_database$Mature2, 
             miRNA_database$Gene_Symbol[which(miRNA_database$Mature2 == Extracted_Signature)], 
             NA_character_),
    
    # If no match is found, leave as NA
    TRUE ~ NA_character_
  ))

# Fill NA values in "Gene_Symbol" with corresponding values from "Extracted_Signature"
df30 <- df30 %>%
  mutate(Gene_Symbol = ifelse(is.na(Gene_Symbol), Extracted_Signature, Gene_Symbol))

# Remove duplicates based on "Gene_Symbol"
df30 <- df30 %>%
  distinct(Gene_Symbol, .keep_all = TRUE)

# Order the dataframe by "Gene_Symbol" alphabetically
df30 <- df30 %>%
  arrange(Gene_Symbol)

# Exclude the "Extracted_Signature" column
df30 <- df30 %>%
  select(-Extracted_Signature)

Gene_Symbols_top_signatures <- df30

# Order the dataframe "Gene_Symbols_top_signatures" alphabetically by the "Gene_Symbol" column
Gene_Symbols_top_signatures <- Gene_Symbols_top_signatures[order(Gene_Symbols_top_signatures$Gene_Symbol), ]

# Remove leading and trailing spaces globally across all columns
Gene_Symbols_top_signatures <- Gene_Symbols_top_signatures %>%
  mutate_all(trimws)

rio::export(Gene_Symbols_top_signatures, "Gene_Symbols_top_signatures.tsv")

####################
####################
#### Analysis of worst signatures #####
####################
####################
####################
df44 <- worst_ranked_signatures

df44 <- df44 %>%
  mutate_all(trimws)

# Step 0: Filter df14 based on matching values in the "Nomenclature" column of df14 and "Signature_nomenclature" column in df44
df14_filtered <- df14 %>%
  filter(Nomenclature %in% df44$Signature_nomenclature)

# Step 1: Copy the filtered columns "Signature" and "Nomenclature" from df14_filtered to df45
df45 <- df14_filtered %>%
  select(Signature, Nomenclature)

# Step 2: Create a new column "Gene_conversion" after "Signature"
df45 <- df45 %>%
  mutate(Gene_conversion = NA)  # Initialize the column with NA values

# Step 3: Populate "Gene_conversion" using left_join for better matching reliability

# For values starting with "ENST", join with the "data" dataframe on "Transcript_ID"
df45 <- df45 %>%
  left_join(data %>% select(Transcript_ID, Display_Name), by = c("Signature" = "Transcript_ID")) %>%
  mutate(Gene_conversion = ifelse(!is.na(Display_Name), Display_Name, Gene_conversion)) %>%
  select(-Display_Name)  # Remove the temporary column

# For values starting with "hsa-miR", join with the "miRNA_database" on "Mature1" and "Mature2"
df45 <- df45 %>%
  left_join(miRNA_database %>% select(Mature1, Gene_Symbol), by = c("Signature" = "Mature1")) %>%
  mutate(Gene_conversion = ifelse(!is.na(Gene_Symbol), Gene_Symbol, Gene_conversion)) %>%
  select(-Gene_Symbol)  # Remove the temporary column

df45 <- df45 %>%
  left_join(miRNA_database %>% select(Mature2, Gene_Symbol), by = c("Signature" = "Mature2")) %>%
  mutate(Gene_conversion = ifelse(!is.na(Gene_Symbol), Gene_Symbol, Gene_conversion)) %>%
  select(-Gene_Symbol)  # Remove the temporary column

# Step 1: Extract gene symbols, transcript names, and miRNA from the "Signature" column
df50 <- df45 %>%
  # Use str_replace_all to remove parentheses and separate by "+" into individual elements
  mutate(Signature_cleaned = str_replace_all(Signature, "[()]", "")) %>%  # Remove parentheses
  separate_rows(Signature_cleaned, sep = " \\+ ") %>%  # Split by " + " into separate rows
  distinct() %>%  # Remove duplicate rows
  filter(str_detect(Signature_cleaned, "^ENST|^hsa-miR|^[A-Za-z]")) %>%  # Keep only gene symbols, ENST transcripts, or miRNA
  arrange(Signature_cleaned)  # Order alphabetically

# Step 2: Create the df50 with the unique and ordered entries
df50 <- df50 %>%
  select(Signature_cleaned) %>%
  rename(Extracted_Signature = Signature_cleaned)

# Step 1: Create a new column "Gene_Symbol" in df50 and populate it based on the mapping logic
df50 <- df50 %>%
  rowwise() %>%
  mutate(Gene_Symbol = case_when(
    
    # For values with prefix "ENST", look up in the "data" dataframe using "Transcript_ID"
    str_detect(Extracted_Signature, "^ENST") ~ 
      ifelse(Extracted_Signature %in% data$Transcript_ID, 
             data$Display_Name[which(data$Transcript_ID == Extracted_Signature)], 
             NA_character_),
    
    # For values with prefix "hsa-miR", look up in "miRNA_database" for Mature1
    str_detect(Extracted_Signature, "^hsa-miR") & Extracted_Signature %in% miRNA_database$Mature1 ~ 
      ifelse(Extracted_Signature %in% miRNA_database$Mature1, 
             miRNA_database$Gene_Symbol[which(miRNA_database$Mature1 == Extracted_Signature)], 
             NA_character_),
    
    # For values with prefix "hsa-miR", look up in "miRNA_database" for Mature2
    str_detect(Extracted_Signature, "^hsa-miR") & Extracted_Signature %in% miRNA_database$Mature2 ~ 
      ifelse(Extracted_Signature %in% miRNA_database$Mature2, 
             miRNA_database$Gene_Symbol[which(miRNA_database$Mature2 == Extracted_Signature)], 
             NA_character_),
    
    # If no match is found, leave as NA
    TRUE ~ NA_character_
  ))

# Fill NA values in "Gene_Symbol" with corresponding values from "Extracted_Signature"
df50 <- df50 %>%
  mutate(Gene_Symbol = ifelse(is.na(Gene_Symbol), Extracted_Signature, Gene_Symbol))

# Remove duplicates based on "Gene_Symbol"
df50 <- df50 %>%
  distinct(Gene_Symbol, .keep_all = TRUE)

# Order the dataframe by "Gene_Symbol" alphabetically
df50 <- df50 %>%
  arrange(Gene_Symbol)

# Exclude the "Extracted_Signature" column
df50 <- df50 %>%
  select(-Extracted_Signature)

Gene_Symbols_worst_signatures <- df50

# Order the dataframe "Gene_Symbols_top_signatures" alphabetically by the "Gene_Symbol" column
Gene_Symbols_worst_signatures <- Gene_Symbols_worst_signatures[order(Gene_Symbols_worst_signatures$Gene_Symbol), ]

# Remove leading and trailing spaces globally across all columns
Gene_Symbols_worst_signatures <- Gene_Symbols_worst_signatures %>%
  mutate_all(trimws)

rio::export(Gene_Symbols_worst_signatures, "Gene_Symbols_worst_signatures.tsv")

#####################
#####################
#####################
### Filtering specific signatures
# Select rows in df14 where "RCD_types" contains "/"
df_filtered_contains <- df14[grep("/", df14$RCD_types), ]

# Select rows in df14 where "RCD_types" does NOT contain "/"
df_filtered_does_not_contain <- df14[grep("/", df14$RCD_types, invert = TRUE), ]

# Select rows in df14 for specific features 
df_filtered_contain_ENST <- df14[grep("ENST", df14$Signature), ]
df_filtered_contain_miRNA <- df14[grep("hsa-miR", df14$Signature), ]

### Calculating the number of interactions; MUST REVISED IN DISCUSSION
number_of_pair_interactions <- ((5913*12*33)+(62090*3*33)+(882*3*33)+(239*3*33) # Radar
                                + (92733*4)+(96488*4)+(2970*4)+(287*4) # Cox
                                +(44641*4) # Survival
                                + (44641*119)) # TIL

### Calculating the number of interactions
number_of_pair_interactions_2 <- ((26766*12*33)+(16244*3*33)+(1470*3*33)+(258*3*33))
df97 <- number_of_pair_interactions
df99 <- number_of_pair_interactions_2
df100 <- df97 + df99
df100
# > df100
# [1] 27238756

############################### RPPA PROTEIN ARRAY MAPPING ##########################
################################ # Define the protein names and their corresponding gene symbols
protein_names <- c("ACC_pS79", "ACC1", "ACETYLATUBULINLYS40", "ACVRL1", "ADAR1", "AKT", "AKT_pS473", "AKT_pT308",
                   "ALPHACATENIN", "AMPKALPHA", "AMPKALPHA_pT172", "ANNEXIN1", "ANNEXINVII", "AR", "ARAF", "ARAF_pS299",
                   "ARID1A", "ASNS", "ATM", "AXL", "BAD_pS112", "BAK", "BAP1C4", "BAX", "BCL2", "BCL2A1", "BCLXL",
                   "BECLIN", "BETACATENIN", "BID", "BIM", "BRAF", "BRAF_pS445", "BRCA2", "BRD4", "CA9", "CABL",
                   "CASPASE3", "CASPASE7CLEAVEDD198", "CASPASE8", "CASPASE9", "CAVEOLIN1", "CD20", "CD26", "CD31",
                   "CD49B", "CDK1", "CDK1_pY15", "CHK1", "CHK1_pS296", "CHK1_pS345", "CHK2", "CHK2_pT68", "CHROMOGRANINANTERM",
                   "CIAP", "CJUN_pS73", "CK5", "CKIT", "CLAUDIN7", "CMET", "CMET_pY1235", "CMYC", "COG3", "COLLAGENVI",
                   "COMPLEXIISUBUNIT30", "CRAF", "CRAF_pS338", "CTLA4", "CYCLINB1", "CYCLIND1", "CYCLINE1", "CYCLINE2",
                   "DIRAS3", "DJ1", "DUSP4", "DVL3", "E2F1", "ECADHERIN", "EEF2", "EEF2K", "EGFR", "EGFR_pY1068", 
                   "EGFR_pY1173", "EIF4E", "EIF4G", "ENY2", "EPPK1", "ERALPHA", "ERALPHA_pS118", "ERCC1", "ERCC5",
                   "ERK2", "ETS1", "EZH2", "FASN", "FIBRONECTIN", "FOXM1", "FOXO3A", "FOXO3A_pS318S321", "G6PD",
                   "GAB2", "GAPDH", "GATA3", "GATA6", "GCN5L2", "GSK3_pS9", "GSK3ALPHABETA", "GSK3ALPHABETA_pS21S9",
                   "GYGGLYCOGENIN1", "GYS", "GYS_pS641", "HER2", "HER2_pY1248", "HER3", "HER3_pY1289", "HEREGULIN",
                   "HIF1ALPHA", "HSP70", "IGF1R_pY1135Y1136", "IGFBP2", "INPP4B", "IRF1", "IRS1", "JAB1", "JAK2", 
                   "JNK_pT183Y185", "JNK2", "KEAP1", "KU80", "LCK", "LCN2A", "LDHA", "LDHB", "LKB1", "MACC1", "MAPK_pT202Y204",
                   "MEK1", "MEK1_pS217S221", "MIG6", "MITOCHONDRIA", "MRE11", "MSH2", "MSH6", "MTOR", "MTOR_pS2448",
                   "MYH11", "MYOSINIIA", "MYOSINIIA_pS1943", "NAPSINA", "NCADHERIN", "NDRG1_pT346", "NF2", "NFKBP65_pS536",
                   "NOTCH1", "NRAS", "NRF2", "OXPHOSCOMPLEXVSUBUNITB", "P16INK4A", "P21", "P27", "P27_pT157", "P27_pT198",
                   "P38_pT180Y182", "P38MAPK", "P53", "P62LCKLIGAND", "P63", "P70S6K_pT389", "P70S6K1", "P90RSK",
                   "P90RSK_pT359S363", "PAI1", "PARP1", "PARPAB3", "PARPCLEAVED", "PAXILLIN", "PCADHERIN", "PCNA",
                   "PDCD1", "PDCD4", "PDK1", "PDK1_pS241", "PDL1", "PEA15", "PEA15_pS116", "PI3KP110ALPHA", "PI3KP85",
                   "PKCALPHA", "PKCALPHA_pS657", "PKCDELTA_pS664", "PKCPANBETAII_pS660", "PKM2", "PR", "PRAS40_pT246",
                   "PRDX1", "PREX1", "PTEN", "PYGB", "PYGBAB2", "PYGL", "PYGM", "RAB11", "RAB25", "RAD50", "RAD51",
                   "RAPTOR", "RB", "RB_pS807S811", "RBM15", "RET_pY905", "RICTOR", "RICTOR_pT1135", "S6", "S6_pS235S236",
                   "S6_pS240S244", "SCD1", "SETD2", "SF2", "SHC_pY317", "SHP2_pY542", "SLC1A5", "SMAC", "SMAD1", 
                   "SMAD3", "SMAD4", "SNAIL", "SRC", "SRC_pY416", "SRC_pY527", "STAT3_pY705", "STAT5ALPHA", "STATHMIN",
                   "SYK", "SYNAPTOPHYSIN", "TAZ", "TFRC", "THYMIDILATESYNTHASE", "TIGAR", "TRANSGLUTAMINASE", "TSC1",
                   "TTF1", "TUBERIN", "TUBERIN_pT1462", "VEGFR2", "X1433BETA", "X1433EPSILON", "X1433ZETA", "X4EBP1",
                   "X4EBP1_pS65", "X4EBP1_pT37T46", "X4EBP1_pT70", "X53BP1", "XBP1", "XRCC1", "YAP", "YAP_pS127", "YB1",
                   "YB1_pS102")

gene_symbols <- c("ACACA", "ACACA", "TUBA1A", "ACVRL1", "ADAR", "AKT1", "AKT1", "AKT1", 
                  "CTNNA1", "PRKAA1", "PRKAA1", "ANXA1", "ANXA7", "AR", "ARAF", "ARAF", 
                  "ARID1A", "ASNS", "ATM", "AXL", "BAD", "BAK1", "BAP1", "BAX", "BCL2", "BCL2A1", "BCL2L1", 
                  "BECN1", "CTNNB1", "BID", "BCL2L11", "BRAF", "BRAF", "BRCA2", "BRD4", "CA9", "ABL1", 
                  "CASP3", "CASP7", "CASP8", "CASP9", "CAV1", "MS4A1", "DPP4", "PECAM1", 
                  "ITGA2", "CDK1", "CDK1", "CHEK1", "CHEK1", "CHEK1", "CHEK2", "CHEK2", "CHGA", 
                  "BIRC2", "JUN", "KRT5", "KIT", "CLDN7", "MET", "MET", "MYC", "COG3", "COL6A1", 
                  "SDHB", "RAF1", "RAF1", "CTLA4", "CCNB1", "CCND1", "CCNE1", "CCNE2", 
                  "DIRAS3", "PARK7", "DUSP4", "DVL3", "E2F1", "CDH1", "EEF2", "EEF2K", "EGFR", "EGFR", 
                  "EGFR", "EIF4E", "EIF4G1", "ENY2", "EPPK1", "ESR1", "ESR1", "ERCC1", "ERCC5", 
                  "MAPK1", "ETS1", "EZH2", "FASN", "FN1", "FOXM1", "FOXO3", "FOXO3", "G6PD", 
                  "GAB2", "GAPDH", "GATA3", "GATA6", "KAT2A", "GSK3B", "GSK3B", "GSK3B", 
                  "GYG1", "GYS1", "GYS1", "ERBB2", "ERBB2", "ERBB3", "ERBB3", "NRG1", 
                  "HIF1A", "HSPA1A", "IGF1R", "IGFBP2", "INPP4B", "IRF1", "IRS1", "COPS5", "JAK2", 
                  "MAPK8", "MAPK9", "KEAP1", "XRCC5", "LCK", "LCN2", "LDHA", "LDHB", "STK11", "MACC1", "MAPK1", 
                  "MAP2K1", "MAP2K1", "ERRFI1", "MT-CO1", "MRE11A", "MSH2", "MSH6", "MTOR", "MTOR", 
                  "MYH11", "MYH9", "MYH9", "NAPSA", "CDH2", "NDRG1", "NF2", "RELA", 
                  "NOTCH1", "NRAS", "NFE2L2", "ATP5B", "CDKN2A", "CDKN1A", "CDKN1B", "CDKN1B", "CDKN1B", 
                  "MAPK14", "MAPK14", "TP53", "SQSTM1", "TP63", "RPS6KB1", "RPS6KB1", "RPS6KA1", 
                  "RPS6KA1", "SERPINE1", "PARP1", "PARP1", "PARP1", "PXN", "CDH3", "PCNA", 
                  "PDCD1", "PDCD4", "PDPK1", "PDPK1", "CD274", "PEA15", "PEA15", "PIK3CA", "PIK3R1", 
                  "PRKCA", "PRKCA", "PRKCD", "PRKCB", "PKM", "PGR", "AKT1S1", 
                  "PRDX1", "PREX1", "PTEN", "PYGB", "PYGB", "PYGL", "PYGM", "RAB11A", "RAB25", "RAD50", "RAD51", 
                  "RPTOR", "RB1", "RB1", "RBM15", "RET", "RICTOR", "RICTOR", "RPS6", "RPS6", 
                  "RPS6", "SCD", "SETD2", "SRSF1", "SHC1", "PTPN11", "SLC1A5", "DIABLO", "SMAD1", 
                  "SMAD3", "SMAD4", "SNAI1", "SRC", "SRC", "SRC", "STAT3", "STAT5A", "STMN1", 
                  "SYK", "SYP", "WWTR1", "TFRC", "TYMS", "TIGAR", "TGM2", "TSC1", 
                  "NKX2-1", "TSC2", "TSC2", "KDR", "YWHAB", "YWHAE", "YWHAZ", "EIF4EBP1", 
                  "EIF4EBP1", "EIF4EBP1", "EIF4EBP1", "TP53BP1", "XBP1", "XRCC1", "YAP1", "YAP1", "YBX1", 
                  "YBX1")

# Create a dataframe
df_protein_mapping <- data.frame(Protein = protein_names, Gene_Symbol = gene_symbols)

# Reorder columns in df3: Swap columns 1 and 3
# Swap columns 1 and 2 in df3
df_protein_mapping  <- df_protein_mapping [, c(2, 1)]

unique_protein_gene_symbols <-  as.data.frame(unique(df_protein_mapping$Gene_Symbol))

rio::export(df_protein_mapping, "TCGA_RPPA_array_protein_mapping_chatGPT.tsv")

###################
###################
## BRCA.antibody_annotation (RPPA)
###################
###################

BRCA_antibody_annotation <- import("BRCA_RPPA_array/BRCA.antibody_annotation.txt")

# Rename the columns in BRCA_antibody_annotation
colnames(BRCA_antibody_annotation)[colnames(BRCA_antibody_annotation) == "Gene Name"] <- "Gene_Symbol"
colnames(BRCA_antibody_annotation)[colnames(BRCA_antibody_annotation) == "Composite Element REF"] <- "Protein"

# Display the updated dataframe
print(colnames(BRCA_antibody_annotation))

# Perform the inner join as before
df32 <- merge(BRCA_antibody_annotation, df_protein_mapping, by.x = "Gene_Symbol", by.y = "Gene_Symbol", all = FALSE)

# Select rows in df_protein_mapping that were not mapped (not included in df32)
unmapped_df <- df_protein_mapping[!df_protein_mapping$Gene_Symbol %in% df32$Gene_Symbol, ]

# Perform an anti join to select rows in df_protein_mapping that were not mapped
# unmapped_df <- anti_join(df_protein_mapping, BRCA_antibody_annotation, by = c("Gene_Symbol" = "Gene name"))

# Perform an anti join on the shared "Gene_Symbol" column
unmapped_df <- anti_join(df_protein_mapping, BRCA_antibody_annotation, by = "Gene_Symbol")

# Select the "Gene_Symbol" and "Protein.y" columns from df32
df_selected <- df32[, c("Gene_Symbol", "Protein.y")]

# Rename the column "Protein.y" to "Protein"
colnames(df_selected)[colnames(df_selected) == "Protein.y"] <- "Protein"

# Ensure both dataframes have the same column names
colnames(unmapped_df) <- colnames(df_selected)

# Combine the two dataframes using rbind
RPPA_array_mapped <- rbind(df_selected, unmapped_df)

unique_RPPA_array_mapped <- as.data.frame(unique(RPPA_array_mapped))

unique_gene_symbols_RPPA_array_mapped  <- as.data.frame(unique(unique_RPPA_array_mapped$Gene_Symbol))

rio::export(unique_RPPA_array_mapped, "unique_RPPA_array_mapped.tsv")

Target_genes <-  import("Target_genes.csv")

# Count the number of unique values in the Gene_Symbol column
unique_gene_count <- length(unique(Target_genes$Gene))

# Create a new dataframe by filtering rows that contain the term "Apoptosis" in Pathway_Categorization
df_apoptosis <- Target_genes[grep("Apoptosis", Target_genes$Pathway_Categorization, ignore.case = TRUE), ]

# Count the number of unique apoptosis-related gene signatures
apoptosis_gene_count <- length(unique(df_apoptosis$Gene))

# Display the results
print(paste("Number of unique genes in the dataset: ", unique_gene_count))
print(paste("Number of unique apoptosis-related genes: ", apoptosis_gene_count))

# Display the filtered dataframe with apoptosis-related genes
df_apoptosis

# Apply trimws() to remove leading and trailing spaces from all columns
Target_genes <- Target_genes %>%
  mutate(across(everything(), ~trimws(.)))  # Applies trimws to all columns

# Rename the column 

colnames(Target_genes)[colnames(Target_genes) == "Gene"] <- "Gene_Symbol"

# Create df33 with rows where 'Gene_Symbol' matches between unique_RPPA_array_mapped and Target_genes
df33 <- merge(unique_RPPA_array_mapped, Target_genes, by = "Gene_Symbol")

RPPA_proteins_in_xena_mapping <- df33

rio::export(RPPA_proteins_in_xena_mapping, "RPPA_proteins_in_xena_mapping.tsv")

# Create df_non_matched with rows in unique_RPPA_array_mapped where 'Gene_symbol' does not match
# These are array protein genes not associated with RCD
# No Boolean Search evidence linking the provided genes to regulated cell death (RCD) 
df_not_matched <- unique_RPPA_array_mapped[!unique_RPPA_array_mapped$Gene_Symbol %in% Target_genes$Gene_Symbol, ]

rio::export(df_not_matched, "RPPA_proteins_NOT_MATCH_in_xena_mapping.tsv")

##################
##################
##################
##################
################## Counting elements within the signatures
##################
##################
################## 
# Adicionar uma coluna que conta o número de elementos em cada assinatura

df101 <-  final_signatures_full_ranking

# Ensure "Ranking" is numeric 
df101 <- df101 %>%
  mutate(Ranking = as.numeric(Ranking))  # Convert Ranking to numeric if it's not already

df101 <- df101 %>%
  mutate(Count_source = sapply(strsplit(Signature, " \\+ "), length)) %>%
  select(1:(which(colnames(df101) == "Nomenclature") - 1), 
         Count_source, 
         everything())

rio::export(df101, "final_signatures_full_ranking_distribution.tsv")

######
######
###### Signature counts per genomic feature
###### 
###### 
# Filtrar para manter apenas linhas onde a coluna genotipica seja "mRNA" or"Transcript"
rank <- df101

rank_filtered <- rank %>%
  filter(
    (Genotype == "mRNA" | Genotype == "Transcript"), 
    !(TNC %in% c("0", "1")),  # Exclude rows where TNC is "1" or "0"
    Type_Cox_OS != "NS", Type_Cox_DSS != "NS", Type_Cox_DFI != "NS", Type_Cox_PFI != "NS",  # Pelo menos uma métrica de Cox diferente de "NS"
    Type_log_rank_OS != "NS", Type_log_rank_DSS != "NS", Type_log_rank_DFI != "NS", Type_log_rank_PFI != "NS",  # Pelo menos uma métrica de log-rank diferente de "NS"
    microenvironment_classification != "NS",
    immune_classification != "NS"
  )

# Get unique genotypes
genotypes <- unique(df101$Genotype)

# Initialize a data frame to store results
results_genotypes <- data.frame(Genotype = character(), Count = integer(), stringsAsFactors = FALSE)

# Create data frames for each genotype and count rows
for (genotype in genotypes) {
  # Create a data frame name by replacing spaces with underscores
  df_name <- paste0("df_", gsub(" ", "_", genotype))
  
  # Assign the subset of df101 to the new data frame
  assign(df_name, df101[df101$Genotype == genotype, ])
  
  # Count the number of rows
  count <- nrow(get(df_name))
  
  # Append results to the results data frame
  results_genotypes <- rbind(results_genotypes, data.frame(Genotype = genotype, Count = count))
}

# Print the total number of counts
total_count <- sum(results_genotypes$Count)
cat("Total count of rows across all genotypes:", total_count, "\n")

rio::export(results_genotypes, "genotype_counts.xlsx")

# Write the results to an Excel file
write.xlsx(results_genotypes, "genotype_counts.xlsx", rowNames = FALSE)

# Calculate mean and median for the "Count_source" variable in df101
mean_value <- mean(df101$Count_source, na.rm = TRUE)
median_value <- median(df101$Count_source, na.rm = TRUE)

# Display the mean and median values
cat("Mean of Count_source:", mean_value, "\n")
cat("Median of Count_source:", median_value, "\n")

# Count the number of rows where "Count_source" is equal to 1
count_equal_1 <- sum(df101$Count_source == 1, na.rm = TRUE)

# Display the result
cat("Number of rows with Count_source equal to 1:", count_equal_1, "\n")
# Number of rows with Count_source equal to 1: 28673 
# > 28673*100/44641
# [1] 64.23019

# Calculate quartiles (0%, 25%, 50%, 75%, 100%)
quantiles <- quantile(df101$Count_source, probs = seq(0, 1, 0.25), na.rm = TRUE)

# Display the quartiles
print(quantiles)
# 0%  25%  50%  75% 100% 
# 1    1    1    2 2052 

# Optional: Summary of basic statistics (includes min, 1st quartile, median, mean, 3rd quartile, and max)
summary_stats <- summary(df101$Count_source)
print(summary_stats)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 1.000    1.000    1.000    4.316    2.000 2052.000 
# 
# Calculate deciles (0%, 10%, 20%, ..., 100%)
deciles <- quantile(df101$Count_source, probs = seq(0, 1, 0.1), na.rm = TRUE)

# Display the deciles
print(deciles)
# 0%  10%  20%  30%  40%  50%  60%  70%  80%  90% 100% 
# 1    1    1    1    1    1    1    2    3    6 2052

# Calculate percentiles (1%, 2%, ..., 100%)
percentiles <- quantile(df101$Count_source, probs = seq(0, 1, 0.01), na.rm = TRUE)

# Display the percentiles
print(percentiles)

# Boxplot to visualize the distribution
boxplot(df101$Count_source, main = "Boxplot of Count_source", ylab = "Count_source")

# Histogram for frequency distribution
hist(df101$Count_source, breaks = 50, main = "Histogram of Count_source", xlab = "Count_source")

# Create a boxplot
boxplot(df101$Count_source, 
        main = "Boxplot of Count_source", 
        ylab = "Count_source", 
        col = "lightblue", 
        outline = TRUE,  # Show outliers
        notch = TRUE)    # Add a notch to emphasize the median

# Calculate quartiles (Q1, Median (Q2), Q3)
quantiles <- quantile(df101$Count_source, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)

# Add horizontal lines for Q1, Median (Q2), and Q3
abline(h = quantiles[1], col = "red", lty = 2)   # Q1 (25th percentile)
abline(h = quantiles[2], col = "green", lty = 2) # Median (50th percentile)
abline(h = quantiles[3], col = "blue", lty = 2)  # Q3 (75th percentile)

# Add text labels for the quartiles
text(x = 1.2, y = quantiles[1], labels = paste("Q1 =", round(quantiles[1], 2)), col = "red", pos = 4)
text(x = 1.2, y = quantiles[2], labels = paste("Median =", round(quantiles[2], 2)), col = "green", pos = 4)
text(x = 1.2, y = quantiles[3], labels = paste("Q3 =", round(quantiles[3], 2)), col = "blue", pos = 4)

# Create a boxplot of Count_source
boxplot(df101$Count_source, 
        main = "Boxplot with Percentile Lines", 
        ylab = "Count_source", 
        col = "lightblue", 
        outline = TRUE,  # Show outliers
        notch = TRUE)    # Add a notch to emphasize the median

dev.off()

# Calculate percentiles (e.g., every 10%)
percentiles <- quantile(df101$Count_source, probs = seq(0, 1, 0.1), na.rm = TRUE)

# Add lines at each percentile
for (i in 1:length(percentiles)) {
  abline(h = percentiles[i], col = "gray", lty = 2)  # Add dashed gray lines for each percentile
}

# Optionally, add text labels for specific percentiles (e.g., 10th, 50th, and 90th)
text(x = 1.2, y = percentiles[2], labels = paste("10th =", round(percentiles[2], 2)), col = "red", pos = 4)
text(x = 1.2, y = percentiles[5], labels = paste("Median =", round(percentiles[5], 2)), col = "green", pos = 4)
text(x = 1.2, y = percentiles[9], labels = paste("90th =", round(percentiles[9], 2)), col = "blue", pos = 4)

# Apply log transformation to "Count_source" to handle the skewed distribution
df101$log_Count_source <- log10(df101$Count_source + 1)  # Adding 1 to avoid log(0)

# Create a boxplot of the log-transformed Count_source
boxplot(df101$log_Count_source, 
        main = "Boxplot of Log-Transformed Count_source with Percentile Lines", 
        ylab = "Log(Count_source)", 
        col = "lightblue", 
        outline = TRUE, 
        notch = TRUE)

dev.off()

# Calculate percentiles on the log-transformed data (e.g., every 10%)
log_percentiles <- quantile(df101$log_Count_source, probs = seq(0, 1, 0.1), na.rm = TRUE)

# Add lines at each percentile
for (i in 1:length(log_percentiles)) {
  abline(h = log_percentiles[i], col = "gray", lty = 2)  # Add dashed gray lines for each percentile
}

# Optionally, add text labels for specific percentiles
text(x = 1.2, y = log_percentiles[2], labels = paste("10th =", round(log_percentiles[2], 2)), col = "red", pos = 4)
text(x = 1.2, y = log_percentiles[5], labels = paste("Median =", round(log_percentiles[5], 2)), col = "green", pos = 4)
text(x = 1.2, y = log_percentiles[9], labels = paste("90th =", round(log_percentiles[9], 2)), col = "blue", pos = 4)

# Select rows where the variable "RCD_types" contains "/"
selected_rows <- df101[grep("/", df101$RCD_types), ]

# Count the number of rows that meet the condition
num_selected_rows <- nrow(selected_rows)

# Display the count
cat("Number of rows where RCD_types contains '/':", num_selected_rows, "\n")

# Total number of rows in df71 (excluding NAs in SCS)
total_rows <- sum(!is.na(df71$SCS))

# Count the number of rows where "SCS" is "N"
count_N <- sum(df71$SCS == "N", na.rm = TRUE)

# Count the number of rows where "SCS" is "P"
count_P <- sum(df71$SCS == "P", na.rm = TRUE)

# Calculate the percentages
percent_N <- (count_N / total_rows) * 100
percent_P <- (count_P / total_rows) * 100

# Display the counts and percentages
cat("Number of rows with 'N' in the 'SCS' column:", count_N, "(", round(percent_N, 2), "%)\n")
cat("Number of rows with 'P' in the 'SCS' column:", count_P, "(", round(percent_P, 2), "%)\n")

# Total number of rows in df71 (excluding NAs in RCD)
total_rows <- sum(!is.na(df71$RCD))

# Count the number of rows where "RCD" is "1"
count_4 <- sum(df71$RCD == "1", na.rm = TRUE)

# Calculate the percentages
percent_4 <- (count_4 / total_rows) * 100

# Display the counts and percentages
cat("Number of rows with 'N=1' in the 'RCD' column:", count_4, "(", round(percent_4, 2), "%)\n")

# Total number of rows in df71 (excluding NAs in TNC)
total_rows <- sum(!is.na(df71$TNC))

# Count the number of rows where "TNC" is "0"
count_0 <- sum(df71$TNC == "0", na.rm = TRUE)

# Count the number of rows where "TNC" is "1"
count_1 <- sum(df71$TNC == "1", na.rm = TRUE)

# Count the number of rows where "TNC" is "2"
count_2 <- sum(df71$TNC == "2", na.rm = TRUE)

# Count the number of rows where "TNC" is "3"
count_3 <- sum(df71$TNC == "3", na.rm = TRUE)

# Calculate the percentages
percent_0 <- (count_0 / total_rows) * 100
percent_1 <- (count_1 / total_rows) * 100
percent_2 <- (count_2 / total_rows) * 100
percent_3 <- (count_3 / total_rows) * 100

# Display the counts and percentages
cat("Number of rows with '0' in the 'TNC' column:", count_0, "(", round(percent_0, 2), "%)\n")
cat("Number of rows with '1' in the 'TNC' column:", count_1, "(", round(percent_1, 2), "%)\n")
cat("Number of rows with '2' in the 'TNC' column:", count_2, "(", round(percent_2, 2), "%)\n")
cat("Number of rows with '3' in the 'TNC' column:", count_3, "(", round(percent_3, 2), "%)\n")

# Total number of rows in df71 (excluding NAs in RCD)
total_rows <- sum(!is.na(df71$RCD))

# Count the number of rows where "RCD" is "1"
count_1 <- sum(df71$RCD == "1", na.rm = TRUE)

# Calculate the percentage
percent_1 <- (count_1 / total_rows) * 100

# Display the count and percentage
cat("Number of rows with 'RCD = 1' in the 'RCD' column:", count_1, "(", round(percent_1, 2), "%)\n")

# Select the rows where "RCD" is "1" and save them to df104
df104 <- df101[df101$RCD == "1", ]

# Count the number of rows where RCD_types is "Apoptosis"
apoptosis_count <- sum(df104$RCD_types == "Apoptosis")

# Display the count
cat("Number of rows with 'Apoptosis' in the 'RCD_types' column:", apoptosis_count, "\n")

# Select the rows where "RCD_types" is "Apoptosis" and save them to df105
df105 <- df104[df104$RCD_types == "Apoptosis", ]

# Ensure "Ranking" is numeric 
df105 <- df105 %>%
  mutate(Ranking = as.numeric(Ranking))  # Convert Ranking to numeric if it's not already

# Select the row with the highest "Ranking" and resolve ties progressively
df_highest_ranking_apoptosis <- df105 %>%
  group_by(Genotype) %>%
  arrange(desc(Ranking), desc(TIC), desc(SMC), desc(TMC), desc(HRC)) %>%  # Order by the tie-breaking variables
  slice_head(n = 1) %>%  # Select the first row after resolving ties
  ungroup()

df106 <- df_highest_ranking_apoptosis

df_highest_ranking_apoptosis_tnc <- df105 %>%
  filter(!(TNC %in% c(0, 1))) %>%  # Exclude rows where TNC_ranking is 0 or 1
  group_by(Genotype) %>%
  arrange(desc(Ranking), desc(TIC), desc(SMC), desc(TMC), desc(HRC)) %>%  # Order by the tie-breaking variables
  slice_head(n = 1) %>%  # Select the first row after resolving ties
  ungroup()

df106_tnc <- df_highest_ranking_apoptosis_tnc

# Select specified columns from df106, rename them, reorder them, and order by descending Rank
df107_tnc <- df106_tnc %>%
  select(Ranking, Nomenclature, Signature, Genotype, RCD_types) %>%
  rename(
    Rank = Ranking,
    `Genomic feature` = Genotype,
    `RCD forms` = RCD_types
  ) %>%
  select(Rank, Nomenclature, Signature, `Genomic feature`, `RCD forms`) %>%  # Reorder the columns
  arrange(desc(Rank))  # Order by descending Rank

rio::export(df107_tnc, "df_highest_ranking_apoptosis_tnc.xlsx")

# Select specified columns from df106, rename them, reorder them, and order by descending Rank
df107 <- df106 %>%
  select(Ranking, Nomenclature, Signature, Genotype, RCD_types) %>%
  rename(
    Rank = Ranking,
    `Genomic feature` = Genotype,
    `RCD forms` = RCD_types
  ) %>%
  select(Rank, Nomenclature, Signature, `Genomic feature`, `RCD forms`) %>%  # Reorder the columns
  arrange(desc(Rank))  # Order by descending Rank

rio::export(df107, "df_highest_ranking_apoptosis.xlsx")

# Ensure "Ranking" is numeric
df105 <- df105 %>%
  mutate(Ranking = as.numeric(Ranking))  # Convert Ranking to numeric if it's not already

# Step 1: Initial filtering and selection
df_highest_ranking_apoptosis_tnc_limit <- df105 %>%
  # Filter rows based on primary conditions
  filter(
    Ranking >= 20,           # "Ranking" must be >= 20
    !TNC %in% c(0, 1),       # "TNC" must not be 0 or 1
    !(TMC %in% c(4, 2)),     # "TMC" must not be 4 or 2
    !(TIC %in% c(4, 2))      # "TIC" must not be 4 or 2
  ) %>%
  # Group by Genotype
  group_by(Genotype) %>%
  # Sort within each group
  arrange(desc(Ranking), desc(TIC), desc(SMC), desc(TMC), desc(HRC)) %>%
  # Select the top-ranked row for each Genotype
  slice_head(n = 1) %>%
  # Ungroup
  ungroup()

# Step 2: Identify missing Genotype values
missing_genotypes <- setdiff(unique(df105$Genotype), unique(df_highest_ranking_apoptosis_tnc_limit$Genotype))

# Step 3: Handle missing Genotype values
if (length(missing_genotypes) > 0) {
  # Filter the original data for rows with the missing Genotype values
  additional_rows <- df105 %>%
    filter(
      Genotype %in% missing_genotypes,  # Include rows with missing Genotype values
      !TNC %in% c(0, 1),               # Maintain original TNC filter
      TMC == "2",             # Maintain original TMC filter
      TIC == "2"  # Allow TIC == "Variable" as an exception
    ) %>%
    # Group by Genotype
    group_by(Genotype) %>%
    # Sort by Ranking and other tie-breaking variables
    arrange(desc(Ranking), desc(TIC == "2"), desc(TIC), desc(SMC), desc(TMC), desc(HRC)) %>%
    # Select the top-ranked row for each missing Genotype
    slice_head(n = 1) %>%
    # Ungroup
    ungroup()
  
  # Combine the original and additional rows
  df_highest_ranking_apoptosis_tnc_limit <- bind_rows(
    df_highest_ranking_apoptosis_tnc_limit,
    additional_rows
  )
}

# Final result with one row per Genotype
df106_tnc_limit <- df_highest_ranking_apoptosis_tnc_limit

# Ensure "Ranking" is numeric
df104 <- df104 %>%
  mutate(Ranking = as.numeric(Ranking))  # Convert Ranking to numeric if it's not already

# Step 1: Initial filtering and selection
df104_df_highest_ranking_apoptosis_tnc_limit <- df104 %>%
  # Filter rows based on primary conditions
  filter(
    Ranking >= 20,           # "Ranking" must be >= 20
    !TNC %in% c(0, 1),       # "TNC" must not be 0 or 1
    !(TMC %in% c(4, 2)),     # "TMC" must not be 4 or 2
    !(TIC %in% c(4, 2))      # "TIC" must not be 4 or 2
  ) %>%
  # Group by RCD_types
  group_by(RCD_types) %>%
  # Sort within each group
  arrange(desc(Ranking), desc(TIC), desc(SMC), desc(TMC), desc(HRC)) %>%
  # Select the top-ranked row for each RCD_types
  slice_head(n = 1) %>%
  # Ungroup
  ungroup()

# Step 2: Identify missing RCD_types values
missing_rcd_types <- setdiff(unique(df104$RCD_types), unique(df104_df_highest_ranking_apoptosis_tnc_limit$RCD_types))

# Step 3: Handle missing RCD_types values
if (length(missing_rcd_types) > 0) {
  # Filter the original data for rows with the missing RCD_types values
  additional_rows <- df104 %>%
    filter(
      RCD_types %in% missing_rcd_types,  # Include rows with missing RCD_types values
      !TNC %in% c(0, 1),                # Maintain original TNC filter
      (TMC == "2" | TIC == "2")         # Relaxed filtering condition for TMC or TIC
    ) %>%
    # Group by RCD_types
    group_by(RCD_types) %>%
    # Sort by Ranking and other tie-breaking variables
    arrange(desc(Ranking), desc(TIC == "2"), desc(TIC), desc(SMC), desc(TMC), desc(HRC)) %>%
    # Select the top-ranked row for each missing RCD_types
    slice_head(n = 1) %>%
    # Ungroup
    ungroup()
  
  # Combine the original and additional rows
  df104_df_highest_ranking_apoptosis_tnc_limit <- bind_rows(
    df104_df_highest_ranking_apoptosis_tnc_limit,
    additional_rows
  )
}

# Deduplicate to ensure one row per RCD_types
df104_df_highest_ranking_apoptosis_tnc_limit <- df104_df_highest_ranking_apoptosis_tnc_limit %>%
  distinct(RCD_types, .keep_all = TRUE)

# Final result with one row per RCD_types
df104_selected_RCD_20 <- df104_df_highest_ranking_apoptosis_tnc_limit

RCD_types <- import("RCD_types.xlsx")

# Update specific values in the "RCD_types" column of RCD_types
RCD_types <- RCD_types %>%
  mutate(RCD_types = case_when(
    RCD_types == "Immunogenic Cell Death" ~ "Immunogenic cell death",
    RCD_types == "Mitotic Catastrophe" ~ "Mitotic catastrophe",
    RCD_types == "Cellular Senescence" ~ "Cellular senescence",
    RCD_types == "Mitochondrial Permeability Transition" ~ "Mitochondrial permeability transition",
    RCD_types == "Lysosome-Dependent Cell Death (LDCD)" ~ "Lysosome-dependent cell death",
    TRUE ~ RCD_types  # Keep other values unchanged
  ))

# Verify the changes
print(unique(RCD_types$RCD_types))

df104_selected_RCD_20 <- df104_selected_RCD_20 %>%
  mutate(RCD_types = str_trim(RCD_types))

# Find values in "RCD_types" of df104 that are absent in df104_selected_RCD_20
absent_rcd_types <- as.data.frame(setdiff(unique(RCD_types$RCD_types), unique(df104_selected_RCD_20$RCD_types)))

# Print the absent values
print(absent_rcd_types)

#######
####### Multi-modular RCD analysis
#######

# Step 1: Select rows where "RCD_types" contains "/"
selected_rows <- df101 %>%
  filter(str_detect(RCD_types, "/"))

# Step 2: Further filter rows where "RCD_types" contains any of the specified terms
RCD_missing_filtered_rows <- selected_rows %>%
  filter(str_detect(RCD_types, "Autosis|Lysosome-dependent cell death|Paraptosis|Erebosis|Methuosis"))

# Ensure "Ranking" is numeric
RCD_missing_filtered_rows <- RCD_missing_filtered_rows %>%
  mutate(Ranking = as.numeric(Ranking))  # Convert Ranking to numeric if it's not already

# Step 1: Initial filtering and selection
filtered_RCD_highest_ranking <- RCD_missing_filtered_rows %>%
  # Filter rows based on primary conditions
  filter(
    Ranking >= 20,           # "Ranking" must be >= 20
    !TNC %in% c(0, 1),       # "TNC" must not be 0 or 1
    !(TMC %in% c(4, 2)),     # "TMC" must not be 4 or 2
    !(TIC %in% c(4, 2))      # "TIC" must not be 4 or 2
  ) %>%
  # Group by RCD_types
  group_by(RCD_types) %>%
  # Sort within each group
  arrange(desc(Ranking), desc(TIC), desc(SMC), desc(TMC), desc(HRC)) %>%
  # Select the top-ranked row for each RCD_types
  slice_head(n = 1) %>%
  # Ungroup
  ungroup()

# Step 2: Identify missing RCD_types values
missing_rcd_types <- setdiff(
  unique(RCD_missing_filtered_rows$RCD_types),
  unique(filtered_RCD_highest_ranking$RCD_types)
)

# Step 3: Handle missing RCD_types values
if (length(missing_rcd_types) > 0) {
  # Filter the original data for rows with the missing RCD_types values
  additional_rows <- RCD_missing_filtered_rows %>%
    filter(
      RCD_types %in% missing_rcd_types,  # Include rows with missing RCD_types values
      Ranking >= 20,                    # Reapply the "Ranking >= 20" condition
      !TNC %in% c(0, 1),                # Maintain original TNC filter
      !(TMC %in% c(4, 2)),              # Exclude TMC equal to 4 or 2
      !(TIC %in% c(4, 2))               # Exclude TIC equal to 4 or 2
    ) %>%
    # Group by RCD_types
    group_by(RCD_types) %>%
    # Sort by Ranking and other tie-breaking variables
    arrange(desc(Ranking), desc(TIC), desc(SMC), desc(TMC), desc(HRC)) %>%
    # Select the top-ranked row for each missing RCD_types
    slice_head(n = 1) %>%
    # Ungroup
    ungroup()
  
  # Combine the original and additional rows
  filtered_RCD_highest_ranking <- bind_rows(
    filtered_RCD_highest_ranking,
    additional_rows
  )
}

# Deduplicate to ensure one row per RCD_types
filtered_RCD_highest_ranking <- filtered_RCD_highest_ranking %>%
  distinct(RCD_types, .keep_all = TRUE)

# Final result with one row per RCD_types
final_RCD_filtered <- filtered_RCD_highest_ranking

# Verification: Ensure no Ranking < 20
ranking_check <- all(final_RCD_filtered$Ranking >= 20)

# Print verification result
print(ranking_check)  # TRUE if no rows with Ranking < 20 exist, FALSE otherwise

# Verification: Ensure expected values are in the final result
expected_terms <- c("Autosis", "Lysosome-dependent cell death", "Paraptosis", "Erebosis", "Methuosis")
verification_result <- all(sapply(expected_terms, function(term) any(str_detect(final_RCD_filtered$RCD_types, term))))

# Print verification result
print(verification_result)  # TRUE if all expected terms are present, FALSE otherwise

# Ensure "Ranking" is numeric (if not already)
final_RCD_filtered <- final_RCD_filtered %>%
  mutate(Ranking = as.numeric(Ranking))  # Ensure "Ranking" is numeric

# Initialize an empty dataframe to store the final result
final_unique_terms <- data.frame()

# Iterate over expected terms and select the higher-ranked row for each term
for (term in expected_terms) {
  # Filter rows containing the specific term
  term_rows <- final_RCD_filtered %>%
    filter(str_detect(RCD_types, term))
  
  # Select the row with the highest ranking for this term
  top_row <- term_rows %>%
    arrange(desc(Ranking)) %>%
    slice_head(n = 1)
  
  # Append the selected row to the final dataframe
  final_unique_terms <- bind_rows(final_unique_terms, top_row)
}

# Ensure there are no overlapping terms in the final result
final_unique_terms <- final_unique_terms %>%
  distinct(RCD_types, .keep_all = TRUE)

#####
#####
# Check if names and order of variables are the same
if (identical(names(df104_selected_RCD_20), names(final_unique_terms))) {
  # If they are the same, combine the dataframes
  RCD_TOP_selected_comprehensive_signatures <- rbind(df104_selected_RCD_20, final_unique_terms)
  # Print success message and preview the combined dataframe
  print("Variable names and order match. Dataframes have been successfully combined.")
  print(head(RCD_TOP_selected_comprehensive_signatures))
} else {
  # If they are not the same, print an error message
  print("Variable names and/or order do not match. Please align them before combining.")
}

rio::export(RCD_TOP_selected_comprehensive_signatures, "RCD_TOP_selected_comprehensive_signatures.xlsx")

#######
#######
#######
####### Integrating cbind full complete signatures
# Select rows from df101 where "Nomenclature" matches the "Signature_nomenclature" in "final_signature_ranking"
complete_final_signature_ranking <- df101 %>%
  filter(Nomenclature %in% final_signature_ranking$Signature_nomenclature)

df1011 <- complete_final_signature_ranking

# Specify the columns to check
columns_to_check <- c(1, 3, 6:8, 10:15)

# Check for "NA" or "NS" values in the specified columns
na_ns_check <- df1011 %>%
  select(all_of(columns_to_check)) %>%
  summarise_all(~ any(. %in% c("NA", "NS", "")))

# Print the results
print(na_ns_check)

# Convert the specified columns to numeric
df1011 <- df1011 %>%
  mutate(across(c(1, 3, 6:8, 10:15), as.numeric))

# Preview the updated dataframe
str(head(df1011))

###################
################### call df71
df71

# Specify the columns to check
columns_to_check <- c(3:5, 7:22)

# Check for "NA" or "NS" values in the specified columns
na_ns_check <- df71 %>%
  select(all_of(columns_to_check)) %>%
  summarise_all(~ any(. %in% c("NA", "NS", "")))

# Print the results
print(na_ns_check)

# Convert the specified columns to numeric
df71 <- df71 %>%
  mutate(across(c(3:5, 7:22), as.numeric))

# Preview the updated dataframe
str(head(df71))

# Order df1011 by "Nomenclature", "CTAB", and "GSI" in descending order
df1011 <- df1011 %>%
  arrange(desc(Nomenclature), desc(CTAB), desc(GSI))

# Order df71 by "Nomenclature", "CTAB", and "GSI" in descending order
df71 <- df71 %>%
  arrange(desc(Signature_nomenclature), desc(CTAB), desc(GSI))

# Select variables 13:22 from df71
df71_subset <- df71 %>%
  select(Signature_nomenclature, all_of(13:22))  # Include the key column for joining

# Join df1011 with df71_subset by the common key
df1011_updated <- df1011 %>%
  left_join(df71_subset, by = c("Nomenclature" = "Signature_nomenclature")) %>%
  relocate(all_of(names(df71_subset)[-1]), .after = 17)  # Place new variables after column 17

####
####
# sub-setting by genotype values
# Create and save separate dataframes for each unique "Genotype"
unique_genotypes <- unique(df1011_updated$Genotype)

# Loop through each unique genotype and save as an .xlsx file
for (genotype in unique_genotypes) {
  # Filter the dataframe for the current genotype
  genotype_df <- df1011_updated %>% filter(Genotype == genotype)
  
  # Create the file name
  file_name <- paste0("Dataset S1_", genotype, ".xlsx")
  
  # Save the dataframe to an .xlsx file
  write.xlsx(genotype_df, file_name)
  
  # Print the name of the file being created
  print(paste("File created:", file_name))
}

# Print message when done
print("All files created successfully!")

#######
#######
####### Table 2 Top-ranked RCD type-specific signatures with comprehensive multi-omic representation
#######
#######
# Ensure each Genotype is represented at least once and the most informative rows for each RCD_type are selected
df_representative <- df104 %>%
  group_by(RCD_types) %>%  # Group by unique values in "RCD_types"
  arrange(desc(Ranking), desc(Count_source), desc(TIC), desc(TMC), desc(SMC), desc(HRC)) %>%  # Order by Ranking, Count_source, TIC, TMC, SMC, HRC
  slice_head(n = 1) %>%  # Select the first row in each group (after resolving ties)
  ungroup()

# Now ensure that each unique "Genotype" is represented at least once
# We will check for any missing "Genotype" and select the highest-ranking row for those Genotypes
missing_genotypes <- setdiff(unique(df104$Genotype), unique(df_representative$Genotype))

# If there are missing genotypes, add one row for each missing genotype
if (length(missing_genotypes) > 0) {
  additional_rows <- df104 %>%
    filter(Genotype %in% missing_genotypes) %>%
    arrange(desc(Ranking), desc(Count_source), desc(TIC), desc(TMC), desc(SMC), desc(HRC)) %>%
    group_by(Genotype) %>%
    slice_head(n = 1) %>%
    ungroup()
  
  # Combine the rows that were already selected with the additional rows
  df_representative <- bind_rows(df_representative, additional_rows)
}

df108 <- df_representative

# Select specified columns from df106, rename them, reorder them, and order by descending Rank
df109 <- df108 %>%
  select(Ranking, Nomenclature, Signature, Genotype, RCD_types) %>%
  rename(
    Rank = Ranking,
    `Genomic feature` = Genotype,
    `RCD forms` = RCD_types
  ) %>%
  select(Rank, Nomenclature, Signature, `Genomic feature`, `RCD forms`) %>%  # Reorder the columns
  arrange(desc(Rank))  # Order by descending Rank


rio::export(df109, "RCD_specific_per_genotype_class_selected_signatures.xlsx")

# Select the rows where "RCD" is NOT "1" ((i.e,, more than one RCD form))
df110 <- df101[df101$RCD != "1", ]

# Ensure each Genotype is represented at least once and the most informative rows for each RCD_type are selected
df_representative_multi_RCD <- df110 %>%
  group_by(RCD_types) %>%  # Group by unique values in "RCD_types"
  arrange(desc(TIC), desc(TMC), desc(SMC), desc(HRC), desc(RCD), desc(Ranking), desc(Count_source)) %>%  # Order by these variables 
  slice_head(n = 1) %>%  # Select the first row in each group (after resolving ties)
  ungroup()

# Now ensure that each unique "Genotype" is represented at least once
# We will check for any missing "Genotype" and select the highest-ranking row for those Genotypes
missing_genotypes_multi_RCD <- setdiff(unique(df110$Genotype), unique(df_representative_multi_RCD$Genotype))

# If there are missing genotypes, add one row for each missing genotype
if (length(missing_genotypes_multi_RCD) > 0) {
  additional_rows_RCD <- df110 %>%
    filter(Genotype %in% missing_genotypes_multi_RCD) %>%
    arrange(desc(TIC), desc(TMC), desc(SMC), desc(HRC), desc(RCD), desc(Ranking), desc(Count_source)) %>% 
    group_by(Genotype) %>%
    slice_head(n = 1) %>%
    ungroup()
  
  # Combine the rows that were already selected with the additional rows
  df_representative_multi_RCD <- bind_rows(df_representative_multi_RCD, additional_rows_RCD)
}

df112 <- df_representative_multi_RCD

##### 
##### Filtering and Formatting Multi-RCD Signatures by Genotype with Ranking Criteria
#####
# Select specified columns, rename, reorder, filter, and sort; ranking >= 20
df113 <- df112 %>%
  select(Ranking, Nomenclature, Signature, Count_source, Genotype, RCD, RCD_types) %>%
  rename(
    Rank = Ranking,
    Elements = Count_source,
    `Genomic feature` = Genotype,
    `RCD count` = RCD,
    `RCD forms` = RCD_types
  ) %>%
  filter(Rank >= 20) %>%  # Filter for rows where Ranking (now Rank) is >= 20
  select(Rank, Nomenclature, Signature, Elements, `Genomic feature`, `RCD count`, `RCD forms`) %>%  # Reorder columns
  arrange(desc(`RCD count`), desc(Elements), desc(Rank))  # Sort in descending order

# rio::export the filtered and processed dataframe to an Excel file
rio::export(df113, "Multi_RCD_specific_per_genotype_class_selected_signatures.xlsx")

# Load required libraries
library(dplyr)
library(openxlsx)

# Order rows by Elements, RCD count, and Rank, ensuring 7 unique Genotypes in the top 20
df113_ordered <- df113 %>%
  arrange(desc(Elements), desc(`RCD count`), desc(Rank)) %>%
  group_by(`Genomic feature`) %>%  # Group by Genotype
  mutate(Genotype_rank = row_number()) %>%  # Assign a rank within each Genotype group
  ungroup() %>%
  arrange(Genotype_rank, desc(Elements), desc(`RCD count`), desc(Rank)) %>%
  mutate(Top20_priority = ifelse(row_number() <= 20, TRUE, FALSE))  # Flag rows in the top 20

# Extract the top 20 rows
df113_top20 <- df113_ordered %>%
  filter(Top20_priority)

# Save the reordered dataframe to an Excel file
rio::export(df113_ordered, "Multi_RCD_specific_per_genotype_class_selected_signatures_with_Priorities.xlsx")

# Save the top 20 rows to a separate Excel file
rio::export(df113_top20, "Top20_Multi_RCD_specific_per_genotype_class_selected_signatures_Genotype_Prioritized_Rows.xlsx")

# Order rows by Elements, RCD count, and Rank, ensuring 7 unique Genotype values are in the top 20
df113_ordered <- df113 %>%
  arrange(desc(Elements), desc(`RCD count`), desc(Rank)) %>%  # Initial ordering
  group_by(`Genomic feature`) %>%  # Group by Genotype
  mutate(Genotype_rank = row_number()) %>%  # Assign a rank within each Genotype group
  ungroup() %>%
  arrange(Genotype_rank, desc(Elements), desc(`RCD count`), desc(Rank)) %>%  # Reorder rows to prioritize Genotypes
  mutate(Top20_priority = ifelse(row_number() <= 20, TRUE, FALSE))  # Flag rows in the top 20

# View the first 20 rows to ensure 7 unique Genotype values are included
df113_top20 <- df113_ordered %>%
  filter(Top20_priority)

#####
#####
#####
##### GLOBAL df1011 DATA FILE -Filtering and rio::exporting Top 20 Genotype-Specific RCD Signatures Based on Priority Metrics"
##### 
#####
#####

##### Manuscript Table 1. top 20 Selected RCD multi-modular signatures with comprehensive multi-omic representation
##### Table 1. Top-ranked multi-modular RCD signatures with comprehensive multi-omic representation

# Filter and order rows where variables 6:15 and 22:25 are different from "0"
df1011_updated_ordered <- df1011_updated %>%
  filter(if_all(c(6:15, 22:25), ~ . != "0")) %>%  # Filter rows where all values in columns 6:15 and 22:25 are not "0"
  arrange(desc(RCD), desc(Count_source), desc(Ranking)) %>%
  group_by(Genotype) %>%
  mutate(Genotype_Ranking = row_number()) %>%  # Assign a Ranking within each Genotype group
  ungroup() %>%
  arrange(Genotype_Ranking, desc(Count_source), desc(RCD), desc(Ranking))

# Print success message
cat("Filtered rows where variables 6:15 and 22:25 are different from '0', and ordered successfully.\n")

# Ensure at least one row per Genotype in the top 20 with Ranking >= 20
df1011_top_genotypes <- df1011_updated_ordered %>%
  filter(Ranking >= 20) %>%
  group_by(Genotype) %>%
  slice_head(n = 1) %>%  # Select the highest-ranked row per Genotype
  ungroup()

# Exclude rows already in df1011_top_genotypes by matching unique identifiers (e.g., Genotype or row contents)
df1011_remaining_rows <- df1011_updated_ordered %>%
  filter(Ranking >= 20, !(Nomenclature %in% df1011_top_genotypes$Nomenclature)) %>%
  head(20 - nrow(df1011_top_genotypes))  # Fill up remaining rows to reach 20

# Combine to create the final top 20 dataframe
df1011_updated_top20 <- bind_rows(df1011_top_genotypes, df1011_remaining_rows) %>%
  arrange(desc(RCD), desc(Count_source), desc(Ranking))  # Reorder final top 20

# Select columns 1:4 and 15:17
df1011_select_updated_top20 <- df1011_updated_top20 %>%
  select(1:4, 15:17)

# Save the top 20 rows to a separate Excel file
rio::export(df1011_select_updated_top20, "Top20_Multi_RCD_specific_per_genotype_class_selected_signatures_Genotype_Prioritized_Rows.xlsx")

#####
#####
##### Table 1. Selected RCD multi-modular signatures with comprehensive multi-omic representation
##### 
##### Revised Table 1. Top-ranked multi-modular RCD signatures with comprehensive multi-omic representation
##### 
##### TOP 20 SIGNATURES FOR TABLE 1
#####
#####
# Step 1: Select at least two rows per Genotype, filtering for RCD >= 2 and RCD != 1
top_ranked_signatures_df1011 <- df1011_updated_ordered %>%
  filter(RCD >= 2 & RCD != 1) %>%  # Apply RCD filter
  group_by(Genotype) %>%
  filter(
    (Genotype %in% c("Protein", "miRNA") & TMC_ranking == 7 & TIC_ranking == 7 & SMC_ranking != 0) |  # Relax SMC_ranking (non-zero) for Protein and miRNA
      (SMC_ranking == 11 & TMC_ranking == 7 & TIC_ranking == 7)  # Strict conditions for all other Genotypes
  ) %>%
  slice_head(n = 2) %>%  # Select up to two rows per Genotype
  ungroup()

# Step 2: Identify missing Genotypes
missing_genotypes <- setdiff(unique(df1011_updated_ordered$Genotype), unique(top_ranked_signatures_df1011$Genotype))

# Step 3: Add rows for missing Genotypes, with the RCD filter applied
missing_genotype_rows <- df1011_updated_ordered %>%
  filter(
    RCD >= 2 & RCD != 1 &  # Apply RCD filter
      Genotype %in% missing_genotypes & 
      ((Genotype %in% c("Protein", "miRNA") & TMC_ranking == 7 & TIC_ranking == 7 & SMC_ranking != 0) |  # Relax SMC_ranking (non-zero) for Protein and miRNA
         (SMC_ranking == 11 & TMC_ranking == 7 & TIC_ranking == 7))  # Strict conditions for others
  ) %>%
  group_by(Genotype) %>%
  slice_head(n = 2) %>%  # Select the top rows for each missing Genotype
  ungroup()

# Combine initial selection and missing Genotypes
combined_df1011 <- bind_rows(top_ranked_signatures_df1011, missing_genotype_rows)

# Filter strictly for filling remaining rows
remaining_rows_df1011 <- df1011_updated_ordered %>%
  filter(
    RCD >= 2 & RCD != 1 &  # Apply RCD filter
      SMC_ranking == 11 & 
      TMC_ranking == 7 & 
      TIC_ranking == 7 & 
      !Nomenclature %in% combined_df1011$Nomenclature  # Avoid duplicates
  ) %>%
  arrange(desc(RCD), desc(Count_source), desc(Ranking))

# Check if strict filling yields 20 rows
if (nrow(combined_df1011) + nrow(remaining_rows_df1011) < 20) {
  # Relax SMC_ranking for additional filling
  relaxed_rows_df1011 <- df1011_updated_ordered %>%
    filter(
      RCD >= 2 & RCD != 1 &  # Apply RCD filter
        SMC_ranking != 0 &     # Relax SMC_ranking but exclude 0
        TMC_ranking == 7 & 
        TIC_ranking == 7 & 
        !Nomenclature %in% combined_df1011$Nomenclature  # Avoid duplicates
    ) %>%
    arrange(desc(RCD), desc(Count_source), desc(Ranking))
  
  # Add rows from relaxed filtering
  remaining_rows_df1011 <- bind_rows(
    remaining_rows_df1011,
    relaxed_rows_df1011[1:(20 - nrow(combined_df1011) - nrow(remaining_rows_df1011)), ]
  )
}

# step 4.Combine all rows to form the final dataframe
final_top_ranked_signatures_df1011 <- bind_rows(combined_df1011, remaining_rows_df1011) %>%
  arrange(desc(RCD), desc(Count_source), desc(Ranking))

cat("Row count after combining initial and missing Genotypes:", nrow(combined_df1011), "\n")
cat("Remaining rows available for filling:", nrow(remaining_rows_df1011), "\n")
cat("Missing Genotypes:\n")
print(setdiff(expected_genotypes, unique(final_top_ranked_signatures_df1011$Genotype)))

# Ensure the total number of rows is exactly 20
stopifnot(nrow(final_top_ranked_signatures_df1011) == 20)

# Ensure all 7 Genotypes are represented
stopifnot(length(unique(final_top_ranked_signatures_df1011$Genotype)) == 7)

# Select columns 1:4 and 15:17
final_selected_columns_df1011 <- final_top_ranked_signatures_df1011 %>%
  select(1:4, 15:17)

# Rename, reorder, and sort the columns in descending order by Rank
final_renamed_ordered_df1011 <- final_selected_columns_df1011 %>%
  rename(
    Rank = Ranking,
    Nomenclature = Nomenclature,
    Signature = Signature,
    Elements = Count_source,
    `Multi-omic feature` = Genotype,
    `RCD count` = RCD,
    `RCD forms` = RCD_types
  ) %>%
  select(Rank, Nomenclature, Signature, Elements, `Multi-omic feature`, `RCD count`, `RCD forms`) %>%  # Reorder columns
  arrange(desc(Rank))  # Sort in descending order by Rank

# Save the renamed and reordered dataframe to an Excel file
rio::export(final_renamed_ordered_df1011, "Renamed_Ordered_Top_20_Ranked_Signatures_df1011.xlsx")

# Print success message
cat("File saved successfully: Renamed_Ordered_Top_20_Ranked_Signatures_df1011.xlsx\n")

#####
#####
###### Table RCD-specific signatures
###### Table 2. Top-ranked signatures by comprehensive RCD type-specific 
###### and multi-omic feature representation.
# Select the rows where "RCD" is "1" and save them to df104
df1044 <- df1011_updated[df1011_updated$RCD == "1", ]

df1044_highest_ranking_tnc <- df1044 %>%
  filter(!(TNC %in% c(0, 1))) %>%  # Exclude rows where TNC_ranking is 0 or 1
  group_by(Genotype) %>%
  arrange(desc(Ranking), desc(TIC), desc(SMC), desc(TMC), desc(HRC)) %>%  # Order by the tie-breaking variables
  ungroup()

df1044_updated_ordered <- df1044_highest_ranking_tnc%>%
  filter(if_all(c(6:15, 22:25), ~ . != "0")) %>%  # Filter rows where all values in columns 6:15 and 22:25 are not "0"
  arrange(desc(RCD), desc(Count_source), desc(Ranking)) %>%
  group_by(Genotype) %>%
  mutate(Genotype_Ranking = row_number()) %>%  # Assign a Ranking within each Genotype group
  ungroup() %>%
  arrange(Genotype_Ranking, desc(Count_source), desc(RCD), desc(Ranking))

# Step 1: Identify the possible values for RCD_types in df1044
all_possible_rcd_types <- unique(df1044$RCD_types)

# Step 2: Identify the RCD_types present in df1044_updated_ordered
present_rcd_types <- unique(df1044_updated_ordered$RCD_types)

# Step 3: Find the missing RCD_types
missing_rcd_types <- setdiff(all_possible_rcd_types, present_rcd_types)
cat("Missing RCD_types:\n")
print(missing_rcd_types)

# Step 4: Filter df1044 to select rows with missing RCD_types
df1045 <- df1044 %>%
  filter(RCD_types %in% missing_rcd_types)

# Print a success message
cat("Filtered rows with missing RCD_types into df1045.\n")

df1045_select <- df1045 %>%
  filter(if_all(c(6:15), ~ . != "0")) %>%  # Filter rows where all values in columns 6:15 and 24:25 are not "0"
  arrange(desc(RCD), desc(Count_source), desc(Ranking)) %>%
  group_by(Genotype) %>%
  mutate(Genotype_Ranking = row_number()) %>%  # Assign a Ranking within each Genotype group
  ungroup() %>%
  arrange(Genotype_Ranking, desc(Count_source), desc(RCD), desc(Ranking))

df1045_select <- df1045 %>%
  filter(if_all(c(6:15), ~ . != "0")) %>%  # Filter rows where all values in columns 6:15 are not "0"
  arrange(desc(RCD), desc(Count_source), desc(Ranking)) %>%
  group_by(Genotype) %>%
  mutate(Genotype_Ranking = row_number()) %>%  # Assign a Ranking within each Genotype group
  ungroup() %>%
  arrange(Genotype_Ranking, desc(Count_source), desc(RCD), desc(Ranking)) %>%
  group_by(RCD_types) %>%
  slice_max(order_by = Ranking, n = 1, with_ties = FALSE) %>%  # Select the row with the highest Ranking per RCD_types
  ungroup()

# Step 1: Check if the variable names are the same and in the same order
if (identical(names(df1044_updated_ordered), names(df1045_select))) {
  # Step 2: Combine the dataframes using rbind
  df1046 <- rbind(df1044_updated_ordered, df1045_select)
  
  # Print success message and display the combined dataframe
  cat("Variable names match and dataframes have been combined into df1046.\n")
  
  # Optional: Save the combined dataframe to an Excel file
  library(openxlsx)
  write.xlsx(df1046, "Combined_df1046.xlsx")
  cat("Combined dataframe saved as 'Combined_df1046.xlsx'.\n")
} else {
  # Step 3: Print error message if variable names do not match
  cat("Error: Variable names do not match or are not in the same order. Cannot combine dataframes.\n")
  
  # Print variable names for debugging
  cat("Names in df1044_updated_ordered:\n")
  print(names(df1044_updated_ordered))
  
  cat("Names in df1045_select:\n")
  print(names(df1045_select))
}

# Step 2: Identify the RCD_types present in df1046
present_rcd_types_df1046 <- unique(df1046$RCD_types)

# Step 3: Find the missing RCD_types
missing_rcd_types_df1046 <- setdiff(all_possible_rcd_types, present_rcd_types_df1046)
cat("Missing RCD_types:\n")
print(missing_rcd_types_df1046)

#######
#######
####### Table 2. Top-ranked RCD type-specific signatures with comprehensive multi-omic representation reordered
#######
#######

# Step 1: Include single-instance Genotype rows
genotype_counts <- df1046 %>%
  count(Genotype)  # Count occurrences of each Genotype

single_instance_genotypes <- genotype_counts %>%
  filter(n == 1) %>%
  pull(Genotype)  # Identify Genotypes with a single occurrence

# Select rows with single-instance Genotypes
single_instance_rows <- df1046 %>%
  filter(Genotype %in% single_instance_genotypes)

# Step 2: Identify missing RCD_types and Genotypes
remaining_rcd_types <- setdiff(unique(df1046$RCD_types), unique(single_instance_rows$RCD_types))
remaining_genotypes <- setdiff(unique(df1046$Genotype), unique(single_instance_rows$Genotype))

# Step 3: Select rows to ensure at least one row per RCD_types
rcd_type_rows <- df1046 %>%
  filter(RCD_types %in% remaining_rcd_types) %>%
  arrange(desc(Ranking)) %>%
  group_by(RCD_types) %>%
  slice_max(order_by = Ranking, n = 1, with_ties = FALSE) %>%
  ungroup()

# Step 4: Select rows to ensure at least one row per Genotype
genotype_rows <- df1046 %>%
  filter(Genotype %in% remaining_genotypes) %>%
  arrange(desc(Ranking)) %>%
  group_by(Genotype) %>%
  slice_max(order_by = Ranking, n = 1, with_ties = FALSE) %>%
  ungroup()

# Combine all selected rows
combined_rows <- bind_rows(single_instance_rows, rcd_type_rows, genotype_rows) %>%
  distinct() %>%
  arrange(desc(Ranking))

# Step 5: Validation
stopifnot(all(unique(df1046$RCD_types) %in% combined_rows$RCD_types))  # Ensure all RCD_types are covered
stopifnot(all(unique(df1046$Genotype) %in% combined_rows$Genotype))  # Ensure all Genotypes are covered

RCD_specific_signatures <- combined_rows

# Select the first 70 columns from both dataframes
filtered_df_RCD_subset <- RCD_specific_signatures[, 1:70]
RCD_specific_signatures_subset <- RCD_specific_signatures[, 1:70]

# Check if the variable names and order are identical for the first 70 columns
if (identical(names(filtered_df_RCD_subset), names(RCD_specific_signatures_subset))) {
  # Perform rbind operation for the first 70 columns
  final_RCD_specific_signatures <- rbind(filtered_df_RCD_subset, RCD_specific_signatures_subset)
  
  # Print success message and show the combined dataframe
  cat("The first 70 columns are identical. Combined dataframe created:\n")
  print(final_RCD_specific_signatures)
} else {
  # Print error message and show differences in variable names or order for the first 70 columns
  cat("The first 70 column names or order are not identical. Differences:\n")
  cat("Names in filtered_df_RCD_subset:\n")
  print(names(filtered_df_RCD_subset))
  cat("Names in RCD_specific_signatures_subset:\n")
  print(names(RCD_specific_signatures_subset))
}

df1047 <- final_RCD_specific_signatures 
# Select variables c(1:4, 15:17) from df1047
selected_columns_df1047 <- df1047[, c(1:4, 15:17)]

# Rename, reorder, and sort the columns in descending order by Rank
selected_columns_df1047_reordered <- selected_columns_df1047 %>%
  rename(
    Rank = Ranking,
    Nomenclature = Nomenclature,
    Signature = Signature,
    Elements = Count_source,
    `Multi-omic feature` = Genotype,
    `RCD count` = RCD,
    `RCD forms` = RCD_types
  ) %>%
  select(Rank, Nomenclature, Signature, Elements, `Multi-omic feature`, `RCD count`, `RCD forms`) %>%  # Reorder columns
  arrange(desc(Rank))  # Sort in descending order by Rank


# rio::export the updated table to Excel
rio::export(selected_columns_df1047_reordered, "Table 2. Top-ranked RCD type-specific signatures with comprehensive multi-omic representation_reordered.xlsx")

#####
#####
#####
##### Table 3. Seven top-ranked signatures by multi-omic feature
#####
#####
##### 

#step 1
top_ranked_signatures_v3 <- df1011_updated[
  df1011_updated$SMC_ranking == 11 & 
    df1011_updated$TMC_ranking == 7 & 
    df1011_updated$TIC_ranking == 7, 
]

# Step 2: Identify the possible values for Genotypes 
all_possible_genotypes <- unique(df1011_updated$Genotype)

# Step 3: Select rows to ensure at least one row per Genotype
genotype_rows_df1011 <- top_ranked_signatures_v3 %>%
  filter(Genotype %in% remaining_genotypes) %>%
  arrange(desc(Ranking)) %>%
  group_by(Genotype) %>%
  slice_max(order_by = Ranking, n = 1, with_ties = FALSE) %>%
  ungroup()

# Step 4: Identify the Genotypes present in genotype_rows_df1011
present_genotypes_types <- unique(genotype_rows_df1011$Genotype)

# Step 5: Find the missing RCD_types
missing_genotypes <- setdiff(all_possible_genotypes, present_genotypes_types)
cat("Missing Genotypes:\n")
print(missing_genotypes)

# Step 6: Select the top row for the Genotype "Protein" with the highest values in HRC, and SMC
genotype_rows_protein <- df1011_updated %>%
  filter(Genotype == "Protein",  # Ensure the Genotype is "Protein"
         TMC_ranking == 7,       # Ensure TMC_ranking is 7
         TIC_ranking == 7) %>%   # Ensure TIC_ranking is 7
  arrange(desc(HRC), desc(SMC)) %>%  # Sort by Ranking first, then HRC, then SMC
  slice_head(n = 1)  # Select the top row

# Step 7: Select the top row for the Genotype "miRNA" with the highest values in Ranking, HRC, and SMC
genotype_rows_miRNA <- df1011_updated %>%
  filter(Genotype == "miRNA",           # Ensure the Genotype is "miRNA"
         !(TNC %in% c(0,1)),           # Exclude rows where TNC is 0 or 1
         TMC_ranking == 7,              # Ensure TMC_ranking is 7
         TIC_ranking == 7) %>%          # Ensure TIC_ranking is 7
  arrange(desc(Ranking), desc(HRC), desc(SMC)) %>%  # Sort by Ranking, HRC, and SMC
  slice_head(n = 1)                     # Select the top row

# step 8. Check if all variable names and order are identical across the three dataframes
# Check if all variable names and order are identical across the three dataframes
if (identical(names(genotype_rows_df1011), names(genotype_rows_miRNA)) &&
    identical(names(genotype_rows_df1011), names(genotype_rows_protein))) {
  
  # Perform rbind operation to combine the dataframes
  combined_genotype_rows <- rbind(genotype_rows_df1011, genotype_rows_miRNA, genotype_rows_protein) %>%
    arrange(desc(Ranking))  # Sort rows by descending Ranking
  
  # Print success message and show the combined dataframe
  cat("Variable names and order are identical. Combined dataframe created and sorted by Ranking:\n")
  print(combined_genotype_rows)
} else {
  # Print error message and show the differences in variable names or order
  cat("Variable names or order are not identical among the dataframes.\n")
  cat("Names in genotype_rows_df1011:\n")
  print(names(genotype_rows_df1011))
  cat("Names in genotype_rows_miRNA:\n")
  print(names(genotype_rows_miRNA))
  cat("Names in genotype_rows_protein:\n")
  print(names(genotype_rows_protein))
}

# Select variables c(1:4, 15:17) from df1047
selected_columns_combined_genotype_rows <- combined_genotype_rows[, c(1:4, 15:17)]

# Rename, reorder, and sort the columns in descending order by Rank
selected_columns_combined_genotype_rows <- selected_columns_combined_genotype_rows %>%
  rename(
    Rank = Ranking,
    Nomenclature = Nomenclature,
    Signature = Signature,
    Elements = Count_source,
    `Multi-omic feature` = Genotype,
    `RCD count` = RCD,
    `RCD forms` = RCD_types
  ) %>%
  select(Rank, Nomenclature, Signature, Elements, `Multi-omic feature`, `RCD count`, `RCD forms`) %>%  # Reorder columns
  arrange(desc(Rank))  # Sort in descending order by Rank

rio::export(selected_columns_combined_genotype_rows, "Table 3. Seven top-ranked signatures by multi-omic feature.xlsx")

#####
#####
#####
#####
# Select rows where "Genotype" is "transcript" and create df114
df114 <- df101 %>%
  filter(Genotype == "Transcript")

# Total number of rows in df71 (excluding NAs in TNC)
total_rows <- sum(!is.na(df114$TNC))

# Count the number of rows where "TNC" is "1"
count_1 <- sum(df114$TNC == "1", na.rm = TRUE)

# Count the number of rows where "TNC" is "2"
count_2 <- sum(df114$TNC == "2", na.rm = TRUE)

# Count the number of rows where "TNC" is "3"
count_3 <- sum(df114$TNC == "3", na.rm = TRUE)

# Calculate the percentages
percent_1 <- (count_1 / total_rows) * 100
percent_2 <- (count_2 / total_rows) * 100
percent_3 <- (count_3 / total_rows) * 100

# Display the counts and percentages
cat("Number of rows with '1' in the 'TNC' column:", count_1, "(", round(percent_1, 2), "%)\n")
cat("Number of rows with '2' in the 'TNC' column:", count_2, "(", round(percent_2, 2), "%)\n")
cat("Number of rows with '3' in the 'TNC' column:", count_3, "(", round(percent_3, 2), "%)\n")

# Calculate mean and median for the "Count_source" variable in df114
mean_value <- mean(df114$Count_source, na.rm = TRUE)
median_value <- median(df114$Count_source, na.rm = TRUE)

# Display the mean and median values
cat("Mean of Count_source:", mean_value, "\n")
cat("Median of Count_source:", median_value, "\n")

# Count the number of rows where "Count_source" is equal to 1
count_equal_1 <- sum(df114$Count_source == 1, na.rm = TRUE)

# Display the result
cat("Number of rows with Count_source equal to 1:", count_equal_1, "\n")

# Calculate quartiles (0%, 25%, 50%, 75%, 100%)
quantiles <- quantile(df114$Count_source, probs = seq(0, 1, 0.25), na.rm = TRUE)

# Display the quartiles
print(quantiles)

# Optional: Summary of basic statistics (includes min, 1st quartile, median, mean, 3rd quartile, and max)
summary_stats <- summary(df114$Count_source)
print(summary_stats)

# Calculate deciles (0%, 10%, 20%, ..., 100%)
deciles <- quantile(df114$Count_source, probs = seq(0, 1, 0.1), na.rm = TRUE)

# Display the deciles
print(deciles)

# Calculate percentiles (1%, 2%, ..., 100%)
percentiles <- quantile(df114$Count_source, probs = seq(0, 1, 0.01), na.rm = TRUE)

# Display the percentiles
print(percentiles)

# Select rows where all specified columns have the value "Risky" and create df114_COX_Risky
df114_COX_Risky <- df114 %>%
  filter(
    Type_Cox_DSS == "Risky" & 
      Type_Cox_DFI == "Risky" & 
      Type_Cox_PFI == "Risky" & 
      Type_Cox_OS == "Risky"
  )

# Calculate the total number of rows in df114
total_rows <- nrow(df114)

# Calculate the number of rows in df114_COX_Risky
risky_rows <- nrow(df114_COX_Risky)

# Calculate the percentage of "Risky" rows over the total rows
percentage_risky <- (risky_rows / total_rows) * 100

# Display the result
cat("Percentage of 'Risky' rows over total rows in df114:", round(percentage_risky, 2), "%\n")

# Display the first few rows of df114_COX_Risky to confirm the selection
head(df114_COX_Risky)


####################
##################### Select rows where all specified columns have the value "Risky" and create df114_COX_Risky
df114_COX_Protective <- df114 %>%
  filter(
    Type_Cox_DSS == "Protective" & 
      Type_Cox_DFI == "Protective" & 
      Type_Cox_PFI == "Protective" & 
      Type_Cox_OS == "Protective"
  )

# Calculate the total number of rows in df114
total_rows <- nrow(df114)

# Calculate the number of rows in df114_COX_Risky
protective_rows <- nrow(df114_COX_Protective)

# Calculate the percentage of "Risky" rows over the total rows
percentage_protective <- (protective_rows / total_rows) * 100

# Display the result
cat("Percentage of 'Protective' rows over total rows in df114:", round(percentage_protective, 2), "%\n")

unique(df114_COX_Risky$GFC)
#[1] 5
unique(df114_COX_Protective$GFC)
#[1] 5

# Estimate the distribution of rows by the unique values under "PFC"
pfc_distribution <- df114_COX_Protective %>%
  group_by(PFC) %>%
  summarise(
    Count = n(),  # Absolute count of each unique value in PFC
    Percentage = (Count / nrow(df114_COX_Protective)) * 100  # Percentage of each unique value in PFC
  )

# Display the distribution with absolute counts and percentages
print(pfc_distribution)

# Estimate the distribution of rows by the unique values under "PFC" for df114_COX_risky
pfc_distribution_risky <- df114_COX_Risky %>%
  group_by(PFC) %>%
  summarise(
    Count = n(),  # Absolute count of each unique value in PFC
    Percentage = (Count / nrow(df114_COX_Risky)) * 100  # Percentage of each unique value in PFC
  )

# Display the distribution with absolute counts and percentages
print(pfc_distribution_risky)

# Estimate the distribution of rows by the unique values under "TNC" for df114_COX_risky
tnc_distribution_risky <- df114_COX_Risky %>%
  group_by(TNC) %>%
  summarise(
    Count = n(),  # Absolute count of each unique value in TNC
    Percentage = (Count / nrow(df114_COX_Risky)) * 100  # Percentage of each unique value in TNC
  )

# Display the distribution with absolute counts and percentages
print(tnc_distribution_risky)

# Estimate the distribution of rows by the unique values under "TNC" for df114_COX_Protective
tnc_distribution_protective <- df114_COX_Protective %>%
  group_by(TNC) %>%
  summarise(
    Count = n(),  # Absolute count of each unique value in TNC
    Percentage = (Count / nrow(df114_COX_Protective)) * 100  # Percentage of each unique value in TNC
  )

# Display the distribution with absolute counts and percentages
print(tnc_distribution_protective)

# Estimate mean value and range (lower and upper values) for "Count_source" in df114_COX_risky
summary_risky <- df114_COX_Risky %>%
  summarise(
    Mean = mean(Count_source, na.rm = TRUE),  # Calculate mean value
    Lower = min(Count_source, na.rm = TRUE),  # Calculate lower (minimum) value
    Upper = max(Count_source, na.rm = TRUE)   # Calculate upper (maximum) value
  )

# Display the summary statistics
print(summary_risky)

# Estimate mean value and range (lower and upper values) for "Count_source" in df114_COX_Protective
summary_protective <- df114_COX_Protective %>%
  summarise(
    Mean = mean(Count_source, na.rm = TRUE),  # Calculate mean value
    Lower = min(Count_source, na.rm = TRUE),  # Calculate lower (minimum) value
    Upper = max(Count_source, na.rm = TRUE)   # Calculate upper (maximum) value
  )

# Display the summary statistics
print(summary_protective)

# Select rows where "Genotype" is "miRNA" and create df115
df115 <- df101 %>%
  filter(Genotype == "miRNA")

# Total number of rows in df71 (excluding NAs in TNC)
total_rows <- sum(!is.na(df115$TNC))

# Count the number of rows where "TNC" is "1"
count_1 <- sum(df115$TNC == "1", na.rm = TRUE)

# Count the number of rows where "TNC" is "2"
count_2 <- sum(df115$TNC == "2", na.rm = TRUE)

# Count the number of rows where "TNC" is "3"
count_3 <- sum(df115$TNC == "3", na.rm = TRUE)

# Calculate the percentages
percent_1 <- (count_1 / total_rows) * 100
percent_2 <- (count_2 / total_rows) * 100
percent_3 <- (count_3 / total_rows) * 100

# Display the counts and percentages
cat("Number of rows with '1' in the 'TNC' column:", count_1, "(", round(percent_1, 2), "%)\n")
cat("Number of rows with '2' in the 'TNC' column:", count_2, "(", round(percent_2, 2), "%)\n")
cat("Number of rows with '3' in the 'TNC' column:", count_3, "(", round(percent_3, 2), "%)\n")

# Calculate mean and median for the "Count_source" variable in df115
mean_value <- mean(df115$Count_source, na.rm = TRUE)
median_value <- median(df115$Count_source, na.rm = TRUE)

# Display the mean and median values
cat("Mean of Count_source:", mean_value, "\n")
cat("Median of Count_source:", median_value, "\n")

# Count the number of rows where "Count_source" is equal to 1
count_equal_1 <- sum(df115$Count_source == 1, na.rm = TRUE)

# Display the result
cat("Number of rows with Count_source equal to 1:", count_equal_1, "\n")

# Overall Risky and Protective miRNA signatures
340*100/1470
#[1] 23.12925
448*100/1470
#[1] 30.47619

  # Calculate quartiles (0%, 25%, 50%, 75%, 100%)
  quantiles <- quantile(df115$Count_source, probs = seq(0, 1, 0.25), na.rm = TRUE)

# Display the quartiles
print(quantiles)

# Optional: Summary of basic statistics (includes min, 1st quartile, median, mean, 3rd quartile, and max)
summary_stats <- summary(df115$Count_source)
print(summary_stats)

# Calculate deciles (0%, 10%, 20%, ..., 100%)
deciles <- quantile(df115$Count_source, probs = seq(0, 1, 0.1), na.rm = TRUE)

# Display the deciles
print(deciles)

# Calculate percentiles (1%, 2%, ..., 100%)
percentiles <- quantile(df115$Count_source, probs = seq(0, 1, 0.01), na.rm = TRUE)

# Display the percentiles
print(percentiles)

# Select rows where all specified columns have the value "Risky" and create df115_COX_Risky
df115_COX_Risky <- df115 %>%
  filter(
    Type_Cox_DSS == "Risky" & 
      Type_Cox_DFI == "Risky" & 
      Type_Cox_PFI == "Risky" & 
      Type_Cox_OS == "Risky"
  )

# Calculate the total number of rows in df115
total_rows <- nrow(df115)

# Calculate the number of rows in df115_COX_Risky
risky_rows <- nrow(df115_COX_Risky)

# Calculate the percentage of "Risky" rows over the total rows
percentage_risky <- (risky_rows / total_rows) * 100

# Display the result
cat("Percentage of 'Risky' rows over total rows in df115:", round(percentage_risky, 2), "%\n")

# Display the first few rows of df115_COX_Risky to confirm the selection
head(df115_COX_Risky)

####################
##################### Select rows where all specified columns have the value "Risky" and create df115_COX_Risky
df115_COX_Protective <- df115 %>%
  filter(
    Type_Cox_DSS == "Protective" & 
      Type_Cox_DFI == "Protective" & 
      Type_Cox_PFI == "Protective" & 
      Type_Cox_OS == "Protective"
  )

# Calculate the total number of rows in df115
total_rows <- nrow(df115)

# Calculate the number of rows in df115_COX_Risky
protective_rows <- nrow(df115_COX_Protective)

# Calculate the percentage of "Risky" rows over the total rows
percentage_protective <- (protective_rows / total_rows) * 100

# Display the result
cat("Percentage of 'Protective' rows over total rows in df115:", round(percentage_protective, 2), "%\n")

unique(df115_COX_Risky$GFC)
#[1] 5
unique(df115_COX_Protective$GFC)
#[1] 5

# Estimate the distribution of rows by the unique values under "PFC"
pfc_distribution <- df115_COX_Protective %>%
  group_by(PFC) %>%
  summarise(
    Count = n(),  # Absolute count of each unique value in PFC
    Percentage = (Count / nrow(df115_COX_Protective)) * 100  # Percentage of each unique value in PFC
  )

# Display the distribution with absolute counts and percentages
print(pfc_distribution)

# Estimate the distribution of rows by the unique values under "PFC" for df115_COX_risky
pfc_distribution_risky <- df115_COX_Risky %>%
  group_by(PFC) %>%
  summarise(
    Count = n(),  # Absolute count of each unique value in PFC
    Percentage = (Count / nrow(df115_COX_Risky)) * 100  # Percentage of each unique value in PFC
  )

# Display the distribution with absolute counts and percentages
print(pfc_distribution_risky)

# Estimate the distribution of rows by the unique values under "TNC" for df115_COX_risky
tnc_distribution_risky <- df115_COX_Risky %>%
  group_by(TNC) %>%
  summarise(
    Count = n(),  # Absolute count of each unique value in TNC
    Percentage = (Count / nrow(df115_COX_Risky)) * 100  # Percentage of each unique value in TNC
  )

# Display the distribution with absolute counts and percentages
print(tnc_distribution_risky)

# Estimate the distribution of rows by the unique values under "TNC" for df115_COX_Protective
tnc_distribution_protective <- df115_COX_Protective %>%
  group_by(TNC) %>%
  summarise(
    Count = n(),  # Absolute count of each unique value in TNC
    Percentage = (Count / nrow(df115_COX_Protective)) * 100  # Percentage of each unique value in TNC
  )

# Display the distribution with absolute counts and percentages
print(tnc_distribution_protective)

# Estimate mean value and range (lower and upper values) for "Count_source" in df115_COX_risky
summary_risky <- df115_COX_Risky %>%
  summarise(
    Mean = mean(Count_source, na.rm = TRUE),  # Calculate mean value
    Lower = min(Count_source, na.rm = TRUE),  # Calculate lower (minimum) value
    Upper = max(Count_source, na.rm = TRUE)   # Calculate upper (maximum) value
  )

# Display the summary statistics
print(summary_risky)

# Estimate mean value and range (lower and upper values) for "Count_source" in df115_COX_Protective
summary_protective <- df115_COX_Protective %>%
  summarise(
    Mean = mean(Count_source, na.rm = TRUE),  # Calculate mean value
    Lower = min(Count_source, na.rm = TRUE),  # Calculate lower (minimum) value
    Upper = max(Count_source, na.rm = TRUE)   # Calculate upper (maximum) value
  )

# Display the summary statistics
print(summary_protective)

# Select rows where "Genotype" is "transcript" and create df116
df116 <- df101 %>%
  filter(Genotype == "Methylation")

# Total number of rows in df71 (excluding NAs in PFC)
total_rows <- sum(!is.na(df116$PFC))

# Count the number of rows where "PFC" is "1"
count_1 <- sum(df116$PFC == "1", na.rm = TRUE)

# Count the number of rows where "PFC" is "2"
count_2 <- sum(df116$PFC == "2", na.rm = TRUE)

# Count the number of rows where "PFC" is "3"
count_3 <- sum(df116$PFC == "3", na.rm = TRUE)

# Calculate the percentages
percent_1 <- (count_1 / total_rows) * 100
percent_2 <- (count_2 / total_rows) * 100
percent_3 <- (count_3 / total_rows) * 100

# Display the counts and percentages
cat("Number of rows with '1' in the 'PFC' column:", count_1, "(", round(percent_1, 2), "%)\n")
cat("Number of rows with '2' in the 'PFC' column:", count_2, "(", round(percent_2, 2), "%)\n")
cat("Number of rows with '3' in the 'PFC' column:", count_3, "(", round(percent_3, 2), "%)\n")

# Create the original df117 by selecting rows where "Genotype" is "Methylation" and "PFC" is "3"
df117 <- df101 %>%
  filter(Genotype == "Methylation", PFC == "3")

# Total row count of df117 for calculating percentages
total_count <- nrow(df117)

# Select rows where "Risky" appears in any column and calculate row count and percentage
df117_Risky <- df117 %>%
  filter(if_any(everything(), ~ . == "Risky"))

risky_count <- nrow(df117_Risky)
risky_percentage <- (risky_count / total_count) * 100
cat("Rows with 'Risky':", risky_count, "(", round(risky_percentage, 2), "%)\n")

# Select rows where both "Protective" and "Hot" appear in any column
df117_Protective_hot <- df117 %>%
  filter(if_any(everything(), ~ . == "Protective") & if_any(everything(), ~ . == "Hot"))

protective_hot_count <- nrow(df117_Protective_hot)
protective_hot_percentage <- (protective_hot_count / total_count) * 100
cat("Rows with 'Protective' and 'Hot':", protective_hot_count, "(", round(protective_hot_percentage, 2), "%)\n")

# Select rows where both "Protective" and "Cold" appear in any column
df117_Protective_cold <- df117 %>%
  filter(if_any(everything(), ~ . == "Protective") & if_any(everything(), ~ . == "Cold"))

protective_cold_count <- nrow(df117_Protective_cold)
protective_cold_percentage <- (protective_cold_count / total_count) * 100
cat("Rows with 'Protective' and 'Cold':", protective_cold_count, "(", round(protective_cold_percentage, 2), "%)\n")

# Select rows where both "Risky" and "Hot" appear in any column
df117_Risky_hot <- df117 %>%
  filter(if_any(everything(), ~ . == "Risky") & if_any(everything(), ~ . == "Hot"))

risky_hot_count <- nrow(df117_Risky_hot)
risky_hot_percentage <- (risky_hot_count / total_count) * 100
cat("Rows with 'Risky' and 'Hot':", risky_hot_count, "(", round(risky_hot_percentage, 2), "%)\n")

# Select rows where both "Risky" and "Cold" appear in any column
df117_Risky_cold <- df117 %>%
  filter(if_any(everything(), ~ . == "Risky") & if_any(everything(), ~ . == "Cold"))

risky_cold_count <- nrow(df117_Risky_cold)
risky_cold_percentage <- (risky_cold_count / total_count) * 100
cat("Rows with 'Risky' and 'Cold':", risky_cold_count, "(", round(risky_cold_percentage, 2), "%)\n")

# Calculate row counts and percentages for other categories
df117_hot <- df117 %>%
  filter(if_any(everything(), ~ . == "Hot"))

hot_count <- nrow(df117_hot)
hot_percentage <- (hot_count / total_count) * 100
cat("Rows with 'Hot':", hot_count, "(", round(hot_percentage, 2), "%)\n")

df117_cold <- df117 %>%
  filter(if_any(everything(), ~ . == "Cold"))

cold_count <- nrow(df117_cold)
cold_percentage <- (cold_count / total_count) * 100
cat("Rows with 'Cold':", cold_count, "(", round(cold_percentage, 2), "%)\n")

df117_variable <- df117 %>%
  filter(if_any(everything(), ~ . == "Variable"))

variable_count <- nrow(df117_variable)
variable_percentage <- (variable_count / total_count) * 100
cat("Rows with 'Variable':", variable_count, "(", round(variable_percentage, 2), "%)\n")

df117_antitumoral <- df117 %>%
  filter(if_any(everything(), ~ . == "anti-tumoral"))

antitumoral_count <- nrow(df117_antitumoral)
antitumoral_percentage <- (antitumoral_count / total_count) * 100
cat("Rows with 'Anti-tumoral':", antitumoral_count, "(", round(antitumoral_percentage, 2), "%)\n")

df117_protumoral <- df117 %>%
  filter(if_any(everything(), ~ . == "pro-tumoral"))

protumoral_count <- nrow(df117_protumoral)
protumoral_percentage <- (protumoral_count / total_count) * 100
cat("Rows with 'Pro-tumoral':", protumoral_count, "(", round(protumoral_percentage, 2), "%)\n")

df117_dual <- df117 %>%
  filter(if_any(everything(), ~ . == "dual"))

dual_count <- nrow(df117_dual)
dual_percentage <- (dual_count / total_count) * 100
cat("Rows with 'Dual':", dual_count, "(", round(dual_percentage, 2), "%)\n")

df117_hot <- df117 %>%
  filter(if_any(everything(), ~ . == "Hot"))

hot_count <- nrow(df117_hot)
hot_percentage <- (hot_count / total_count) * 100
cat("Rows with 'Hot':", hot_count, "(", round(hot_percentage, 2), "%)\n")

df117_cold <- df117 %>%
  filter(if_any(everything(), ~ . == "Cold"))

cold_count <- nrow(df117_cold)
cold_percentage <- (cold_count / total_count) * 100
cat("Rows with 'Cold':", cold_count, "(", round(cold_percentage, 2), "%)\n")

# Calculate mean and median for the "Count_source" variable in df116
mean_value <- mean(df116$Count_source, na.rm = TRUE)
median_value <- median(df116$Count_source, na.rm = TRUE)

# Display the mean and median values
cat("Mean of Count_source:", mean_value, "\n")
cat("Median of Count_source:", median_value, "\n")

# Count the number of rows where "Count_source" is equal to 1
count_equal_1 <- sum(df116$Count_source == 1, na.rm = TRUE)

# Display the result
cat("Number of rows with Count_source equal to 1:", count_equal_1, "\n")

#Overall Risky and Protetive Methylation signatures
340*100/1470
#[1] 23.12925
448*100/1470
#30.47619

# Calculate quartiles (0%, 25%, 50%, 75%, 100%)
quantiles <- quantile(df116$Count_source, probs = seq(0, 1, 0.25), na.rm = TRUE)

# Display the quartiles
print(quantiles)

# Optional: Summary of basic statistics (includes min, 1st quartile, median, mean, 3rd quartile, and max)
summary_stats <- summary(df116$Count_source)
print(summary_stats)

# Calculate deciles (0%, 10%, 20%, ..., 100%)
deciles <- quantile(df116$Count_source, probs = seq(0, 1, 0.1), na.rm = TRUE)

# Display the deciles
print(deciles)

# Calculate percentiles (1%, 2%, ..., 100%)
percentiles <- quantile(df116$Count_source, probs = seq(0, 1, 0.01), na.rm = TRUE)

# Display the percentiles
print(percentiles)

# Select rows where all specified columns have the value "Risky" and create df116_COX_Risky
df116_COX_Risky <- df116 %>%
  filter(
    Type_Cox_DSS == "Risky" & 
      Type_Cox_DFI == "Risky" & 
      Type_Cox_PFI == "Risky" & 
      Type_Cox_OS == "Risky"
  )

# Calculate the total number of rows in df116
total_rows <- nrow(df116)

# Calculate the number of rows in df116_COX_Risky
risky_rows <- nrow(df116_COX_Risky)

# Calculate the percentage of "Risky" rows over the total rows
percentage_risky <- (risky_rows / total_rows) * 100

# Display the result
cat("Percentage of 'Risky' rows over total rows in df116:", round(percentage_risky, 2), "%\n")

# Display the first few rows of df116_COX_Risky to confirm the selection
head(df116_COX_Risky)

####################
##################### Select rows where all specified columns have the value "Risky" and create df116_COX_Risky
df116_COX_Protective <- df116 %>%
  filter(
    Type_Cox_DSS == "Protective" & 
      Type_Cox_DFI == "Protective" & 
      Type_Cox_PFI == "Protective" & 
      Type_Cox_OS == "Protective"
  )

# Calculate the total number of rows in df116
total_rows <- nrow(df116)

# Calculate the number of rows in df116_COX_Risky
protective_rows <- nrow(df116_COX_Protective)

# Calculate the percentage of "Risky" rows over the total rows
percentage_protective <- (protective_rows / total_rows) * 100

# Display the result
cat("Percentage of 'Protective' rows over total rows in df116:", round(percentage_protective, 2), "%\n")

unique(df116_COX_Risky$GFC)
#[1] 5
unique(df116_COX_Protective$GFC)
#[1] 5

# Estimate the distribution of rows by the unique values under "PFC"
pfc_distribution <- df116_COX_Protective %>%
  group_by(PFC) %>%
  summarise(
    Count = n(),  # Absolute count of each unique value in PFC
    Percentage = (Count / nrow(df116_COX_Protective)) * 100  # Percentage of each unique value in PFC
  )

# Display the distribution with absolute counts and percentages
print(pfc_distribution)

# Estimate the distribution of rows by the unique values under "PFC" for df116_COX_risky
pfc_distribution_risky <- df116_COX_Risky %>%
  group_by(PFC) %>%
  summarise(
    Count = n(),  # Absolute count of each unique value in PFC
    Percentage = (Count / nrow(df116_COX_Risky)) * 100  # Percentage of each unique value in PFC
  )

# Display the distribution with absolute counts and percentages
print(pfc_distribution_risky)

# Estimate the distribution of rows by the unique values under "PFC" for df116_COX_risky
PFC_distribution_risky <- df116_COX_Risky %>%
  group_by(PFC) %>%
  summarise(
    Count = n(),  # Absolute count of each unique value in PFC
    Percentage = (Count / nrow(df116_COX_Risky)) * 100  # Percentage of each unique value in PFC
  )

# Display the distribution with absolute counts and percentages
print(PFC_distribution_risky)

# Estimate the distribution of rows by the unique values under "PFC" for df116_COX_Protective
PFC_distribution_protective <- df116_COX_Protective %>%
  group_by(PFC) %>%
  summarise(
    Count = n(),  # Absolute count of each unique value in PFC
    Percentage = (Count / nrow(df116_COX_Protective)) * 100  # Percentage of each unique value in PFC
  )

# Display the distribution with absolute counts and percentages
print(PFC_distribution_protective)

# Estimate mean value and range (lower and upper values) for "Count_source" in df116_COX_risky
summary_risky <- df116_COX_Risky %>%
  summarise(
    Mean = mean(Count_source, na.rm = TRUE),  # Calculate mean value
    Lower = min(Count_source, na.rm = TRUE),  # Calculate lower (minimum) value
    Upper = max(Count_source, na.rm = TRUE)   # Calculate upper (maximum) value
  )

# Display the summary statistics
print(summary_risky)

# Estimate mean value and range (lower and upper values) for "Count_source" in df116_COX_Protective
summary_protective <- df116_COX_Protective %>%
  summarise(
    Mean = mean(Count_source, na.rm = TRUE),  # Calculate mean value
    Lower = min(Count_source, na.rm = TRUE),  # Calculate lower (minimum) value
    Upper = max(Count_source, na.rm = TRUE)   # Calculate upper (maximum) value
  )

# Display the summary statistics
print(summary_protective)

dfBia <- complete_top_signatures

# Seleciona todas as linhas de dfBia usando as colunas "Signature", "CTAB", e "RCD_types"
selected_rows_Bia <- dfBia[, c("Signature", "CTAB", "RCD_types")]

dfBia_2 <-  selected_rows_Bia 
rio::export(dfBia_2, "Bia_top_signatures.xlsx")

# Exibir o resultado
print(selected_rows_Bia)

# Extract the "RCD_types" column
rcd_values <- dfBia[["RCD_types"]]

# Split the values by "/", keeping solo values as they are
split_rcd_values <- lapply(rcd_values, function(x) unlist(strsplit(x, split = "/")))

# Create a dataframe where each row contains the individual value, keeping the original structure
df_split_rcd <- data.frame(RCD_types = unlist(split_rcd_values))

# Remove duplicate values based on unique rows
df_unique_rcd <- unique(df_split_rcd)

rio::export(df_unique_rcd, "RCD_top_signatures.xlsx")

# Exibir o resultado
print(df_unique_rcd)

dfBia_worst <- complete_worst_signatures

# Seleciona todas as linhas de dfBia_worst usando as colunas "Signature", "CTAB", e "RCD_types"
selected_rows_Bia_worst <- dfBia_worst[, c("Signature", "CTAB", "RCD_types")]

dfBia_worst_2 <-  selected_rows_Bia_worst 
rio::export(dfBia_worst_2, "Bia_worst_signatures.xlsx")

# Extract the "RCD_types" column
rcd_values_worst <- dfBia_worst[["RCD_types"]]

# Split the values by "/", keeping solo values as they are
split_rcd_values_worst <- lapply(rcd_values_worst, function(x) unlist(strsplit(x, split = "/")))

# Create a dataframe where each row contains the individual value, keeping the original structure
df_split_rcd_worst <- data.frame(RCD_types = unlist(split_rcd_values_worst))

# Remove duplicate values based on unique rows
df_unique_rcd_worst <- unique(df_split_rcd_worst)

rio::export(df_unique_rcd_worst, "RCD_worst_signatures.xlsx")

# Filtrar para manter apenas linhas onde a coluna genotipica seja "mRNA" e "Transcript" (Emanuell)
rank <-  df101

rank_filtered <- rank %>%
  filter(
    (Genotype == "mRNA" | Genotype == "Transcript"), 
    TNC != "1",
    Type_Cox_OS != "NS", Type_Cox_DSS != "NS", Type_Cox_DFI != "NS", Type_Cox_PFI != "NS", # Pelo menos uma métrica de Cox diferente de "NS"
    Type_log_rank_OS != "NS", Type_log_rank_DSS != "NS", Type_log_rank_DFI != "NS", Type_log_rank_PFI != "NS", # Pelo menos uma métrica de log-rank diferente de "NS"
    microenvironment_classification != "NS",
    immune_classification != "NS"
  )

rio::export(rank_filtered, "Emanuell_Everest_top_signatures.xlsx")

######
###### Distribution of protein-specific signatures
###### 
# Create df118 by selecting rows where "Genotype" is "Protein"
df118 <- df101 %>%
  filter(Genotype == "Protein")

# Create df118_PFC_3 by selecting rows where "Genotype" is "Protein" and "PFC" is "3"
df118_PFC_3 <- df101 %>%
  filter(Genotype == "Protein", PFC == "3")

# Total row count of df118_PFC_3
total_count_PFC_3 <- nrow(df118_PFC_3)

# Total row count of df118
total_count_df118 <- nrow(df118)

# Calculate the percentage of rows in df118_PFC_3 relative to df118
percentage_PFC_3 <- (total_count_PFC_3 / total_count_df118) * 100

# Print the percentage
print(paste("Percentage of rows in df118_PFC_3 relative to df118:", round(percentage_PFC_3, 2), "%"))

# Count rows in df118_PFC_3 where "SCS" is "P"
df118_PFC_3_SCS_P <- df118_PFC_3 %>%
  filter(SCS == "P")

total_count_SCS_P <- nrow(df118_PFC_3_SCS_P)

# Calculate the percentage of rows with "SCS" as "P" relative to df118_PFC_3
percentage_SCS_P <- (total_count_SCS_P / total_count_PFC_3) * 100

# Print the percentage for "SCS" == "P" in df118_PFC_3
print(paste("Percentage of rows with SCS as 'P' in df118_PFC_3:", round(percentage_SCS_P, 2), "%"))

# Count rows in df118_PFC_3 where "SCS" is "N"
df118_PFC_3_SCS_N <- df118_PFC_3 %>%
  filter(SCS == "N")

total_count_SCS_N <- nrow(df118_PFC_3_SCS_N)

# Calculate the percentage of rows with "SCS" as "N" relative to df118_PFC_3
percentage_SCS_N <- (total_count_SCS_N / total_count_PFC_3) * 100

# Print the percentage for "SCS" == "N" in df118_PFC_3
print(paste("Percentage of rows with SCS as 'N' in df118_PFC_3:", round(percentage_SCS_N, 2), "%"))

# Select rows where "Risky" appears in any column and calculate row count and percentage
df118_Risky <- df118 %>%
  filter(if_any(everything(), ~ . == "Risky"))

risky_count <- nrow(df118_Risky)
risky_percentage <- (risky_count / total_count_df118) * 100
cat("Rows with 'Risky':", risky_count, "(", round(risky_percentage, 2), "%)\n")

# Select rows where both "Protective" and "Hot" appear in any column
df118_Protective <- df118 %>%
  filter(if_any(everything(), ~ . == "Protective")) #& if_any(everything(), ~ . == "Hot"))

protective_count <- nrow(df118_Protective)
protective_percentage <- (protective_count / total_count_df118) * 100
cat("Rows with 'Protective':", protective_count, "(", round(protective_percentage, 2), "%)\n")

# Select rows where both "Protective" and "Cold" appear in any column
df118_cold <- df118 %>%
  filter(if_any(everything(), ~ . == "Cold"))

df118_hot <- df118 %>%
  filter(if_any(everything(), ~ . == "Hot"))

cold_count <- nrow(df118_cold)
cold_percentage <- (cold_count / total_count_df118) * 100
cat("Rows with 'Cold':", cold_count, "(", round(cold_percentage, 2), "%)\n")

hot_count <- nrow(df118_hot)
hot_percentage <- (hot_count / total_count_df118) * 100
cat("Rows with 'Hot':", hot_count, "(", round(hot_percentage, 2), "%)\n")

df118_variable <- df118 %>%
  filter(if_any(everything(), ~ . == "Variable"))
variable_count <- nrow(df118_variable)
variable_percentage <- (variable_count / total_count_df118) * 100
cat("Rows with 'Variable':", variable_count, "(", round(variable_percentage, 2), "%)\n")

df118_pro_tumoral <- df118 %>%
  filter(if_any(everything(), ~ . == "pro-tumoral"))
pro_tumoral_count <- nrow(df118_pro_tumoral)
pro_tumoral_percentage <- (pro_tumoral_count / total_count_df118) * 100
cat("Rows with 'pro_tumoral':", pro_tumoral_count, "(", round(pro_tumoral_percentage, 2), "%)\n")

df118_dual <- df118 %>%
  filter(if_any(everything(), ~ . == "dual"))
dual_count <- nrow(df118_dual)
dual_percentage <- (dual_count / total_count_df118) * 100
cat("Rows with 'dual':", dual_count, "(", round(dual_percentage, 2), "%)\n")

# Select rows where both "Risky" appear in any column
df118_Risky_hot <- df118 %>%
  filter(if_any(everything(), ~ . == "Risky") & if_any(everything(), ~ . == "Hot"))

risky_hot_count <- nrow(df118_Risky_hot)
risky_hot_percentage <- (risky_hot_count / total_count) * 100
cat("Rows with 'Risky' and 'Hot':", risky_hot_count, "(", round(risky_hot_percentage, 2), "%)\n")

risky_hot_count <- nrow(df118_Risky_hot)
risky_hot_percentage <- (risky_hot_count / total_count) * 100
cat("Rows with 'Risky' and 'Hot':", risky_hot_count, "(", round(risky_hot_percentage, 2), "%)\n")

# Select rows where both "Risky" and "Cold" appear in any column
df118_Risky_cold <- df118 %>%
  filter(if_any(everything(), ~ . == "Risky") & if_any(everything(), ~ . == "Cold"))

risky_cold_count <- nrow(df118_Risky_cold)
risky_cold_percentage <- (risky_cold_count / total_count) * 100
cat("Rows with 'Risky' and 'Cold':", risky_cold_count, "(", round(risky_cold_percentage, 2), "%)\n")

# Calculate row counts and percentages for other categories
df118_hot <- df118 %>%
  filter(if_any(everything(), ~ . == "Hot"))

hot_count <- nrow(df118_hot)
hot_percentage <- (hot_count / total_count) * 100
cat("Rows with 'Hot':", hot_count, "(", round(hot_percentage, 2), "%)\n")

df118_cold <- df118 %>%
  filter(if_any(everything(), ~ . == "Cold"))

cold_count <- nrow(df118_cold)
cold_percentage <- (cold_count / total_count_df118) * 100
cat("Rows with 'Cold':", cold_count, "(", round(cold_percentage, 2), "%)\n")

df118_variable <- df118 %>%
  filter(if_any(everything(), ~ . == "Variable"))

variable_count <- nrow(df118_variable)
variable_percentage <- (variable_count / total_count) * 100
cat("Rows with 'Variable':", variable_count, "(", round(variable_percentage, 2), "%)\n")

df118_antitumoral <- df118 %>%
  filter(if_any(everything(), ~ . == "anti-tumoral"))

antitumoral_count <- nrow(df118_antitumoral)
antitumoral_percentage <- (antitumoral_count / total_count_df118) * 100
cat("Rows with 'Anti-tumoral':", antitumoral_count, "(", round(antitumoral_percentage, 2), "%)\n")

df118_protumoral <- df118 %>%
  filter(if_any(everything(), ~ . == "pro-tumoral"))

protumoral_count <- nrow(df118_protumoral)
protumoral_percentage <- (protumoral_count / total_count) * 100
cat("Rows with 'Pro-tumoral':", protumoral_count, "(", round(protumoral_percentage, 2), "%)\n")

df118_dual <- df118 %>%
  filter(if_any(everything(), ~ . == "dual"))

dual_count <- nrow(df118_dual)
dual_percentage <- (dual_count / total_count) * 100
cat("Rows with 'Dual':", dual_count, "(", round(dual_percentage, 2), "%)\n")

df118_hot <- df118 %>%
  filter(if_any(everything(), ~ . == "Hot"))

hot_count <- nrow(df118_hot)
hot_percentage <- (hot_count / total_count) * 100
cat("Rows with 'Hot':", hot_count, "(", round(hot_percentage, 2), "%)\n")

df118_cold <- df118 %>%
  filter(if_any(everything(), ~ . == "Cold"))

cold_count <- nrow(df118_cold)
cold_percentage <- (cold_count / total_count) * 100
cat("Rows with 'Cold':", cold_count, "(", round(cold_percentage, 2), "%)\n")

# Calculate mean and median for the "Count_source" variable in df118
mean_value <- mean(df118$Count_source, na.rm = TRUE)
median_value <- median(df118$Count_source, na.rm = TRUE)

# Display the mean and median values
cat("Mean of Count_source:", mean_value, "\n")
cat("Median of Count_source:", median_value, "\n")

######
###### Distribution of Mutation-specific signatures
###### 

# Create df119 and its subsets based on "Genotype" and "PFC" values
df119 <- df101 %>% filter(Genotype == "Mutation")
df119_PFC_1 <- df119 %>% filter(PFC == "1")
df119_PFC_2 <- df119 %>% filter(PFC == "2")
df119_PFC_3 <- df119 %>% filter(PFC == "3")

# Count the total rows for each subset
total_count_PFC_1 <- nrow(df119_PFC_1)
total_count_PFC_2 <- nrow(df119_PFC_2)
total_count_PFC_3 <- nrow(df119_PFC_3)

# Total row count of df119
total_count_df119 <- nrow(df119)

# Calculate the percentage of rows in each subset relative to df119
percentage_PFC_1 <- (total_count_PFC_1 / total_count_df119) * 100
percentage_PFC_2 <- (total_count_PFC_2 / total_count_df119) * 100
percentage_PFC_3 <- (total_count_PFC_3 / total_count_df119) * 100

# Display the counts and percentages in the console
cat("Total count of PFC = 1:", total_count_PFC_1, 
    "- Percentage:", round(percentage_PFC_1, 2), "%\n")
cat("Total count of PFC = 2:", total_count_PFC_2, 
    "- Percentage:", round(percentage_PFC_2, 2), "%\n")
cat("Total count of PFC = 3:", total_count_PFC_3, 
    "- Percentage:", round(percentage_PFC_3, 2), "%\n")

# Count rows in df119_PFC_3 where "SCS" is "P"
df119_PFC_1_SCS_P <- df119_PFC_1 %>%
  filter(SCS == "P")

# Count rows in df119_PFC_3 where "SCS" is "P"
df119_PFC_2_SCS_P <- df119_PFC_2 %>%
  filter(SCS == "P")

# Count rows in df119_PFC_3 where "SCS" is "P"
df119_PFC_3_SCS_P <- df119_PFC_3 %>%
  filter(SCS == "P")
# Count rows in df119_PFC_3 where "SCS" is "P"
df119_PFC_3_SCS_N <- df119_PFC_3 %>%
  filter(SCS == "N")

# Select rows where "Risky" appears in any column and calculate row count and percentage
df119_Risky <- df119 %>%
  filter(if_any(everything(), ~ . == "Risky"))

risky_count <- nrow(df119_Risky)
risky_percentage <- (risky_count / total_count_df119) * 100
cat("Rows with 'Risky':", risky_count, "(", round(risky_percentage, 2), "%)\n")

# Select rows where both "Protective" appear in any column
df119_Protective <- df119 %>%
  filter(if_any(everything(), ~ . == "Protective")) #& if_any(everything(), ~ . == "Hot"))

protective_count <- nrow(df119_Protective)
protective_percentage <- (protective_count / total_count_df119) * 100
cat("Rows with 'Protective':", protective_count, "(", round(protective_percentage, 2), "%)\n")

# Select rows where both "Protective" and "Cold" appear in any column
df119_cold <- df119 %>%
  filter(if_any(everything(), ~ . == "Cold"))

df119_hot <- df119 %>%
  filter(if_any(everything(), ~ . == "Hot"))

cold_count <- nrow(df119_cold)
cold_percentage <- (cold_count / total_count_df119) * 100
cat("Rows with 'Cold':", cold_count, "(", round(cold_percentage, 2), "%)\n")

hot_count <- nrow(df119_hot)
hot_percentage <- (hot_count / total_count_df119) * 100
cat("Rows with 'Hot':", hot_count, "(", round(hot_percentage, 2), "%)\n")

df119_variable <- df119 %>%
  filter(if_any(everything(), ~ . == "Variable"))
variable_count <- nrow(df119_variable)
variable_percentage <- (variable_count / total_count_df119) * 100
cat("Rows with 'Variable':", variable_count, "(", round(variable_percentage, 2), "%)\n")

df119_pro_tumoral <- df119 %>%
  filter(if_any(everything(), ~ . == "pro-tumoral"))
pro_tumoral_count <- nrow(df119_pro_tumoral)
pro_tumoral_percentage <- (pro_tumoral_count / total_count_df119) * 100
cat("Rows with 'pro_tumoral':", pro_tumoral_count, "(", round(pro_tumoral_percentage, 2), "%)\n")

df119_dual <- df119 %>%
  filter(if_any(everything(), ~ . == "dual"))
dual_count <- nrow(df119_dual)
dual_percentage <- (dual_count / total_count_df119) * 100
cat("Rows with 'dual':", dual_count, "(", round(dual_percentage, 2), "%)\n")

# Count the rows where microenvironment_classification is "anti-tumoral"
total_anti_tumoral <- sum(df119$microenvironment_classification == "anti-tumoral")
percentage_anti_tumoral <- (total_anti_tumoral / total_count_df119) * 100

df119_anti_tumoral <- df119 %>%
  filter(if_any(everything(), ~ . == "anti-tumoral"))
anti_tumoral_count <- nrow(df119_anti_tumoral)
anti_tumoral_percentage <- (anti_tumoral_count / total_count_df119) * 100
cat("Rows with 'anti_tumoral':", anti_tumoral_count, "(", round(anti_tumoral_percentage, 2), "%)\n")

#####################
#####################
#####################
#####################
# Select rows where both "Risky" appear in any column
df119_Risky_hot <- df119 %>%
  filter(if_any(everything(), ~ . == "Risky") & if_any(everything(), ~ . == "Hot"))

risky_hot_count <- nrow(df119_Risky_hot)
risky_hot_percentage <- (risky_hot_count / total_count) * 100
cat("Rows with 'Risky' and 'Hot':", risky_hot_count, "(", round(risky_hot_percentage, 2), "%)\n")

risky_hot_count <- nrow(df119_Risky_hot)
risky_hot_percentage <- (risky_hot_count / total_count) * 100
cat("Rows with 'Risky' and 'Hot':", risky_hot_count, "(", round(risky_hot_percentage, 2), "%)\n")

# Select rows where both "Risky" and "Cold" appear in any column
df119_Risky_cold <- df119 %>%
  filter(if_any(everything(), ~ . == "Risky") & if_any(everything(), ~ . == "Cold"))

risky_cold_count <- nrow(df119_Risky_cold)
risky_cold_percentage <- (risky_cold_count / total_count) * 100
cat("Rows with 'Risky' and 'Cold':", risky_cold_count, "(", round(risky_cold_percentage, 2), "%)\n")

# Calculate row counts and percentages for other categories
df119_hot <- df119 %>%
  filter(if_any(everything(), ~ . == "Hot"))

hot_count <- nrow(df119_hot)
hot_percentage <- (hot_count / total_count) * 100
cat("Rows with 'Hot':", hot_count, "(", round(hot_percentage, 2), "%)\n")

df119_cold <- df119 %>%
  filter(if_any(everything(), ~ . == "Cold"))

cold_count <- nrow(df119_cold)
cold_percentage <- (cold_count / total_count_df119) * 100
cat("Rows with 'Cold':", cold_count, "(", round(cold_percentage, 2), "%)\n")

df119_variable <- df119 %>%
  filter(if_any(everything(), ~ . == "Variable"))

variable_count <- nrow(df119_variable)
variable_percentage <- (variable_count / total_count) * 100
cat("Rows with 'Variable':", variable_count, "(", round(variable_percentage, 2), "%)\n")

df119_antitumoral <- df119 %>%
  filter(if_any(everything(), ~ . == "anti-tumoral"))

antitumoral_count <- nrow(df119_antitumoral)
antitumoral_percentage <- (antitumoral_count / total_count_df119) * 100
cat("Rows with 'Anti-tumoral':", antitumoral_count, "(", round(antitumoral_percentage, 2), "%)\n")

df119_protumoral <- df119 %>%
  filter(if_any(everything(), ~ . == "pro-tumoral"))

protumoral_count <- nrow(df119_protumoral)
protumoral_percentage <- (protumoral_count / total_count) * 100
cat("Rows with 'Pro-tumoral':", protumoral_count, "(", round(protumoral_percentage, 2), "%)\n")

df119_dual <- df119 %>%
  filter(if_any(everything(), ~ . == "dual"))

dual_count <- nrow(df119_dual)
dual_percentage <- (dual_count / total_count) * 100
cat("Rows with 'Dual':", dual_count, "(", round(dual_percentage, 2), "%)\n")

df119_hot <- df119 %>%
  filter(if_any(everything(), ~ . == "Hot"))

hot_count <- nrow(df119_hot)
hot_percentage <- (hot_count / total_count) * 100
cat("Rows with 'Hot':", hot_count, "(", round(hot_percentage, 2), "%)\n")

df119_cold <- df119 %>%
  filter(if_any(everything(), ~ . == "Cold"))

cold_count <- nrow(df119_cold)
cold_percentage <- (cold_count / total_count) * 100
cat("Rows with 'Cold':", cold_count, "(", round(cold_percentage, 2), "%)\n")

# Calculate mean and median for the "Count_source" variable in df119
mean_value <- mean(df119$Count_source, na.rm = TRUE)
median_value <- median(df119$Count_source, na.rm = TRUE)

# Display the mean and median values
cat("Mean of Count_source:", mean_value, "\n")
cat("Median of Count_source:", median_value, "\n")

#####################
#####################
###########################
###### Distribution of Mutation-specific signatures
#####################
#####################

# Create df120 and its subsets based on "Genotype" and "PFC" values
df120 <- df101 %>% filter(Genotype == "CNV")
df120_PFC_1 <- df120 %>% filter(PFC == "1")
df120_PFC_2 <- df120 %>% filter(PFC == "2")
df120_PFC_3 <- df120 %>% filter(PFC == "3")

# Count the total rows for each subset
total_count_PFC_1 <- nrow(df120_PFC_1)
total_count_PFC_2 <- nrow(df120_PFC_2)
total_count_PFC_3 <- nrow(df120_PFC_3)

# Total row count of df120
total_count_df120 <- nrow(df120)

# Calculate the percentage of rows in each subset relative to df120
percentage_PFC_1 <- (total_count_PFC_1 / total_count_df120) * 100
percentage_PFC_2 <- (total_count_PFC_2 / total_count_df120) * 100
percentage_PFC_3 <- (total_count_PFC_3 / total_count_df120) * 100

# Display the counts and percentages in the console
cat("Total count of PFC = 1:", total_count_PFC_1, 
    "- Percentage:", round(percentage_PFC_1, 2), "%\n")
cat("Total count of PFC = 2:", total_count_PFC_2, 
    "- Percentage:", round(percentage_PFC_2, 2), "%\n")
cat("Total count of PFC = 3:", total_count_PFC_3, 
    "- Percentage:", round(percentage_PFC_3, 2), "%\n")

# Count and filter rows for SCS == "P" and SCS == "N" in each subset
df120_PFC_1_SCS_P <- df120_PFC_1 %>% filter(SCS == "P")
df120_PFC_2_SCS_P <- df120_PFC_2 %>% filter(SCS == "P")
df120_PFC_3_SCS_P <- df120_PFC_3 %>% filter(SCS == "P")

df120_PFC_1_SCS_N <- df120_PFC_1 %>% filter(SCS == "N")
df120_PFC_2_SCS_N <- df120_PFC_2 %>% filter(SCS == "N")
df120_PFC_3_SCS_N <- df120_PFC_3 %>% filter(SCS == "N")

# Total row counts for each subset
total_count_PFC_1 <- nrow(df120_PFC_1)
total_count_PFC_2 <- nrow(df120_PFC_2)
total_count_PFC_3 <- nrow(df120_PFC_3)

# Count rows for "P" and "N" in each subset
count_PFC_1_SCS_P <- nrow(df120_PFC_1_SCS_P)
count_PFC_2_SCS_P <- nrow(df120_PFC_2_SCS_P)
count_PFC_3_SCS_P <- nrow(df120_PFC_3_SCS_P)

count_PFC_1_SCS_N <- nrow(df120_PFC_1_SCS_N)
count_PFC_2_SCS_N <- nrow(df120_PFC_2_SCS_N)
count_PFC_3_SCS_N <- nrow(df120_PFC_3_SCS_N)

# Calculate percentages for "P" and "N" in each subset
percentage_PFC_1_SCS_P <- (count_PFC_1_SCS_P / total_count_PFC_1) * 100
percentage_PFC_2_SCS_P <- (count_PFC_2_SCS_P / total_count_PFC_2) * 100
percentage_PFC_3_SCS_P <- (count_PFC_3_SCS_P / total_count_PFC_3) * 100

percentage_PFC_1_SCS_N <- (count_PFC_1_SCS_N / total_count_PFC_1) * 100
percentage_PFC_2_SCS_N <- (count_PFC_2_SCS_N / total_count_PFC_2) * 100
percentage_PFC_3_SCS_N <- (count_PFC_3_SCS_N / total_count_PFC_3) * 100

# Display counts and percentages in the console
cat("PFC = 1 | SCS = 'P':", count_PFC_1_SCS_P, "- Percentage:", round(percentage_PFC_1_SCS_P, 2), "%\n")
cat("PFC = 2 | SCS = 'P':", count_PFC_2_SCS_P, "- Percentage:", round(percentage_PFC_2_SCS_P, 2), "%\n")
cat("PFC = 3 | SCS = 'P':", count_PFC_3_SCS_P, "- Percentage:", round(percentage_PFC_3_SCS_P, 2), "%\n")

cat("PFC = 1 | SCS = 'N':", count_PFC_1_SCS_N, "- Percentage:", round(percentage_PFC_1_SCS_N, 2), "%\n")
cat("PFC = 2 | SCS = 'N':", count_PFC_2_SCS_N, "- Percentage:", round(percentage_PFC_2_SCS_N, 2), "%\n")
cat("PFC = 3 | SCS = 'N':", count_PFC_3_SCS_N, "- Percentage:", round(percentage_PFC_3_SCS_N, 2), "%\n")

# Select rows where "Risky" appears in any column and calculate row count and percentage
df120_Risky <- df120 %>%
  filter(if_any(everything(), ~ . == "Risky"))

risky_count <- nrow(df120_Risky)
risky_percentage <- (risky_count / total_count_df120) * 100
cat("Rows with 'Risky':", risky_count, "(", round(risky_percentage, 2), "%)\n")

# Select rows where both "Protective" appear in any column
df120_Protective <- df120 %>%
  filter(if_any(everything(), ~ . == "Protective")) #& if_any(everything(), ~ . == "Hot"))

protective_count <- nrow(df120_Protective)
protective_percentage <- (protective_count / total_count_df120) * 100
cat("Rows with 'Protective':", protective_count, "(", round(protective_percentage, 2), "%)\n")

# Select rows where both "Protective" and "Cold" appear in any column
df120_cold <- df120 %>%
  filter(if_any(everything(), ~ . == "Cold"))

df120_hot <- df120 %>%
  filter(if_any(everything(), ~ . == "Hot"))

cold_count <- nrow(df120_cold)
cold_percentage <- (cold_count / total_count_df120) * 100
cat("Rows with 'Cold':", cold_count, "(", round(cold_percentage, 2), "%)\n")

hot_count <- nrow(df120_hot)
hot_percentage <- (hot_count / total_count_df120) * 100
cat("Rows with 'Hot':", hot_count, "(", round(hot_percentage, 2), "%)\n")

df120_variable <- df120 %>%
  filter(if_any(everything(), ~ . == "Variable"))
variable_count <- nrow(df120_variable)
variable_percentage <- (variable_count / total_count_df120) * 100
cat("Rows with 'Variable':", variable_count, "(", round(variable_percentage, 2), "%)\n")

df120_pro_tumoral <- df120 %>%
  filter(if_any(everything(), ~ . == "pro-tumoral"))
pro_tumoral_count <- nrow(df120_pro_tumoral)
pro_tumoral_percentage <- (pro_tumoral_count / total_count_df120) * 100
cat("Rows with 'pro_tumoral':", pro_tumoral_count, "(", round(pro_tumoral_percentage, 2), "%)\n")

df120_dual <- df120 %>%
  filter(if_any(everything(), ~ . == "dual"))
dual_count <- nrow(df120_dual)
dual_percentage <- (dual_count / total_count_df120) * 100
cat("Rows with 'dual':", dual_count, "(", round(dual_percentage, 2), "%)\n")

# Count the rows where microenvironment_classification is "anti-tumoral"
total_anti_tumoral <- sum(df120$microenvironment_classification == "anti-tumoral")
percentage_anti_tumoral <- (total_anti_tumoral / total_count_df120) * 100

df120_anti_tumoral <- df120 %>%
  filter(if_any(everything(), ~ . == "anti-tumoral"))
anti_tumoral_count <- nrow(df120_anti_tumoral)
anti_tumoral_percentage <- (anti_tumoral_count / total_count_df120) * 100
cat("Rows with 'anti_tumoral':", anti_tumoral_count, "(", round(anti_tumoral_percentage, 2), "%)\n")

#####################
#####################
#####################
#####################
# Select rows where both "Risky" appear in any column
df120_Risky_hot <- df120 %>%
  filter(if_any(everything(), ~ . == "Risky") & if_any(everything(), ~ . == "Hot"))

risky_hot_count <- nrow(df120_Risky_hot)
risky_hot_percentage <- (risky_hot_count / total_count) * 100
cat("Rows with 'Risky' and 'Hot':", risky_hot_count, "(", round(risky_hot_percentage, 2), "%)\n")

risky_hot_count <- nrow(df120_Risky_hot)
risky_hot_percentage <- (risky_hot_count / total_count) * 100
cat("Rows with 'Risky' and 'Hot':", risky_hot_count, "(", round(risky_hot_percentage, 2), "%)\n")

# Select rows where both "Risky" and "Cold" appear in any column
df120_Risky_cold <- df120 %>%
  filter(if_any(everything(), ~ . == "Risky") & if_any(everything(), ~ . == "Cold"))

risky_cold_count <- nrow(df120_Risky_cold)
risky_cold_percentage <- (risky_cold_count / total_count) * 100
cat("Rows with 'Risky' and 'Cold':", risky_cold_count, "(", round(risky_cold_percentage, 2), "%)\n")

# Calculate row counts and percentages for other categories
df120_hot <- df120 %>%
  filter(if_any(everything(), ~ . == "Hot"))

hot_count <- nrow(df120_hot)
hot_percentage <- (hot_count / total_count) * 100
cat("Rows with 'Hot':", hot_count, "(", round(hot_percentage, 2), "%)\n")

df120_cold <- df120 %>%
  filter(if_any(everything(), ~ . == "Cold"))

cold_count <- nrow(df120_cold)
cold_percentage <- (cold_count / total_count_df120) * 100
cat("Rows with 'Cold':", cold_count, "(", round(cold_percentage, 2), "%)\n")

df120_variable <- df120 %>%
  filter(if_any(everything(), ~ . == "Variable"))

variable_count <- nrow(df120_variable)
variable_percentage <- (variable_count / total_count) * 100
cat("Rows with 'Variable':", variable_count, "(", round(variable_percentage, 2), "%)\n")

df120_antitumoral <- df120 %>%
  filter(if_any(everything(), ~ . == "anti-tumoral"))

antitumoral_count <- nrow(df120_antitumoral)
antitumoral_percentage <- (antitumoral_count / total_count_df120) * 100
cat("Rows with 'Anti-tumoral':", antitumoral_count, "(", round(antitumoral_percentage, 2), "%)\n")

df120_protumoral <- df120 %>%
  filter(if_any(everything(), ~ . == "pro-tumoral"))

protumoral_count <- nrow(df120_protumoral)
protumoral_percentage <- (protumoral_count / total_count) * 100
cat("Rows with 'Pro-tumoral':", protumoral_count, "(", round(protumoral_percentage, 2), "%)\n")

df120_dual <- df120 %>%
  filter(if_any(everything(), ~ . == "dual"))

dual_count <- nrow(df120_dual)
dual_percentage <- (dual_count / total_count) * 100
cat("Rows with 'Dual':", dual_count, "(", round(dual_percentage, 2), "%)\n")

df120_hot <- df120 %>%
  filter(if_any(everything(), ~ . == "Hot"))

hot_count <- nrow(df120_hot)
hot_percentage <- (hot_count / total_count) * 100
cat("Rows with 'Hot':", hot_count, "(", round(hot_percentage, 2), "%)\n")

df120_cold <- df120 %>%
  filter(if_any(everything(), ~ . == "Cold"))

cold_count <- nrow(df120_cold)
cold_percentage <- (cold_count / total_count) * 100
cat("Rows with 'Cold':", cold_count, "(", round(cold_percentage, 2), "%)\n")

# Calculate mean and median for the "Count_source" variable in df120
mean_value <- mean(df120$Count_source, na.rm = TRUE)
median_value <- median(df120$Count_source, na.rm = TRUE)

# Display the mean and median values
cat("Mean of Count_source:", mean_value, "\n")
cat("Median of Count_source:", median_value, "\n")


#####################
#####################
#####################
#####################
# Filter spcefic rows from df101 based on the specified conditions, ensuring multiple terms match exactly within a structured string
df_filtered <- df101[df101$PFC == "3" & 
                       df101$SCS == "P" & 
                       df101$Genotype == "mRNA" & 
                       grepl("\\bPOU5F1\\b|\\bSOX2\\b|\\bNANOG\\b", df101$Signature), ]

# Display the filtered dataframe
df_filtered

#####################
#####################
# Function to extract unique terms from the Signature column
extract_terms <- function(signature_value) {
  # Remove any unnecessary symbols (e.g., parentheses or plus signs) and split the string into individual terms
  terms <- unlist(strsplit(gsub("[()]", "", signature_value), " \\+ "))
  # Return unique terms
  return(unique(terms))
}

# Step 1: Extract unique terms from the Signature column for each row
terms_list <- unique(unlist(lapply(df101$Signature, extract_terms)))

# Step 2: Filter the dataframe based on PFC, SCS, and Genotype conditions
df_filtered_conditions <- df101[df101$PFC == "3" & 
                                  df101$SCS == "P" & 
                                  df101$Genotype == "mRNA", ]

# Step 3: Initialize an empty logical vector to store row matches
row_matches <- rep(FALSE, nrow(df_filtered_conditions))

# Step 4: Iteratively search for each term in the filtered dataframe
for (term in terms_list) {
  row_matches <- row_matches | grepl(paste0("\\b", term, "\\b"), df_filtered_conditions$Signature)
}

# Step 5: Filter the dataframe based on row matches
df_final_filtered_mRNA <- df_filtered_conditions[row_matches, ]

######
######
######
######
# Group by 'Signature' and count the number of rows per 'Signature'
signature_distribution_mRNA <- df_final_filtered_mRNA %>%
  group_by(Signature) %>%
  summarise(count = n())

# Display the distribution summary
summary(signature_distribution_mRNA$count)

# Optional: Plot the distribution for better visualization
ggplot(signature_distribution_mRNA, aes(x = count)) +
  geom_histogram(binwidth = 1, fill = 'blue', color = 'black') +
  labs(title = "Distribution of Rows per 'Signature' Value",
       x = "Number of Rows per Signature",
       y = "Frequency") +
  theme_minimal()

dev.off()

######
######
######
######
# Function to extract unique terms from the Signature column
extract_terms <- function(signature_value) {
  # Remove any unnecessary symbols (e.g., parentheses or plus signs) and split the string into individual terms
  terms <- unlist(strsplit(gsub("[()]", "", signature_value), " \\+ "))
  # Return unique terms
  return(unique(terms))
}

# Step 1: Extract unique terms from the Signature column for each row
terms_list <- unique(unlist(lapply(df101$Signature, extract_terms)))

# Step 2: Filter the dataframe based on PFC, SCS, and Genotype conditions
df_filtered_conditions <- df101[df101$PFC == "3" & 
                                  df101$SCS == "P" & 
                                  df101$Genotype == "miRNA", ]

# Step 3: Initialize an empty logical vector to store row matches
row_matches <- rep(FALSE, nrow(df_filtered_conditions))

# Step 4: Iteratively search for each term in the filtered dataframe
for (term in terms_list) {
  row_matches <- row_matches | grepl(paste0("\\b", term, "\\b"), df_filtered_conditions$Signature)
}

# Step 5: Filter the dataframe based on row matches
df_final_filtered_miRNA <- df_filtered_conditions[row_matches, ]

######
######
######
######
######
# Group by 'Signature' and count the number of rows per 'Signature'
signature_distribution_miRNA <- df_final_filtered_miRNA %>%
  group_by(Signature) %>%
  summarise(count = n())

# Display the distribution summary
summary(signature_distribution_miRNA$count)

# Optional: Plot the distribution for better visualization
ggplot(signature_distribution_miRNA, aes(x = count)) +
  geom_histogram(binwidth = 1, fill = 'blue', color = 'black') +
  labs(title = "Distribution of Rows per 'Signature' Value",
       x = "Number of Rows per Signature",
       y = "Frequency") +
  theme_minimal()

dev.off()

#####
#####
#####
#####

# Function to extract unique terms from the Signature column
extract_terms <- function(signature_value) {
  # Remove any unnecessary symbols (e.g., parentheses or plus signs) and split the string into individual terms
  terms <- unlist(strsplit(gsub("[()]", "", signature_value), " \\+ "))
  # Return unique terms
  return(unique(terms))
}

# Step 1: Extract unique terms from the Signature column for each row for PFC == 3 and SCS =="P"
terms_list <- unique(unlist(lapply(df101$Signature, extract_terms)))

# Define the Genotype values
genotype_values <- c("Mutation", "mRNA", "Methylation", "Transcript", "Protein", "CNV", "miRNA")

# Initialize a list to store results
results_list <- list()

# Step 2: Loop through each genotype and perform the analysis
for (genotype in genotype_values) {
  # Step 2.1: Filter the dataframe based on PFC, SCS, and the current Genotype
  df_filtered_conditions <- df101[df101$PFC == "3" & 
                                    df101$SCS == "P" & 
                                    df101$Genotype == genotype, ]
  
  # Step 2.2: Initialize an empty logical vector to store row matches
  row_matches <- rep(FALSE, nrow(df_filtered_conditions))
  
  # Step 2.3: Iteratively search for each term in the filtered dataframe
  for (term in terms_list) {
    row_matches <- row_matches | grepl(paste0("\\b", term, "\\b"), df_filtered_conditions$Signature)
  }
  
  # Step 2.4: Filter the dataframe based on row matches
  df_final_filtered <- df_filtered_conditions[row_matches, ]
  
  # Step 2.5: Group by 'Signature' and count the number of rows per 'Signature'
  signature_distribution <- df_final_filtered %>%
    group_by(Signature) %>%
    summarise(count = n())
  
  # Step 2.6: Save the results in the list, using the Genotype as a key
  results_list[[genotype]] <- list(
    filtered_data = df_final_filtered,
    signature_distribution = signature_distribution
  )
  
  # Step 2.7: Display the summary for the current Genotype
  cat("\nSummary for Genotype:", genotype, "\n")
  print(summary(signature_distribution$count))
  
  # Step 2.8: Plot the distribution for better visualization
  p <- ggplot(signature_distribution, aes(x = count)) +
    geom_histogram(binwidth = 1, fill = 'blue', color = 'black') +
    labs(title = paste("Distribution of Rows per 'Signature' for", genotype),
         x = "Number of Rows per Signature",
         y = "Frequency") +
    theme_minimal()
  
  # Display the plot
  print(p)

dev.off() 

  # Save the plot as an image file
  ggsave(paste0("signature_distribution_", genotype, ".png"), plot = p)
}

# The list 'results_list' contains dataframes and summary results for each genotype

# Convert the filtered_data under miRNA into a dataframe
df_Protein_2 <- as.data.frame(results_list[["Protein"]][["filtered_data"]])

# Convert the filtered_data under miRNA into a dataframe
df_Mutation_2 <- as.data.frame(results_list[["Mutation"]][["filtered_data"]])

# Convert the filtered_data under miRNA into a dataframe
df_CNV_2 <- as.data.frame(results_list[["CNV"]][["filtered_data"]])

# Convert the filtered_data under miRNA into a dataframe
df_miRNA_2 <- as.data.frame(results_list[["miRNA"]][["filtered_data"]])

# Convert the filtered_data under miRNA into a dataframe
df_Transcript_2 <- as.data.frame(results_list[["Transcript"]][["filtered_data"]])

# Convert the filtered_data under miRNA into a dataframe
df_mRNA_2 <- as.data.frame(results_list[["mRNA"]][["filtered_data"]])

# Convert the filtered_data under miRNA into a dataframe
df_Methylation_2 <- as.data.frame(results_list[["Methylation"]][["filtered_data"]])

##########
##########
##########
# Extract rows where the 'Genotype' column is equal to "mRNA"
# Using dplyr to filter the dataframe
df101_mRNA <- df101 %>%
  filter(Genotype == "mRNA")

# Function to extract unique terms from the Signature column
extract_terms_4 <- function(signature_value) {
  # Remove any unnecessary symbols (e.g., parentheses or plus signs) and split the string into individual terms
  terms <- unlist(strsplit(gsub("[()]", "", signature_value), " \\+ "))
  # Return unique terms
  return(unique(terms))
}

# Step 1: Extract unique terms from the Signature column for each row for PFC == 3 and SCS =="P"
terms_list_4 <- unique(unlist(lapply(df101_mRNA$Signature, extract_terms_4)))

# Apply trim function to remove any leading or trailing spaces from the terms

terms_list_4_trimmed <- trimws(terms_list_4)

terms_list_4 <- terms_list_4_trimmed 

# Convert the terms list to a dataframe and remove any duplicates
terms_df_4 <- data.frame(Terms = unique(terms_list_4_trimmed))

# Step 1: Remove backticks from the 'Term' column
terms_df_4$Terms <- gsub("`", "", terms_df_4$Terms)

# Step 3: Rename the 'Term' column to 'Genes'
colnames(terms_df_4)[colnames(terms_df_4) == "Terms"] <- "Genes"

# Convert the terms list to a dataframe and remove any duplicates
terms_df_4 <- data.frame(Genes = unique(terms_df_4))

# Step 4: Remove backticks from the 'Genes' column
terms_df_4$Genes <- gsub("`", "", terms_df_4$Genes)

# Step 5: Sort the 'Genes' column alphabetically
terms_df_4 <- terms_df_4 %>% arrange(Genes)

Effective_target_genes <- terms_df_4

# rio::export the dataframe to an Excel file
rio::export(Effective_target_genes, "Effective_target_genes.xlsx")

#######
#######
#######
# Extract rows where the 'Genotype' column is equal to "mRNA", "Methylation", "CNV", or "Mutation"
df101_MMCM <- df101 %>%
  filter(Genotype %in% c("mRNA", "Methylation", "CNV", "Mutation"))

# Function to extract unique terms from the Signature column
extract_terms_5 <- function(signature_value) {
  # Remove any unnecessary symbols (e.g., parentheses or plus signs) and split the string into individual terms
  terms <- unlist(strsplit(gsub("[()]", "", signature_value), " \\+ "))
  # Return unique terms
  return(unique(terms))
}

# Step 1: Extract unique terms from the Signature column for each row
terms_list_5 <- unique(unlist(lapply(df101_MMCM$Signature, extract_terms_5)))

# Apply trim function to remove any leading or trailing spaces from the terms

terms_list_5_trimmed <- trimws(terms_list_5)

terms_list_5 <- terms_list_5_trimmed 

# Convert the terms list to a dataframe and remove any duplicates
terms_df_5 <- data.frame(Terms = unique(terms_list_5_trimmed))

# Step 1: Remove backticks from the 'Term' column
terms_df_5$Terms <- gsub("`", "", terms_df_5$Terms)

# Step 3: Rename the 'Term' column to 'Genes'
colnames(terms_df_5)[colnames(terms_df_5) == "Terms"] <- "Genes"

# Convert the terms list to a dataframe and remove any duplicates
terms_df_5 <- data.frame(Genes = unique(terms_df_5))

# Step 4: Remove backticks from the 'Genes' column
terms_df_5$Genes <- gsub("`", "", terms_df_5$Genes)

# Step 5: Sort the 'Genes' column alphabetically
terms_df_5 <- terms_df_5 %>% arrange(Genes)

# Convert the terms list to a dataframe and remove any duplicates
terms_df_5 <- data.frame(Genes = unique(terms_df_5))

Effective_target_genes_5 <- terms_df_5

  # rio::export the dataframe to an Excel file
  rio::export(Effective_target_genes, "Effective_target_genes_5.xlsx")
  
  df140 <- import("Itherapy_target_genes.xlsx")
  
  # Step 1: Trim leading and trailing spaces in the 'Gene_symbol_EMA' column
  df140$Gene_symbol_EMA <- trimws(df140$Gene_symbol_EMA)
  
  # Step 2: Create a list "IT_gene_symbols" with the values under the "Gene_Symbol_EMA" column
  IT_gene_symbols <- as.list(df140$Gene_symbol_EMA)
  
  # Step 3: Check which values from "IT_gene_symbols" are present in the "Genes" column of df "Effective_target_genes"
  common_genes_5 <- IT_gene_symbols[IT_gene_symbols %in% Effective_target_genes_5$Genes]
  
  # Calculate the percentage of common values
  total_genes_5 <- length(IT_gene_symbols)  # Total number of values in IT_gene_symbols
  common_genes_count_5 <- length(common_genes_5)  # Number of common values
  
  # Estimate the percentage of common values
  percentage_common_5 <- (common_genes_count_5 / total_genes_5) * 100
  
  # Print the percentage
  percentage_common_5

# Extract rows with values equal to "Hot", and "anti-tumoral
df190_MHP <- df_Methylation %>%
  filter(immune_classification %in% "Hot", 
         microenvironment_classification %in% "anti-tumoral") 

# Extract rows where 'immune_classification' is "Hot" and 'microenvironment_classification' is "anti-tumoral"
df190_MHP <- df_Methylation %>%
  filter(immune_classification == "Hot", 
         microenvironment_classification == "anti-tumoral",
         SCS == "N",
         PFC != "3")

rio::export(df190_MHP,"df190_MHP.xlsx")  

# Extract rows where 'immune_classification' is "Hot" and 'microenvironment_classification' is "anti-tumoral"
df191_MRHP <- df_mRNA %>%
  filter(immune_classification == "Hot", 
         microenvironment_classification == "anti-tumoral",
         SCS == "N",
         PFC != "3")

rio::export(df191_MRHP,"df191_MRHP.xlsx")  

# Extract rows where 'immune_classification' is "Hot" and 'microenvironment_classification' is "anti-tumoral"
df192_PN3HAT <- df_Protein %>%
  filter(immune_classification == "Hot", 
         microenvironment_classification == "anti-tumoral",
         SCS == "N",
         PFC == "3")

# Extract rows where 'immune_classification' is "Hot" and 'microenvironment_classification' is "anti-tumoral"
df192_PN2HAT <- df_Protein %>%
  filter(immune_classification == "Hot", 
         microenvironment_classification == "anti-tumoral",
         SCS == "N",
         PFC == "2")

rio::export(df192_PN2HAT,"df192_PN2HAT.xlsx")  

# Extract rows where 'immune_classification' is "Hot" and 'microenvironment_classification' is "anti-tumoral"
df192_PN1HAT <- df_Protein %>%
  filter(immune_classification == "Hot", 
         microenvironment_classification == "anti-tumoral",
         SCS == "N",
         PFC == "1")

rio::export(df192_PN1HAT,"df192_PN1HAT.xlsx")  

######
######
######
# Extract rows where 'immune_classification' is "Hot" and 'microenvironment_classification' is "anti-tumoral"
df192_PP3HAT <- df_Protein %>%
  filter(immune_classification == "Hot", 
         microenvironment_classification == "anti-tumoral",
         SCS == "P",
         PFC == "3")

# Extract rows where 'immune_classification' is "Hot" and 'microenvironment_classification' is "anti-tumoral"
df192_PP2HAT <- df_Protein %>%
  filter(immune_classification == "Hot", 
         microenvironment_classification == "anti-tumoral",
         SCS == "P",
         PFC == "2")

rio::export(df192_PP2HAT,"df192_PP2HAT.xlsx")  

# Extract rows where 'immune_classification' is "Hot" and 'microenvironment_classification' is "anti-tumoral"
df192_PP1HAT <- df_Protein %>%
  filter(immune_classification == "Hot", 
         microenvironment_classification == "anti-tumoral",
         SCS == "P",
         PFC == "1")

rio::export(df192_PP1HAT,"df192_PP1HAT.xlsx") 

# # Get all objects in the environment
# objs <- ls()
# 
# # Filter the objects to find those that start with 'df192_'
# df_to_remove <- objs[grep("^df192_", objs)]
# 
# # Remove the filtered data frames
# rm(list = df_to_remove)

# Extract rows where 'immune_classification' is "Hot" and 'microenvironment_classification' is "anti-tumoral"
df193_MN3HAT <- df_Mutation %>%
  filter(immune_classification == "Hot", 
         microenvironment_classification == "anti-tumoral",
         SCS == "N",
         PFC == "3")

rio::export(df193_MN3HAT,"df193_MN3HAT.xlsx") 

# Extract rows where 'immune_classification' is "Hot" and 'microenvironment_classification' is "anti-tumoral"
df193_MN2HAT <- df_Mutation %>%
  filter(immune_classification == "Hot", 
         microenvironment_classification == "anti-tumoral",
         SCS == "N",
         PFC == "2")

rio::export(df193_MN2HAT,"df193_MN2HAT.xlsx")  

# Extract rows where 'immune_classification' is "Hot" and 'microenvironment_classification' is "anti-tumoral"
df193_MN1HAT <- df_Mutation %>%
  filter(immune_classification == "Hot", 
         microenvironment_classification == "anti-tumoral",
         SCS == "N",
         PFC == "1")

rio::export(df193_MN1HAT,"df193_MN1HAT.xlsx")  

######
######
######
# Extract rows where 'immune_classification' is "Hot" and 'microenvironment_classification' is "anti-tumoral"
df193_MP3HAT <- df_Mutation %>%
  filter(immune_classification == "Hot", 
         microenvironment_classification == "anti-tumoral",
         SCS == "P",
         PFC == "3")

# Extract rows where 'immune_classification' is "Hot" and 'microenvironment_classification' is "anti-tumoral"
df193_MP2HAT <- df_Mutation %>%
  filter(immune_classification == "Hot", 
         microenvironment_classification == "anti-tumoral",
         SCS == "P",
         PFC == "2")

rio::export(df193_MP2HAT,"df193_MP2HAT.xlsx")  

# Extract rows where 'immune_classification' is "Hot" and 'microenvironment_classification' is "anti-tumoral"
df193_MP1HAT <- df_Mutation %>%
  filter(immune_classification == "Hot", 
         microenvironment_classification == "anti-tumoral",
         SCS == "P",
         PFC == "1")

rio::export(df193_MP1HAT,"df193_MP1HAT.xlsx") 

######
######
######
# Extract rows where 'immune_classification' is "Hot" and 'microenvironment_classification' is "anti-tumoral"
df194_CNVP3HAT <- df_CNV %>%
  filter(immune_classification == "Hot", 
         microenvironment_classification == "anti-tumoral",
         SCS == "P",
         PFC == "3")

# Extract rows where 'immune_classification' is "Hot" and 'microenvironment_classification' is "anti-tumoral"
df194_CNVP2HAT <- df_CNV %>%
  filter(immune_classification == "Hot", 
         microenvironment_classification == "anti-tumoral",
         SCS == "P",
         PFC == "2")

rio::export(df194_CNVP2HAT,"df194_CNVP2HAT.xlsx")  

# Extract rows where 'immune_classification' is "Hot" and 'microenvironment_classification' is "anti-tumoral"
df194_CNVP1HAT <- df_CNV %>%
  filter(immune_classification == "Hot", 
         microenvironment_classification == "anti-tumoral",
         SCS == "P",
         PFC == "1")

rio::export(df194_CNVP1HAT,"df194_CNVP1HAT.xlsx") 

# Extract rows where 'immune_classification' is "Hot" and 'microenvironment_classification' is "anti-tumoral"
df195 <- df101 %>%
  filter(Signature == "TP53", 
         PFC == "3",
         SCS == "N"
  )

rio::export(df194_CNVP1HAT,"df194_CNVP1HAT.xlsx")

#####
#####
#####"Counting and Calculating Percentage of RCD_split Occurrences in RCD_types with String Separation
# Split the values in the "RCD_types" column of df101 into individual components
df101_split <- df101 %>%
  separate_rows(RCD_types, sep = "/")

df196 <- df_RCD_frequencies

# Count the occurrences of each "RCD_split" value in "RCD_types"
df_counts <- df196 %>%
  select(RCD_split) %>% # Focus on "RCD_split" column from df196
  left_join(
    df101_split %>%
      count(RCD_types) %>%
      rename(Count = n),
    by = c("RCD_split" = "RCD_types")
  )

# Replace NA with 0 where counts were not found
df_counts$Count[is.na(df_counts$Count)] <- 0

# Calculate the percentage
total_rows <- nrow(df101) # Total number of rows in df101
df_counts <- df_counts %>%
  mutate(Percentage = (Count / total_rows) * 100)

# Create the final dataframe with "RCD_split", "Count", and "Percentage"
df_final <- df_counts %>%
  select(RCD_split, Count, Percentage)

rio::export(df_final, "Counting and Calculating Percentage of RCD_split Occurrences in RCD_types with String Separation.xlsx")

# Normalize the 'Count' column to avoid outlier effects
df_final$Normalized_Count <- df_final$Count / sum(df_final$Count)

# Reorder the RCD_split factor levels based on the normalized count in descending order
df_final$RCD_split <- factor(df_final$RCD_split, levels = df_final$RCD_split[order(df_final$Normalized_Count, decreasing = TRUE)])

# Create the histogram with the reordered RCD types
ggplot(df_final, aes(x = RCD_split, y = Normalized_Count)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(title = "Normalized Frequency of Instances of Each RCD Type",
       x = "RCD Types",
       y = "Normalized Frequency") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()

# Save the plot as an image (optional)
ggsave("RCD_Types_Normalized_Histogram_Descending.png", width = 8, height = 6)

# Load necessary libraries
library(ggplot2)

# Apply log10 normalization to the 'Count' column to avoid outlier effects
df_final$Log10_Count <- log10(df_final$Count + 1)  # Add 1 to avoid log(0) issues

# Reorder the RCD_split factor levels based on the log10 count in descending order
df_final$RCD_split <- factor(df_final$RCD_split, levels = df_final$RCD_split[order(df_final$Log10_Count, decreasing = TRUE)])

# Create the histogram with the reordered RCD types based on log10 normalization
ggplot(df_final, aes(x = RCD_split, y = Log10_Count)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(title = "Log10 Normalized Frequency of Instances of Each RCD Type",
       x = "RCD Types",
       y = "Log10 Normalized Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()

# Save the plot as an image (optional)
ggsave("RCD_Types_Log10_Normalized_Histogram_Descending.png", width = 8, height = 6)

### Global estimation of the distribution of gene and RCD forms in the final signatures

# Step 0: Load the full ranking dataframe
df200 <- final_signatures_full_ranking

# Step 1: Split the values in "Signature" and "RCD_types" in df200
df200 <- df200 %>%
  dplyr::mutate(Signature_Split = strsplit(as.character(Signature), " \\+ "),    # Split Signature by " + "
                RCD_types_Split = strsplit(as.character(RCD_types), "/"))        # Split RCD_types by "/"

# Step 2: Unnest the split values to get individual gene symbols and RCD types into separate rows
df200_split <- df200 %>%
  tidyr::unnest(Signature_Split) %>%           # Unnest Signature to have individual genes
  tidyr::unnest(RCD_types_Split)               # Unnest RCD_types to have individual RCD forms

# Step 3: Filter to keep only rows where gene symbols in df200 match Target_genes$Gene_Symbol
# Trim spaces and remove unwanted characters without changing case
matching_genes_df <- df200_split %>%
  dplyr::mutate(Signature_Split = trimws(gsub("[()`]", "", Signature_Split))) %>%  # Remove unwanted characters
  dplyr::filter(trimws(gsub("[()`]", "", Signature_Split)) %in% trimws(gsub("[()`]", "", Target_genes$Gene_Symbol)))

# Step 4: Estimate the distribution of genes by RCD form
# Count how often each gene appears for each RCD type
gene_distribution <- matching_genes_df %>%
  dplyr::group_by(Signature_Split, RCD_types_Split) %>%
  dplyr::summarize(count = n(), .groups = 'drop')   # Count the number of occurrences per gene and RCD form

# Step 5: Extract unique gene symbols and store them in a dataframe
unique_df200_matching_genes <- gene_distribution %>%
  dplyr::distinct(Signature_Split) %>%           # Extract unique gene symbols
  dplyr::rename(Genes = Signature_Split)         # Rename the column to "Genes"

# Step 6: Summarize the overall distribution of genes by RCD form
rcd_distribution <- gene_distribution %>%
  dplyr::group_by(RCD_types_Split) %>%
  dplyr::summarize(total_genes = n_distinct(Signature_Split),   # Count distinct genes per RCD type
                   total_occurrences = sum(count))              # Sum the occurrences of genes per RCD form

# Step 7: Count the total number of unique gene symbols in the gene_distribution
unique_gene_count <- dplyr::n_distinct(gene_distribution$Signature_Split)

# Display the results as a list
list(unique_df200_matching_genes = unique_df200_matching_genes, 
     rcd_distribution = rcd_distribution, 
     unique_gene_count = unique_gene_count)

df300 <- final_signatures_full_ranking[final_signatures_full_ranking$Nomenclature == "HNSC-308.5.3.N.3.0.0.3.2.3", ]

#### Visualizing Cancer Type and Genomic Feature Distributions cocnentric circle
# Sample Data Preprocessing
# Count occurrences of each Genotype class within each CTAB group
df_plot <- df101 %>%
  group_by(CTAB, Genotype) %>%
  summarise(Frequency = n()) %>%
  ungroup()

# Plotting Concentric Circle Plot
ggplot(df_plot, aes(x0 = 0, y0 = 0, r = Frequency, fill = Genotype)) +
  geom_arc_bar(aes(start = 0, end = 2 * pi, r0 = as.numeric(factor(CTAB)), r = as.numeric(factor(CTAB)) + 0.8)) +
  scale_fill_brewer(palette = "Set3") +
  theme_void() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    title = "Concentric Circle Plot: Distribution of Genotype Classes by CTAB",
    fill = "Genotype"
  )

dev.off()
#####
#####
#####
#####
##### Visualizing Cancer Type and Genomic Feature Distributions with Okabe-Ito Palette for Enhanced Accessibility in Data Interpretation
#####
#####
##########
# Load required packages
library(ggplot2)
library(dplyr)

# Define Okabe-Ito palette (color-blind friendly)
okabe_ito_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                       "#0072B2", "#D55E00", "#CC79A7")

# Data Preprocessing with Renamed Columns
df_plot <- df101 %>%
  group_by(`Cancer Type Abbreviation` = CTAB, `Genomic Feature` = Genotype) %>%
  summarise(Absolute_Count = n()) %>%
  ungroup()

# Plotting Accumulated Histogram with Absolute Counts and color-blind friendly palette
signature_distribution_plot <- ggplot(df_plot, aes(x = `Cancer Type Abbreviation`, y = Absolute_Count, fill = `Genomic Feature`)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = okabe_ito_palette) +  # Using Okabe-Ito color-blind palette
  theme_classic() + # Ensures a white background
  theme(
    panel.background = element_rect(fill = "white", color = "white"), # Set background to white
    plot.background = element_rect(fill = "white", color = "white"),   # Set plot background to white
    axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5), 
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 14)
  ) +
  labs(
    x = "Cancer Type Abbreviation",
    y = "Absolute Count (Accumulated)",
    fill = "Genomic Feature"
  )
signature_distribution_plot

dev.off()

# Corrected ggsave function usage
ggsave("Signature_Distribution.tif", plot = signature_distribution_plot, 
       width = 29, height = 21, units = "cm", dpi = 600)

#####
#####
#####
#####
# Creating seldom dfs by Genotyp    e values 
# Assuming df101 is already loaded into the environment
# Get the unique values in the "Genotype" column
genotype_values <- unique(df101$Genotype)

# Loop through each unique genotype value to create separate dataframes
for (genotype in genotype_values) {
  # Create a subset of df101 where Genotype equals the current genotype value
  subset_df <- subset(df101, Genotype == genotype)
  
  # Assign the subset to a new dataframe named after the genotype value
  # Replace any spaces or special characters in the genotype value with underscores
  assign(paste0("df101_", gsub("[^[:alnum:]]", "_", genotype)), subset_df)
}

##### Analysis of df_mRNA

# Total number of rows in df101_mRNA for calculating percentages
df101_total_count <- nrow(df101_mRNA)

# Filter for SCS == "N" and calculate count and percentage
df101_mRNA_N <- df101_mRNA %>%
  filter(SCS == "N") %>%
  summarise(count_N = n(),
            percentage_N = (count_N / df101_total_count) * 100)

# Filter for SCS == "P" and calculate count and percentage
df101_mRNA_P <- df101_mRNA %>%
  filter(SCS == "P") %>%
  summarise(count_P = n(),
            percentage_P = (count_P / df101_total_count) * 100)

# Display the results
df101_mRNA_N
df101_mRNA_P

df101_mRNA_N1 <- df101_mRNA %>%
  filter(SCS == "N",
         PFC == "1") %>%
  summarise(count_N1 = n(),
            percentage_N1 = (count_N1 / df101_total_count) * 100)

df101_mRNA_N2 <- df101_mRNA %>%
  filter(SCS == "N",
         PFC == "2") %>%
  summarise(count_N2 = n(),
            percentage_N2 = (count_N2 / df101_total_count) * 100)

df101_mRNA_N3 <- df101_mRNA %>%
  filter(SCS == "N",
         PFC == "3") %>%
  summarise(count_N3 = n(),
            percentage_N3 = (count_N3 / df101_total_count) * 100)

##### Analysis of df_mRNA

# Calculate the total count of rows in df101_mRNA for percentage calculation
df101_total_count <- nrow(df101_mRNA)

# Filter for SCS == "N", PFC == "3", and any "Risky" in the specified columns
df101_mRNA_N3R <- df101_mRNA %>%
  filter(SCS == "N",
         PFC == "3",
         Type_Cox_DSS == "Risky" | Type_Cox_DFI == "Risky" | Type_Cox_PFI == "Risky" | Type_Cox_OS == "Risky") %>%
  summarise(count_N3R = n(),
            percentage_N3R = (count_N3R / df101_total_count) * 100)

# Display the result
df101_mRNA_N3R

# Filter for SCS == "N", PFC == "3", and any "Protective" in the specified columns
df101_mRNA_N3P <- df101_mRNA %>%
  filter(SCS == "N",
         PFC == "3",
         Type_Cox_DSS == "Protective" | Type_Cox_DFI == "Protective" | Type_Cox_PFI == "Protective" | Type_Cox_OS == "Protective") %>%
  summarise(count_N3P = n(),
            percentage_N3P = (count_N3P / df101_total_count) * 100)

# Display the result
df101_mRNA_N3P

#### Protumoral antitumor mRNA signatures

df101_mRNA_antitumoral <- df101_mRNA %>%
  filter(microenvironment_classification == "anti-tumoral") %>%
  summarise(count_AT = n(),
            percentage_AT = (count_AT / df101_total_count) * 100)

df101_mRNA_protumoral <- df101_mRNA %>%
  filter(microenvironment_classification == "pro-tumoral") %>%
  summarise(count_PT = n(),
            percentage_PT = (count_PT / df101_total_count) * 100)

df101_mRNA_dual <- df101_mRNA %>%
  filter(microenvironment_classification == "dual") %>%
  summarise(count_DU = n(),
            percentage_DU = (count_DU/ df101_total_count) * 100)

df101_mRNA_hot <- df101_mRNA %>%
  filter(immune_classification == "Hot") %>%
  summarise(count_hot = n(),
            percentage_hot = (count_hot/ df101_total_count) * 100)

df101_mRNA_cold <- df101_mRNA %>%
  filter(immune_classification == "Cold") %>%
  summarise(count_cold = n(),
            percentage_cold = (count_cold/ df101_total_count) * 100)

df101_mRNA_Variable <- df101_mRNA %>%
  filter(immune_classification == "Variable") %>%
  summarise(count_Variable = n(),
            percentage_Variable = (count_Variable/ df101_total_count) * 100)

#####
#####
# Load necessary libraries
library(dplyr)
library(tidyr)

# Assuming 'df101' is already loaded with a column named 'Signature'
df101 <-  final_signatures_full_ranking

# Step 1: Calculate the number of elements in each row of the "Signature" column without altering it
df101 <- df101 %>%
  mutate(
    # Clean the content of Signature without modifying the original column
    cleaned_signature = gsub("[()]", "", Signature),         # Remove parentheses
    cleaned_signature = gsub("`", "", cleaned_signature),    # Remove backticks
    cleaned_signature = str_split(cleaned_signature, "\\s*\\+\\s*"), # Split by "+" with optional spaces
    Member_counts = sapply(cleaned_signature, length)        # Count the number of elements
  ) %>%
  select(-cleaned_signature) %>%                             # Remove the temporary cleaned column
  relocate(Member_counts, .after = Signature)                # Place Member_counts immediately after Signature

# Assuming 'df101' already has the variables 'Member_counts' and 'Genotype'

# Calculate the minimum and maximum range of Member_counts for each unique value in Genotype
df101_range_summary_per_genotype <- df101 %>%
  group_by(Genotype) %>%
  summarise(
    Min_Member_counts = min(Member_counts, na.rm = TRUE),  # Minimum value per Genotype
    Max_Member_counts = max(Member_counts, na.rm = TRUE)   # Maximum value per Genotype
  ) %>%
  ungroup()

# View the result
print(df101_range_summary_per_genotype)

rio::export(df101_range_summary_per_genotype, "Signature member counts per genotype.xlsx")

# install.packages("devtools")

# # Install release from GitHub:
# devtools::install_github("stefanedwards/lemon", ref='v0.3.1')
# 
# # Or get the lastest development version from GitHub:
# devtools::install_github("stefanedwards/lemon")
# 
# 
###
###
### Making Figure S2. Distribution of signature member element counts across omic features.
### The element counts are displayed in order of range, from smallest to largest, 
### within each omic feature
# Define the Okabe-Ito color-blind friendly palette
okabe_ito_colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Step 1: Calculate the range of Member_counts for each Genotype
genotype_order <- df101 %>%
  group_by(Genotype) %>%
  summarise(
    range_count = max(Member_counts, na.rm = TRUE) - min(Member_counts, na.rm = TRUE)
  ) %>%
  arrange(range_count) %>%
  pull(Genotype)

# Step 2: Reorder Genotype factor levels based on calculated order
df101$Genotype <- factor(df101$Genotype, levels = genotype_order)

# Step 3: Create the plot
plot <- ggplot(df101, aes(x = Genotype, y = Member_counts, color = Genotype)) +
  geom_boxplot(outlier.shape = NA) +                               # Boxplot without outlier shapes
  geom_jitter(width = 0.2, size = 2, alpha = 0.7) +                # Jitter for individual data points
  scale_color_manual(values = okabe_ito_colors) +                 # Apply the color-blind friendly palette
  labs(
    title = "Distribution of Signature Element Counts Across Omic Features",
    x = "",                                # Restore x-axis title
    y = "Signature Element Counts (Absolute)",              # Retain y-axis label
    color = "Omic Feature"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 24),             # Centered and enlarged title
    axis.title.x = element_text(size = 18, margin = margin(t = 10)), # Restore x-axis title with spacing
    axis.title.y = element_text(size = 18),                       # Enlarged y-axis title
    axis.text.x = element_blank(),                                # Hide default x-axis text
    axis.text.y = element_text(size = 12),                        # Enlarged y-axis text
    legend.title = element_text(size = 14),                       # Enlarged legend title
    legend.text = element_text(size = 12),                        # Enlarged legend text
    legend.position = c(0.05, 0.9),                               # Place legend inside plot (top-left corner)
    legend.justification = c(0, 1),                               # Adjust legend anchor point
    panel.grid.major = element_line(size = 0.5, color = "grey80"),# Subtle major gridlines
    panel.grid.minor = element_blank(),                           # No minor gridlines for clarity
    panel.background = element_rect(fill = "white", color = NA),  # Ensure white background
    plot.background = element_rect(fill = "white", color = NA)    # Ensure white plot area
  ) +
  coord_cartesian(clip = "off")                                   # Allow annotations outside the plot

print(plot)

dev.off()

# Step 4: Add genomic feature labels at y = -200
genomic_features <- data.frame(
  Genotype = levels(df101$Genotype),
  x = 1:length(levels(df101$Genotype)),
  y = -200  # Position at -200 relative to the y-axis
)

plot <- plot +
  geom_text(
    data = genomic_features,
    aes(x = x, y = y, label = Genotype),                          # Add feature labels
    inherit.aes = FALSE,                                          # Prevent overriding plot aesthetics
    size = 5,                                                     # Adjust text size
    hjust = 0.5                                                   # Center align the labels
  )

print(plot)

dev.off()

# Step 5: Add brackets **below the x-axis and above the labels**
bracket_positions <- data.frame(
  x = 1:length(levels(df101$Genotype)),
  xmin = 1:length(levels(df101$Genotype)) - 0.4,
  xmax = 1:length(levels(df101$Genotype)) + 0.4,
  y = -100  # Bracket position remains at -100
)

for (i in seq_len(nrow(bracket_positions))) {
  plot <- plot +
    # Horizontal line for the bracket
    annotate(
      "segment",
      x = bracket_positions$xmin[i], xend = bracket_positions$xmax[i],
      y = bracket_positions$y[i], yend = bracket_positions$y[i], # Horizontal line
      colour = "black", size = 0.8  # Thinner bracket
    ) +
    # Left vertical line of the bracket
    annotate(
      "segment",
      x = bracket_positions$xmin[i], xend = bracket_positions$xmin[i],
      y = bracket_positions$y[i], yend = bracket_positions$y[i] + 25, # Vertical segment
      colour = "black", size = 0.8  # Thinner bracket
    ) +
    # Right vertical line of the bracket
    annotate(
      "segment",
      x = bracket_positions$xmax[i], xend = bracket_positions$xmax[i],
      y = bracket_positions$y[i], yend = bracket_positions$y[i] + 25, # Vertical segment
      colour = "black", size = 0.8  # Thinner bracket
    )
}

# Display the plot
print(plot)

dev.off()

# Step 6: Save the plot as .pdf and .tiff
ggsave("Figure S2_Member_Counts_Per_Genotype_With_Features_XTitle.pdf", plot = plot, device = "pdf", width = 11.69, height = 8.27, dpi = 600)
ggsave("Figure S2_Member_Counts_Per_Genotype_With_Features_XTitle.tiff", plot = plot, device = "tiff", width = 11.69, height = 8.27, dpi = 600, bg = "white")

#####
#####
#####
#####
# Filter for rows where Signature == "CDH1" and TNC == "3" and create df_CDH1_3N
df_CDH1_3N <- df101 %>%
  filter(Signature == "CDH1", TNC == "3",
         SCS == "N")

# Calculate the count and percentage for CDH1 entries with TNC == "3"
df101_CDH1_3N_solo <- df_CDH1_3N %>%
  summarise(
    count_CDH1_3N_solo = n(),
    percentage_CDH1_3N_solo = (count_CDH1_3N_solo / nrow(df101)) * 100
  )

# Display the results
df_CDH1_3N
df101_CDH1_3N_solo

######
######
######
######
# "Filtering Rows Containing Target Gene with Flexible Boundaries in Complex Expressions"
# Filter for rows where the Signature column contains "SLC7A11" either standalone or within parentheses/plus signs
df_SLC7A11 <- df101 %>%
  filter(grepl("(^|\\(|\\+| )SLC7A11(\\)|\\+| |$)", Signature))

#########
#####
#####
df101_final <- df101

rio::export(df101_final, "df101_final.tsv")

# Get unique values in the "Genotype" column
genotype_values <- unique(df101_final$Genotype)

# Create separate dataframes for each unique value in "Genotype" and assign them to the environment
for (value in genotype_values) {
  # Filter df101_final based on the current genotype value
  temp_df <- df101_final %>% filter(Genotype == value)
  
  # Create a unique name for each dataframe based on the Genotype value
  assign(paste0("df101_final_", value), temp_df, envir = .GlobalEnv)
}

# Check the created dataframes in the environment
ls(pattern = "df101_final_")

#####
##### Observation: 91 out of 150 widely recognized immunological targets in cancer research are linked to RCD
datasets1k <- import("Dataset S1K.xlsx")
dfTarget_genes <- import("Target_genes.csv")

# Use trimws to remove any leading or trailing spaces in both dataframes
datasets1k$`Gene Symbol` <- trimws(datasets1k$`Gene Symbol`)
dfTarget_genes$Gene <- trimws(dfTarget_genes$Gene)

# Create a new column in dataset1k to indicate membership in RCD signature database
datasets1k$`Member of RCD signature database` <- ifelse(
  datasets1k$`Gene Symbol` %in% dfTarget_genes$Gene, "Yes", "No"
)

# Display the result
datasets1k

# Count the occurrences of "Yes" and "No" in the new column
counts <- table(datasets1k$`Member of RCD signature database`)

# Calculate the percentages
percentages <- prop.table(counts) * 100

# Combine counts and percentages in a data frame for easy display
results <- data.frame(
  Membership = names(counts),
  Count = as.integer(counts),
  Percentage = round(as.numeric(percentages), 2)
)

# Display the results
results
rio::export(datasets1k, "Dataset S1K corrected.xlsx")

#####
#####
##### plotting distribuition of signatures per genomic feature
# Load necessary libraries
library(gganimate)
library(png)

# Create the dataframe
data <- data.frame(
  Genotype = c("Transcript", "Mutation", "mRNA", "Methylation", "CNV", "miRNA", "Protein"),
  Min_Member_counts = c(1, 1, 1, 1, 1, 1, 1),
  Max_Member_counts = c(2052, 487, 477, 423, 124, 58, 4)
)

# Transform data for ggplot
data_long <- data %>%
  pivot_longer(cols = c(Min_Member_counts, Max_Member_counts), 
               names_to = "Count_Type", values_to = "Counts")

# Create the circular bar plot
ggplot(data_long, aes(x = Genotype, y = Counts, fill = Count_Type)) +
  geom_bar(stat = "identity", position = "identity", alpha = 0.7) +
  scale_fill_manual(values = c("Min_Member_counts" = "skyblue", "Max_Member_counts" = "steelblue")) +
  coord_polar(start = 0) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 10, face = "bold"),
    legend.position = "bottom"
  ) +
  labs(
    title = "Circular Bar Plot of Member Counts by Genotype",
    fill = "Count Type"
  )

dev.off()

# Create the data frame
data <- data.frame(
  Genotype = c("Transcript", "Mutation", "mRNA", "Methylation", "CNV", "miRNA", "Protein"),
  Min_Member_counts = c(1, 1, 1, 1, 1, 1, 1),
  Max_Member_counts = c(2052, 487, 477, 423, 124, 58, 4)
)

# Transform data for ggplot
data_long <- data %>%
  pivot_longer(cols = c(Min_Member_counts, Max_Member_counts), 
               names_to = "Count_Type", values_to = "Counts")

# Create the connected bubble plot
ggplot(data_long, aes(x = Genotype, y = Counts, color = Count_Type, group = Genotype)) +
  geom_line(aes(group = Genotype), color = "grey", size = 0.8) +  # Connecting line
  geom_point(aes(size = Counts), alpha = 0.7) +  # Bubbles
  scale_size(range = c(3, 15)) +  # Adjust bubble size range for better visualization
  scale_color_manual(values = c("Min_Member_counts" = "skyblue", "Max_Member_counts" = "steelblue")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 8),
    legend.position = "bottom"
  ) +
  labs(
    title = "Connected Bubble Plot of Member Counts by Genotype",
    x = "Genotype",
    y = "Member Counts",
    color = "Count Type",
    size = "Counts"
  )

dev.off()

# Load necessary libraries
library(plotly)
library(dplyr)

# Create the data frame
data <- data.frame(
  Genotype = c("Transcript", "Mutation", "mRNA", "Methylation", "CNV", "miRNA", "Protein"),
  Min_Member_counts = c(1, 1, 1, 1, 1, 1, 1),
  Max_Member_counts = c(2052, 487, 477, 423, 124, 58, 4)
)

# Transform data for plotly
data_long <- data %>%
  pivot_longer(cols = c(Min_Member_counts, Max_Member_counts), 
               names_to = "Count_Type", values_to = "Counts")

# Create the interactive plot with plotly
fig <- plot_ly(data_long, x = ~Genotype, y = ~Counts, type = 'scatter', mode = 'lines+markers',
               color = ~Count_Type, colors = c("skyblue", "steelblue"),
               text = ~paste("Genotype:", Genotype, "<br>Count Type:", Count_Type, "<br>Counts:", Counts),
               hoverinfo = "text") %>%
  layout(
    title = "Interactive Range Slider Plot of Member Counts by Genotype",
    xaxis = list(title = "Genotype"),
    yaxis = list(title = "Member Counts", rangeselector = list(
      buttons = list(
        list(count = 500, label = "500", step = "value", stepmode = "backward"),
        list(count = 1000, label = "1000", step = "value", stepmode = "backward"),
        list(count = 1500, label = "1500", step = "value", stepmode = "backward"),
        list(step = "all", label = "All")
      ))),
    showlegend = TRUE
  )

# Add range slider to the y-axis
fig <- fig %>% layout(yaxis = list(rangeslider = list(visible = TRUE)))

# Display the plot
fig

dev.off()

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(gganimate)

# Create the data frame
data <- data.frame(
  Genotype = c("Transcript", "Mutation", "mRNA", "Methylation", "CNV", "miRNA", "Protein"),
  Min_Member_counts = c(1, 1, 1, 1, 1, 1, 1),
  Max_Member_counts = c(2052, 487, 477, 423, 124, 58, 4)
)

# Transform data for ggplot
data_long <- data %>%
  pivot_longer(cols = c(Min_Member_counts, Max_Member_counts), 
               names_to = "Count_Type", values_to = "Counts")

# Create the animated plot
p <- ggplot(data_long, aes(x = Genotype, y = Counts, group = Genotype, color = Count_Type)) +
  geom_line(color = "grey", size = 1) +  # Line connecting Min and Max for each Genotype
  geom_point(aes(size = Counts), alpha = 0.7) +  # Points for Min and Max values
  scale_size(range = c(3, 10)) +  # Adjust bubble size
  scale_color_manual(values = c("Min_Member_counts" = "skyblue", "Max_Member_counts" = "steelblue")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 8),
    legend.position = "bottom"
  ) +
  labs(
    title = "Animated Connected Scatter Plot of Member Counts by Genotype",
    x = "Genotype",
    y = "Member Counts",
    color = "Count Type",
    size = "Counts"
  ) +
  transition_states(Count_Type, transition_length = 2, state_length = 1) +  # Animation transition
  enter_fade() +  # Fade-in effect
  exit_fade()     # Fade-out effect

# Save the animation as a GIF file
anim_save("connected_scatter_animation.gif", animation = animate(p, duration = 5, fps = 10, width = 800, height = 600))

dev.off()

# The GIF will be saved in your working directory as "connected_scatter_animation.gif".
# You can open it directly from your file system.

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(gganimate)

# Create the data frame
data <- data.frame(
  Genotype = c("Transcript", "Mutation", "mRNA", "Methylation", "CNV", "miRNA", "Protein"),
  Min_Member_counts = c(1, 1, 1, 1, 1, 1, 1),
  Max_Member_counts = c(2052, 487, 477, 423, 124, 58, 4)
)

# Transform data for ggplot
data_long <- data %>%
  pivot_longer(cols = c(Min_Member_counts, Max_Member_counts), 
               names_to = "Count_Type", values_to = "Counts")

# Define Okabe-Ito color palette
okabe_ito_colors <- c("Min_Member_counts" = "#E69F00", "Max_Member_counts" = "#56B4E9")

# Create the animated plot with larger bubble sizes and colorblind-friendly colors
p <- ggplot(data_long, aes(x = Genotype, y = Counts, group = Genotype, color = Count_Type)) +
  geom_line(color = "grey", size = 1) +  # Line connecting Min and Max for each Genotype
  geom_point(aes(size = Counts), alpha = 0.7) +  # Points for Min and Max values
  scale_size(range = c(9, 30)) +  # Increased bubble size range by 3x
  scale_color_manual(values = okabe_ito_colors) +  # Apply Okabe-Ito colorblind-friendly colors
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 14),
    legend.position = "bottom"
  ) +
  labs(
    title = "Animated Connected Scatter Plot of Member Counts by Genotype",
    x = "Genotype",
    y = "Member Counts",
    color = "Count Type",
    size = "Counts"
  ) +
  transition_states(Count_Type, transition_length = 2, state_length = 1) +  # Animation transition
  enter_fade() +  # Fade-in effect
  exit_fade()     # Fade-out effect

# Save the animation as a GIF file with colorblind-friendly colors and larger bubbles
anim_save("connected_scatter_animation_colorblind_friendly.gif", animation = animate(p, duration = 5, fps = 10, width = 800, height = 600))

dev.off()

######
######
library(openxlsx)

# Assuming df101_final_CNV is already loaded with dimensions 2442 x 58
df <- df101_final_CNV  # Replace this with your actual dataframe

# Define the chunk size limit for "Signature"
chunk_limit <- 30000

# Check for values in "Signature" that exceed the limit and count them
long_values_count <- sum(nchar(df$Signature) > chunk_limit)

  # Create a new workbook
  wb <- createWorkbook()

# Create the main sheet and write headers
addWorksheet(wb, "Main_Data")
writeData(wb, "Main_Data", df[0, ], startRow = 1, colNames = TRUE)  # Write only headers initially

if (long_values_count > 0) {
  # Notify about the count of values that exceed the character limit
  message(paste("There are", long_values_count, "values in the 'Signature' column exceeding", chunk_limit, "characters."))
  
  # Loop through each row to check the "Signature" column
  for (i in 1:nrow(df)) {
    if (nchar(df$Signature[i]) > chunk_limit) {
      # Split the "Signature" value into chunks of 30,000 characters to be safer
      signature_chunks <- substring(df$Signature[i], 
                                    seq(1, nchar(df$Signature[i]), 30000), 
                                    seq(30000, nchar(df$Signature[i]), 30000))
      
      # Create a separate sheet for this row's chunks
      chunk_sheet_name <- paste0("Signature_", i)
      addWorksheet(wb, chunk_sheet_name)
      
      # Write headers in the new chunk sheet
      writeData(wb, chunk_sheet_name, df[0, ], startRow = 1, colNames = TRUE)
      
      # Write the chunked "Signature" values in the new sheet, starting from the second row
      for (j in seq_along(signature_chunks)) {
        writeData(wb, chunk_sheet_name, signature_chunks[j], startRow = j + 1, startCol = which(names(df) == "Signature"))
      }
      
      # Write a text reference in the main sheet pointing to the chunk sheet
      reference_text <- paste("See", chunk_sheet_name)
      writeData(wb, "Main_Data", df[i, -which(names(df) == "Signature")], startRow = i + 1, colNames = FALSE)  # Write other columns
      writeData(wb, "Main_Data", reference_text, startRow = i + 1, startCol = which(names(df) == "Signature"))
      
    } else {
      # Write rows directly to main sheet if "Signature" is within the limit
      writeData(wb, "Main_Data", df[i, ], startRow = i + 1, colNames = FALSE)
    }
  }
  
} else {
  # Notify if no values exceed the limit
  message("No values in the 'Signature' column exceed the character limit.")
  
  # Write the entire dataframe directly to the main sheet
  writeData(wb, "Main_Data", df, startRow = 1, colNames = TRUE)
}

# Save the workbook regardless of conditions
saveWorkbook(wb, "df101_final_CNV_safe.xlsx", overwrite = TRUE)

###
library(openxlsx)

# Assuming df101_final_Methylation is already loaded with its respective dimensions
df <- df101_final_Methylation  # Use the target dataframe

# Define the chunk size limit for "Signature"
chunk_limit <- 30000

# Check for values in "Signature" that exceed the limit and count them
long_values_count <- sum(nchar(df$Signature) > chunk_limit)

# Create a new workbook
wb <- createWorkbook()

# Create the main sheet and write headers
addWorksheet(wb, "Main_Data")
writeData(wb, "Main_Data", df[0, ], startRow = 1, colNames = TRUE)  # Write only headers initially

if (long_values_count > 0) {
  # Notify about the count of values that exceed the character limit
  message(paste("There are", long_values_count, "values in the 'Signature' column exceeding", chunk_limit, "characters."))
  
  # Loop through each row to check the "Signature" column
  for (i in 1:nrow(df)) {
    if (nchar(df$Signature[i]) > chunk_limit) {
      # Split the "Signature" value into chunks of 30,000 characters
      signature_chunks <- substring(df$Signature[i], 
                                    seq(1, nchar(df$Signature[i]), chunk_limit), 
                                    seq(chunk_limit, nchar(df$Signature[i]), chunk_limit))
      
      # Create a separate sheet for this row's chunks
      chunk_sheet_name <- paste0("Signature_", i)
      addWorksheet(wb, chunk_sheet_name)
      
      # Write headers in the new chunk sheet
      writeData(wb, chunk_sheet_name, df[0, ], startRow = 1, colNames = TRUE)
      
      # Write the chunked "Signature" values in the new sheet, starting from the second row
      for (j in seq_along(signature_chunks)) {
        writeData(wb, chunk_sheet_name, signature_chunks[j], startRow = j + 1, startCol = which(names(df) == "Signature"))
      }
      
      # Write a text reference in the main sheet pointing to the chunk sheet
      reference_text <- paste("See", chunk_sheet_name)
      writeData(wb, "Main_Data", df[i, -which(names(df) == "Signature")], startRow = i + 1, colNames = FALSE)  # Write other columns
      writeData(wb, "Main_Data", reference_text, startRow = i + 1, startCol = which(names(df) == "Signature"))
      
    } else {
      # Write rows directly to main sheet if "Signature" is within the limit
      writeData(wb, "Main_Data", df[i, ], startRow = i + 1, colNames = FALSE)
    }
  }
  
} else {
  # Notify if no values exceed the limit
  message("No values in the 'Signature' column exceed the character limit.")
  
  # Write the entire dataframe directly to the main sheet
  writeData(wb, "Main_Data", df, startRow = 1, colNames = TRUE)
}

# Save the workbook regardless of conditions
saveWorkbook(wb, "df101_final_Methylation_safe.xlsx", overwrite = TRUE)

##
library(openxlsx)

# Assuming df101_final_miRNA is already loaded with its respective dimensions
df <- df101_final_miRNA  # Use the target dataframe

# Define the chunk size limit for "Signature"
chunk_limit <- 30000

# Check for values in "Signature" that exceed the limit and count them
long_values_count <- sum(nchar(df$Signature) > chunk_limit)

# Create a new workbook
wb <- createWorkbook()

# Create the main sheet and write headers
addWorksheet(wb, "Main_Data")
writeData(wb, "Main_Data", df[0, ], startRow = 1, colNames = TRUE)  # Write only headers initially

if (long_values_count > 0) {
  # Notify about the count of values that exceed the character limit
  message(paste("There are", long_values_count, "values in the 'Signature' column exceeding", chunk_limit, "characters."))
  
  # Loop through each row to check the "Signature" column
  for (i in 1:nrow(df)) {
    if (nchar(df$Signature[i]) > chunk_limit) {
      # Split the "Signature" value into chunks of 30,000 characters
      signature_chunks <- substring(df$Signature[i], 
                                    seq(1, nchar(df$Signature[i]), chunk_limit), 
                                    seq(chunk_limit, nchar(df$Signature[i]), chunk_limit))
      
      # Create a separate sheet for this row's chunks
      chunk_sheet_name <- paste0("Signature_", i)
      addWorksheet(wb, chunk_sheet_name)
      
      # Write headers in the new chunk sheet
      writeData(wb, chunk_sheet_name, df[0, ], startRow = 1, colNames = TRUE)
      
      # Write the chunked "Signature" values in the new sheet, starting from the second row
      for (j in seq_along(signature_chunks)) {
        writeData(wb, chunk_sheet_name, signature_chunks[j], startRow = j + 1, startCol = which(names(df) == "Signature"))
      }
      
      # Write a text reference in the main sheet pointing to the chunk sheet
      reference_text <- paste("See", chunk_sheet_name)
      writeData(wb, "Main_Data", df[i, -which(names(df) == "Signature")], startRow = i + 1, colNames = FALSE)  # Write other columns
      writeData(wb, "Main_Data", reference_text, startRow = i + 1, startCol = which(names(df) == "Signature"))
      
    } else {
      # Write rows directly to main sheet if "Signature" is within the limit
      writeData(wb, "Main_Data", df[i, ], startRow = i + 1, colNames = FALSE)
    }
  }
  
} else {
  # Notify if no values exceed the limit
  message("No values in the 'Signature' column exceed the character limit.")
  
  # Write the entire dataframe directly to the main sheet
  writeData(wb, "Main_Data", df, startRow = 1, colNames = TRUE)
}

# Save the workbook regardless of conditions
saveWorkbook(wb, "df101_final_miRNA_safe.xlsx", overwrite = TRUE)

##
library(openxlsx)

# Assuming df101_final_mRNA is already loaded with its respective dimensions
df <- df101_final_mRNA  # Use the target dataframe

# Define the chunk size limit for "Signature"
chunk_limit <- 30000

# Check for values in "Signature" that exceed the limit and count them
long_values_count <- sum(nchar(df$Signature) > chunk_limit)

# Create a new workbook
wb <- createWorkbook()

# Create the main sheet and write headers
addWorksheet(wb, "Main_Data")
writeData(wb, "Main_Data", df[0, ], startRow = 1, colNames = TRUE)  # Write only headers initially

if (long_values_count > 0) {
  # Notify about the count of values that exceed the character limit
  message(paste("There are", long_values_count, "values in the 'Signature' column exceeding", chunk_limit, "characters."))
  
  # Loop through each row to check the "Signature" column
  for (i in 1:nrow(df)) {
    if (nchar(df$Signature[i]) > chunk_limit) {
      # Split the "Signature" value into chunks of 30,000 characters
      signature_chunks <- substring(df$Signature[i], 
                                    seq(1, nchar(df$Signature[i]), chunk_limit), 
                                    seq(chunk_limit, nchar(df$Signature[i]), chunk_limit))
      
      # Create a separate sheet for this row's chunks
      chunk_sheet_name <- paste0("Signature_", i)
      addWorksheet(wb, chunk_sheet_name)
      
      # Write headers in the new chunk sheet
      writeData(wb, chunk_sheet_name, df[0, ], startRow = 1, colNames = TRUE)
      
      # Write the chunked "Signature" values in the new sheet, starting from the second row
      for (j in seq_along(signature_chunks)) {
        writeData(wb, chunk_sheet_name, signature_chunks[j], startRow = j + 1, startCol = which(names(df) == "Signature"))
      }
      
      # Write a text reference in the main sheet pointing to the chunk sheet
      reference_text <- paste("See", chunk_sheet_name)
      writeData(wb, "Main_Data", df[i, -which(names(df) == "Signature")], startRow = i + 1, colNames = FALSE)  # Write other columns
      writeData(wb, "Main_Data", reference_text, startRow = i + 1, startCol = which(names(df) == "Signature"))
      
    } else {
      # Write rows directly to main sheet if "Signature" is within the limit
      writeData(wb, "Main_Data", df[i, ], startRow = i + 1, colNames = FALSE)
    }
  }
  
} else {
  # Notify if no values exceed the limit
  message("No values in the 'Signature' column exceed the character limit.")
  
  # Write the entire dataframe directly to the main sheet
  writeData(wb, "Main_Data", df, startRow = 1, colNames = TRUE)
}

# Save the workbook regardless of conditions
saveWorkbook(wb, "df101_final_mRNA_safe.xlsx", overwrite = TRUE)

##
library(openxlsx)

# Assuming df101_final_Mutation is already loaded with its respective dimensions
df <- df101_final_Mutation  # Use the target dataframe

# Define the chunk size limit for "Signature"
chunk_limit <- 30000

# Check for values in "Signature" that exceed the limit and count them
long_values_count <- sum(nchar(df$Signature) > chunk_limit)

# Create a new workbook
wb <- createWorkbook()

# Create the main sheet and write headers
addWorksheet(wb, "Main_Data")
writeData(wb, "Main_Data", df[0, ], startRow = 1, colNames = TRUE)  # Write only headers initially

if (long_values_count > 0) {
  # Notify about the count of values that exceed the character limit
  message(paste("There are", long_values_count, "values in the 'Signature' column exceeding", chunk_limit, "characters."))
  
  # Loop through each row to check the "Signature" column
  for (i in 1:nrow(df)) {
    if (nchar(df$Signature[i]) > chunk_limit) {
      # Split the "Signature" value into chunks of 30,000 characters
      signature_chunks <- substring(df$Signature[i], 
                                    seq(1, nchar(df$Signature[i]), chunk_limit), 
                                    seq(chunk_limit, nchar(df$Signature[i]), chunk_limit))
      
      # Create a separate sheet for this row's chunks
      chunk_sheet_name <- paste0("Signature_", i)
      addWorksheet(wb, chunk_sheet_name)
      
      # Write headers in the new chunk sheet
      writeData(wb, chunk_sheet_name, df[0, ], startRow = 1, colNames = TRUE)
      
      # Write the chunked "Signature" values in the new sheet, starting from the second row
      for (j in seq_along(signature_chunks)) {
        writeData(wb, chunk_sheet_name, signature_chunks[j], startRow = j + 1, startCol = which(names(df) == "Signature"))
      }
      
      # Write a text reference in the main sheet pointing to the chunk sheet
      reference_text <- paste("See", chunk_sheet_name)
      writeData(wb, "Main_Data", df[i, -which(names(df) == "Signature")], startRow = i + 1, colNames = FALSE)  # Write other columns
      writeData(wb, "Main_Data", reference_text, startRow = i + 1, startCol = which(names(df) == "Signature"))
      
    } else {
      # Write rows directly to main sheet if "Signature" is within the limit
      writeData(wb, "Main_Data", df[i, ], startRow = i + 1, colNames = FALSE)
    }
  }
  
} else {
  # Notify if no values exceed the limit
  message("No values in the 'Signature' column exceed the character limit.")
  
  # Write the entire dataframe directly to the main sheet
  writeData(wb, "Main_Data", df, startRow = 1, colNames = TRUE)
}

# Save the workbook regardless of conditions
saveWorkbook(wb, "df101_final_Mutation_safe.xlsx", overwrite = TRUE)

####
library(openxlsx)

# Assuming df101_final_Protein is already loaded with its respective dimensions
df <- df101_final_Protein  # Use the target dataframe

# Define the chunk size limit for "Signature"
chunk_limit <- 30000

# Check for values in "Signature" that exceed the limit and count them
long_values_count <- sum(nchar(df$Signature) > chunk_limit)

# Create a new workbook
wb <- createWorkbook()

# Create the main sheet and write headers
addWorksheet(wb, "Main_Data")
writeData(wb, "Main_Data", df[0, ], startRow = 1, colNames = TRUE)  # Write only headers initially

if (long_values_count > 0) {
  # Notify about the count of values that exceed the character limit
  message(paste("There are", long_values_count, "values in the 'Signature' column exceeding", chunk_limit, "characters."))
  
  # Loop through each row to check the "Signature" column
  for (i in 1:nrow(df)) {
    if (nchar(df$Signature[i]) > chunk_limit) {
      # Split the "Signature" value into chunks of 30,000 characters
      signature_chunks <- substring(df$Signature[i], 
                                    seq(1, nchar(df$Signature[i]), chunk_limit), 
                                    seq(chunk_limit, nchar(df$Signature[i]), chunk_limit))
      
      # Create a separate sheet for this row's chunks
      chunk_sheet_name <- paste0("Signature_", i)
      addWorksheet(wb, chunk_sheet_name)
      
      # Write headers in the new chunk sheet
      writeData(wb, chunk_sheet_name, df[0, ], startRow = 1, colNames = TRUE)
      
      # Write the chunked "Signature" values in the new sheet, starting from the second row
      for (j in seq_along(signature_chunks)) {
        writeData(wb, chunk_sheet_name, signature_chunks[j], startRow = j + 1, startCol = which(names(df) == "Signature"))
      }
      
      # Write a text reference in the main sheet pointing to the chunk sheet
      reference_text <- paste("See", chunk_sheet_name)
      writeData(wb, "Main_Data", df[i, -which(names(df) == "Signature")], startRow = i + 1, colNames = FALSE)  # Write other columns
      writeData(wb, "Main_Data", reference_text, startRow = i + 1, startCol = which(names(df) == "Signature"))
      
    } else {
      # Write rows directly to main sheet if "Signature" is within the limit
      writeData(wb, "Main_Data", df[i, ], startRow = i + 1, colNames = FALSE)
    }
  }
  
} else {
  # Notify if no values exceed the limit
  message("No values in the 'Signature' column exceed the character limit.")
  
  # Write the entire dataframe directly to the main sheet
  writeData(wb, "Main_Data", df, startRow = 1, colNames = TRUE)
}

# Save the workbook regardless of conditions
saveWorkbook(wb, "df101_final_Protein_safe.xlsx", overwrite = TRUE)

##
##
##### Optimized Signature Data rio::export with Reduced Loop Write Operations #####

# Load required packages
library(openxlsx)

# Assuming `import` is from `rio` package

# Load data
df165 <- df101_final_Transcript

# Create a new workbook
wb <- createWorkbook()

# Add the main worksheet for primary data
addWorksheet(wb, "MainData")

# Write main data (excluding Signature) to MainData sheet in one go
writeData(wb, "MainData", df165[, -which(names(df165) == "Signature")])

# Define optimized function to split long Signature strings into chunks under Excel's character limit
split_signature <- function(signature_entry, max_length = 30000) {
  ids <- strsplit(signature_entry, " \\+ ")[[1]]
  current_chunk <- ids[1]
  chunks <- vector("list")  # Predefine the list size
  
  for (i in seq_along(ids)[-1]) {  # Start loop from the second element
    if (nchar(paste(current_chunk, ids[i], sep = " + ")) <= max_length) {
      current_chunk <- paste(current_chunk, ids[i], sep = " + ")
    } else {
      chunks[[length(chunks) + 1]] <- current_chunk
      current_chunk <- ids[i]
    }
  }
  chunks[[length(chunks) + 1]] <- current_chunk  # Final chunk
  
  return(unlist(chunks))
}

# Check if overflow is needed, and if so, prepare TranscriptDetails
overflow_needed <- any(nchar(df165$Signature) > 30000)
if (overflow_needed) {
  addWorksheet(wb, "TranscriptDetails")
  writeData(wb, "TranscriptDetails", data.frame(Row_Index = "Row_Index", Nomenclature = "Nomenclature", Signature_Chunk = "Signature_Chunk"), startRow = 1)
}

# Prepare to write data in bulk by storing rows in lists before writing them all at once
main_data_hyperlinks <- vector("list", nrow(df165))
transcript_details_data <- list()
start_row <- 2  # Row in TranscriptDetails to start data (skip headers)

# Loop over each row in df165
for (i in 1:nrow(df165)) {
  if (nchar(df165$Signature[i]) > 30000) {
    chunks <- split_signature(df165$Signature[i])
    
    # Add row information for each chunk to transcript details in a bulk list
    for (chunk in chunks) {
      transcript_details_data[[length(transcript_details_data) + 1]] <- data.frame(
        Row_Index = i + 1,
        Nomenclature = df165$Nomenclature[i],
        Signature_Chunk = chunk
      )
    }
    
    # Prepare hyperlink text for MainData
    summary_text <- paste0("n=", df165$Member_counts[i], " transcript isoforms")
    main_data_hyperlinks[[i]] <- data.frame(
      row = i + 1,
      col = which(names(df165) == "Signature"),
      text = summary_text,
      formula = paste0("HYPERLINK(\"#TranscriptDetails!A", start_row, "\", \"See full list\")")
    )
    
    start_row <- start_row + length(chunks)
  } else {
    # If within limit, write Signature entry directly to MainData
    main_data_hyperlinks[[i]] <- data.frame(
      row = i + 1,
      col = which(names(df165) == "Signature"),
      text = df165$Signature[i],
      formula = NA  # No formula required for short entries
    )
  }
}

# Bulk write transcript details data to reduce function calls
if (overflow_needed && length(transcript_details_data) > 0) {
  transcript_details_df <- do.call(rbind, transcript_details_data)
  writeData(wb, "TranscriptDetails", transcript_details_df, startRow = 2, colNames = FALSE)
}

# Write hyperlinks and text data to MainData in one go
for (item in main_data_hyperlinks) {
  writeData(wb, "MainData", item$text, startRow = item$row, startCol = item$col)
  if (!is.na(item$formula)) {
    writeFormula(wb, "MainData", item$formula, startRow = item$row, startCol = item$col)
  }
}

# Save workbook
saveWorkbook(wb, "Transcript_Optimized_Supplementary_Data_with_Selective_Hyperlinks.xlsx", overwrite = TRUE)

# Import saved workbook to df166 for verification
df166 <- import("Transcript_Optimized_Supplementary_Data_with_Selective_Hyperlinks.xlsx")

##### Optimized Signature Data  rio::export with Reduced Loop Write Operations #####

# Load required packages
library(openxlsx)
  # Assuming `import` is from `rio` package

# Load data
df165 <- df101_final_Transcript

# Relocate "Count_source" immediately after "Ranking"
df165 <- df165 %>%
  relocate(Member_counts, .after = Ranking)

# Create a new workbook
wb <- createWorkbook()

# Add the main worksheet for primary data
addWorksheet(wb, "MainData")

# Write main data (excluding Signature) to MainData sheet in one go
writeData(wb, "MainData", df165[, -which(names(df165) == "Signature")])

# Define optimized function to split long Signature strings into chunks under Excel's character limit
split_signature <- function(signature_entry, max_length = 30000) {
  ids <- strsplit(signature_entry, " \\+ ")[[1]]
  current_chunk <- ids[1]
  chunks <- vector("list")  # Predefine the list size
  
  for (i in seq_along(ids)[-1]) {  # Start loop from the second element
    if (nchar(paste(current_chunk, ids[i], sep = " + ")) <= max_length) {
      current_chunk <- paste(current_chunk, ids[i], sep = " + ")
    } else {
      chunks[[length(chunks) + 1]] <- current_chunk
      current_chunk <- ids[i]
    }
  }
  chunks[[length(chunks) + 1]] <- current_chunk  # Final chunk
  
  return(unlist(chunks))
}

# Check if overflow is needed, and if so, prepare TranscriptDetails
overflow_needed <- any(nchar(df165$Signature) > 30000)
if (overflow_needed) {
  addWorksheet(wb, "TranscriptDetails")
  writeData(wb, "TranscriptDetails", data.frame(Row_Index = "Row_Index", Nomenclature = "Nomenclature", Signature_Chunk = "Signature_Chunk"), startRow = 1)
}

# Prepare to write data in bulk by storing rows in lists before writing them all at once
main_data_hyperlinks <- vector("list", nrow(df165))
transcript_details_data <- list()
start_row <- 2  # Row in TranscriptDetails to start data (skip headers)

# Loop over each row in df165
for (i in 1:nrow(df165)) {
  if (nchar(df165$Signature[i]) > 30000) {
    chunks <- split_signature(df165$Signature[i])
    
    # Add row information for each chunk to transcript details in a bulk list
    for (chunk in chunks) {
      transcript_details_data[[length(transcript_details_data) + 1]] <- data.frame(
        Row_Index = i + 1,
        Nomenclature = df165$Nomenclature[i],
        Signature_Chunk = chunk
      )
    }
    
    # Prepare hyperlink text for MainData
    summary_text <- paste0("n=", df165$Member_counts[i], " transcript isoforms")
    main_data_hyperlinks[[i]] <- data.frame(
      row = i + 1,
      col = which(names(df165) == "Signature"),
      text = summary_text,
      formula = paste0("HYPERLINK(\"#TranscriptDetails!A", start_row, "\", \"See full list\")")
    )
    
    start_row <- start_row + length(chunks)
  } else {
    # If within limit, write Signature entry directly to MainData
    main_data_hyperlinks[[i]] <- data.frame(
      row = i + 1,
      col = which(names(df165) == "Signature"),
      text = df165$Signature[i],
      formula = NA  # No formula required for short entries
    )
  }
}

# Bulk write transcript details data to reduce function calls
if (overflow_needed && length(transcript_details_data) > 0) {
  transcript_details_df <- do.call(rbind, transcript_details_data)
  writeData(wb, "TranscriptDetails", transcript_details_df, startRow = 2, colNames = FALSE)
}

# Write hyperlinks and text data to MainData in one go
for (item in main_data_hyperlinks) {
  writeData(wb, "MainData", item$text, startRow = item$row, startCol = item$col)
  if (!is.na(item$formula)) {
    writeFormula(wb, "MainData", item$formula, startRow = item$row, startCol = item$col)
  }
}

# Save workbook
saveWorkbook(wb, "Transcript_Optimized_Supplementary_Data_with_Selective_Hyperlinks_v02.xlsx", overwrite = TRUE)

# Import saved workbook to df166 for verification
df166_v02 <- import("Transcript_Optimized_Supplementary_Data_with_Selective_Hyperlinks_v02.xlsx")

######
######
##### Optimized Signature Data rio::export with Reduced Loop Write Operations #####
######
######
######
# Load required packages
library(openxlsx)
# Assuming `import` is from `rio` package

# Load data
df165 <- df101_final_Transcript

# Relocate "Member_counts" immediately after "Ranking"
df165 <- df165 %>%
  relocate(Member_counts, .after = Ranking)

# Create a new workbook
wb <- createWorkbook()

# Add the main worksheet for primary data
addWorksheet(wb, "MainData")

# Write main data (excluding Signature) to MainData sheet in one go
writeData(wb, "MainData", df165[, -which(names(df165) == "Signature")])

# Define optimized function to split long Signature strings into chunks under Excel's character limit
split_signature <- function(signature_entry, max_length = 30000) {
  ids <- strsplit(signature_entry, " \\+ ")[[1]]
  current_chunk <- ids[1]
  chunks <- vector("list")  # Predefine the list size
  
  for (i in seq_along(ids)[-1]) {  # Start loop from the second element
    if (nchar(paste(current_chunk, ids[i], sep = " + ")) <= max_length) {
      current_chunk <- paste(current_chunk, ids[i], sep = " + ")
    } else {
      chunks[[length(chunks) + 1]] <- current_chunk
      current_chunk <- ids[i]
    }
  }
  chunks[[length(chunks) + 1]] <- current_chunk  # Final chunk
  
  return(unlist(chunks))
}

# Check if overflow is needed, and if so, prepare TranscriptDetails
overflow_needed <- any(nchar(df165$Signature) > 30000)
if (overflow_needed) {
  addWorksheet(wb, "TranscriptDetails")
  writeData(wb, "TranscriptDetails", data.frame(Row_Index = "Row_Index", Nomenclature = "Nomenclature", Signature_Chunk = "Signature_Chunk"), startRow = 1)
}

# Prepare to write data in bulk by storing rows in lists before writing them all at once
main_data_hyperlinks <- vector("list", nrow(df165))
transcript_details_data <- list()
start_row <- 2  # Row in TranscriptDetails to start data (skip headers)

# Loop over each row in df165
for (i in 1:nrow(df165)) {
  if (nchar(df165$Signature[i]) > 30000) {
    chunks <- split_signature(df165$Signature[i])
    
    # Add row information for each chunk to transcript details in a bulk list
    for (chunk in chunks) {
      transcript_details_data[[length(transcript_details_data) + 1]] <- data.frame(
        Row_Index = i + 1,
        Nomenclature = df165$Nomenclature[i],
        Signature_Chunk = chunk
      )
    }
    
    # Prepare hyperlink text for MainData
    summary_text <- paste0("n=", df165$Member_counts[i], " transcript isoforms")
    main_data_hyperlinks[[i]] <- data.frame(
      row = i + 1,
      col = which(names(df165) == "Signature"),
      text = summary_text,
      formula = paste0("HYPERLINK(\"#TranscriptDetails!A", start_row, "\", \"See full list\")")
    )
    
    start_row <- start_row + length(chunks)
  } else {
    # If within limit, write Signature entry directly to MainData
    main_data_hyperlinks[[i]] <- data.frame(
      row = i + 1,
      col = which(names(df165) == "Signature"),
      text = df165$Signature[i],
      formula = NA  # No formula required for short entries
    )
  }
}

# Bulk write transcript details data to reduce function calls
if (overflow_needed && length(transcript_details_data) > 0) {
  transcript_details_df <- do.call(rbind, transcript_details_data)
  writeData(wb, "TranscriptDetails", transcript_details_df, startRow = 2, colNames = FALSE)
}

# Write hyperlinks and text data to MainData in one go
for (item in main_data_hyperlinks) {
  writeData(wb, "MainData", item$text, startRow = item$row, startCol = item$col)
  if (!is.na(item$formula)) {
    writeFormula(wb, "MainData", item$formula, startRow = item$row, startCol = item$col)
  }
}

# Save workbook
saveWorkbook(wb, "Transcript_Optimized_Supplementary_Data_with_Selective_Hyperlinks_v02.xlsx", overwrite = TRUE)

# Import saved workbook to df166 for verification
df166_v03 <- import("Transcript_Optimized_Supplementary_Data_with_Selective_Hyperlinks_v02.xlsx")

######
######
###### top protein signatures
df101_final_Protein_TOP <- df101_final_Protein[df101_final_Protein$TMC == 1 & 
                                                 df101_final_Protein$TIC == 1, ]

# Define the correspondence table
PFC_to_Phenotype <- data.frame(
  Phenotype = c("TMB", "MSI", "TSM"),
  PFC = c(1, 2, 3)
)

# Analyzing df1011_updated

data_df1011 <- df1011_updated

# Ensure PFC in data_df1011 is numeric for proper joining
data_df1011$PFC <- as.numeric(data_df1011$PFC)

# Create the new "Phenotype" variable in df101_final based on the values in "PFC"
data_df1011<- data_df1011 %>%
  left_join(PFC_to_Phenotype, by = "PFC") %>%  # Add the "Phenotype" column based on "PFC" values
  relocate(Phenotype, .after = Genotype)       # Move "Phenotype" right after "Genotype"

# Define the correspondence table
TNC_to_Expression <- data.frame(
  Expression = c("No_data", "Unchanged", "Underexpression", "Overexpression"),  # Original TNC values
  TNC = c(0, 1, 2, 3)  # The corresponding IDs for each TNC value
)
# Ensure TNC in data_df1011 is numeric for proper joining
data_df1011$TNC <- as.numeric(data_df1011$TNC)

# Create the new "Expression" variable in data_df1011 based on the values in "TNC"
data_df1011 <- data_df1011 %>%
  left_join(TNC_to_Expression, by = "TNC") %>%  # Add the "Expression" column based on "TNC" values
  relocate(Expression, .after = Phenotype)       # Move "Expression" right after "Genotype"

# rio::export the modified dataframe
rio::export(data_df1011, "data_df1011.tsv")

### Meaningful mutation -specific signatures by immunotheray potential
meaningful_filtered_signatures_Mutation_pfc1 <- final_signatures %>%
  filter(HRC_series != "1A2A3A4A", SMC_series != "1A2A3A4A", immune_classification != "NS", 
         microenvironment_classification != "NS", 
         Genotype == "Mutation", PFC == '1',
         TNC != "1")

rio::export(meaningful_filtered_signatures_Mutation_pfc1, "Meaningful_filtered_signatures_mutation_pfc1.tsv") 

# Define the Okabe-Ito color palette
okabe_ito_colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

# Extend the Okabe-Ito color palette to match the required number of colors
extended_okabe_ito_colors <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
  "#0072B2", "#D55E00", "#CC79A7", "#999999", 
  "#CC3333" # Adding an extra distinct color
)

df210 <- meaningful_filtered_signatures_Mutation_pfc1

# Define the TMC_TIC mapping
TMC_TIC_map_2 <- data.frame(
  TMC = c("1", "2", "3"),
  TIC = c("1", "2", "3"),
  TMC_Rank = c("Anti-tumoral", "Dual", "Pro-tumoral"),
  TIC_Rank = c("Hot", "Variable", "Cold")
)

# Assuming df210 has variables CTAB, TIC, and TMC
# Map the ranks based on TMC and TIC values
df210_mapped <- df210 %>%
  mutate(
    TMC_Rank = TMC_TIC_map_2$TMC_Rank[match(TMC, TMC_TIC_map_2$TMC)],
    TIC_Rank = TMC_TIC_map_2$TIC_Rank[match(TIC, TMC_TIC_map_2$TIC)]
  )

# Combine the mapped TMC_Rank and TIC_Rank into a single variable for clarity
df210_mapped <- df210_mapped %>%
  mutate(
    Combined_Rank = paste(TMC_Rank, TIC_Rank, sep = " & ")
  )

# Reshape the data for plotting
df_long <- df210_mapped %>%
  select(CTAB, Combined_Rank) %>%
  group_by(CTAB, Combined_Rank) %>%
  summarize(Count = n(), .groups = "drop")

# Create the stacked histogram
plot <- ggplot(df_long, aes(x = CTAB, y = Count, fill = Combined_Rank)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = extended_okabe_ito_colors) +
  labs(
    #title = "Distribution of Mutation-specific signatures by immunotherapy potential",
    x = "Cancer type",
    y = "Signature counts (Absolute)",
    fill = "Combined TMC & TIC"
  ) +
  theme_minimal(base_size = 14) 

print(plot)

dev.off()

# Step 4: Save the plot as .pdf and .tiff with 600 DPI in A4 landscape orientation
ggsave("Distribution of Mutation-specific signatures by immunotherapic potential.pdf", plot = plot, device = "pdf", width = 11.69, height = 8.27, dpi = 600)
ggsave("Distribution of Mutation-specific signatures by immunotherapic potential.tiff", plot = plot, device = "tiff", width = 11.69, height = 8.27, dpi = 600, bg = "white")

rio::export(df210_mapped, "Mutation-specific signatures by immunotherapic potential.xlsx")

###
###
###
### Meaningful mutation-specific signatures by immunotherapy potential based on df1157
### Figure 8. Accumulated histogram illustrating the distribution of mutation-specific signatures
###  by meaningful immunotherapy potential across cancer types.
df1157 <- data_df1011

df1157_meaningful_filtered_signatures_Mutation_pfc1 <- df1157 %>%
  filter(HRC_series != "1A2A3A4A", SMC_series != "1A2A3A4A", immune_classification != "NS", 
         microenvironment_classification != "NS", 
         Genotype == "Mutation", PFC == '1',
         TNC != "1")

rio::export(df1157_meaningful_filtered_signatures_Mutation_pfc1, "df1157_Meaningful_filtered_signatures_mutation_pfc1.tsv") 

# Define the Okabe-Ito color palette
okabe_ito_colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

# Extend the Okabe-Ito color palette to match the required number of colors
extended_okabe_ito_colors <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
  "#0072B2", "#D55E00", "#CC79A7", "#999999", 
  "#CC3333" # Adding an extra distinct color
)


###
### 
### Validation in independent database: The Yifan_Wang_cell_death_genes example - PART A
### Which gene symbols from the 13-gene-based RCD signatures by 
### 13 gene-based cell death signature Yifan_Wang_cell_death_genes
### (SFRP1, CDO1, HGF, SETD7, IRAK3, STEAP4, CD22, C4A, VIM, TUBB6, MFN2, FOXO3, and YAP1) are in oir database
# 13 gene-based cell death signature reported in thta 15-optosis model
### (SFRP1, CDO1, HGF, SETD7, IRAK3, STEAP4, CD22, C4A, VIM, TUBB6, MFN2, FOXO3, and YAP1)
# 
# Open the PubMed link in the default browser
# Leveraging programmed cell death signature to predict clinical outcome and immunotherapy benefits 
# in postoperative bladder cancer
# shell.exec("https://www.ncbi.nlm.nih.gov/pubmed/39363008")

# Define the list of gene symbols for the 13 gene-based cell death signature
genes <- c("SFRP1", "CDO1", "HGF", "SETD7", "IRAK3", "STEAP4", 
           "CD22", "C4A", "VIM", "TUBB6", "MFN2", "FOXO3", "YAP1")

# Subset df101 to select rows containing any of the gene symbols in any column
df101_selected_rows <- df101[apply(df101, 1, function(row) any(row %in% genes)), ]

# Create a dataframe with the unique signatures of the selected rows
Yifan_Wang_cell_death_genes <- as.data.frame(unique(df101_selected_rows$Signature))

# Rename the single column to the desired name, e.g., "Gene_Symbols"
colnames(Yifan_Wang_cell_death_genes) <- "Yifan_Wang_cell_death_genes"

# rio::export the dataframe to an Excel file
rio::export(Yifan_Wang_cell_death_genes, "Yifan_Wang_cell_death_genes.xlsx")

#####
#####
#### Validation in independent database: The Yifan_Wang_cell_death_genes example - PART B
#####
##### Filtering Yifan_Wang_cell_death_genes from df1157 in BLCA
# Load necessary library
library(stringr)

# Define the list of target genes
Yifan_Wan_target_genes <- c("SFRP1", "CDO1", "HGF", "SETD7", "IRAK3", "STEAP4", 
                            "CD22", "C4A", "VIM", "TUBB6", "MFN2", "FOXO3", "YAP1")
# Create a regex pattern for whole-word matching
pattern <- paste0("\\b(", paste(Yifan_Wan_target_genes, collapse = "|"), ")\\b")

# Filter the dataframe for rows where CTAB is 'BLCA' and Signature contains any target genes
Yifan_Wan_filtered_df <- df1157[df1157$CTAB == "BLCA" & grepl(pattern, df1157$Signature, perl = TRUE), ]

# Extract the 'Signature' column from the dataframe
signature_data <- Yifan_Wan_filtered_df$Signature

# Initialize a vector to store all extracted genes
extracted_genes <- c()

# Loop through each entry in the 'Signature' column
for (entry in signature_data) {
  # Split the entry into individual genes based on delimiters (adjust as needed)
  genes <- unlist(str_split(entry, pattern = "[,;\\s]+"))
  # Add the extracted genes to the vector
  extracted_genes <- c(extracted_genes, genes)
}

# Remove any leading or trailing whitespace from gene names
extracted_genes <- str_trim(extracted_genes)

# Get unique gene names from the extracted data
unique_extracted_genes <- unique(extracted_genes)

# Determine which target genes are present in the extracted genes
present_genes <- intersect(Yifan_Wan_target_genes, unique_extracted_genes)
# [1] "SFRP1" "HGF"   "CD22"  "VIM"  
# Determine which target genes are absent in the extracted genes
absent_genes <- setdiff(Yifan_Wan_target_genes, unique_extracted_genes)
# [1] "CDO1"   "SETD7"  "IRAK3"  "STEAP4" "C4A"    "TUBB6"  "MFN2"   "FOXO3"  "YAP1" 

# Output the results
cat("Present genes in Yifan_Wan_filtered_df$Signature:\n", paste(present_genes, collapse = ", "), "\n\n")
cat("Absent genes in Yifan_Wan_filtered_df$Signature:\n", paste(absent_genes, collapse = ", "), "\n")

df211 <- df1157_meaningful_filtered_signatures_Mutation_pfc1

# Define the TMC_TIC mapping
TMC_TIC_map_2 <- data.frame(
  TMC = c("1", "2", "3"),
  TIC = c("1", "2", "3"),
  TMC_Rank = c("Anti-tumoral", "Dual", "Pro-tumoral"),
  TIC_Rank = c("Hot", "Variable", "Cold")
)

# Assuming df211 has variables CTAB, TIC, and TMC
# Map the ranks based on TMC and TIC values
df211_mapped <- df211 %>%
  mutate(
    TMC_Rank = TMC_TIC_map_2$TMC_Rank[match(TMC, TMC_TIC_map_2$TMC)],
    TIC_Rank = TMC_TIC_map_2$TIC_Rank[match(TIC, TMC_TIC_map_2$TIC)]
  )

# Combine the mapped TMC_Rank and TIC_Rank into a single variable for clarity
df211_mapped <- df211_mapped %>%
  mutate(
    Combined_Rank = paste(TMC_Rank, TIC_Rank, sep = " & ")
  )

# Reshape the data for plotting
df_long <- df211_mapped %>%
  select(CTAB, Combined_Rank) %>%
  group_by(CTAB, Combined_Rank) %>%
  summarize(Count = n(), .groups = "drop")

# Create the stacked histogram
plot <- ggplot(df_long, aes(x = CTAB, y = Count, fill = Combined_Rank)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = extended_okabe_ito_colors) +
  labs(
    #title = "Distribution of Mutation-specific signatures by immunotherapy potential",
    x = "Cancer type",
    y = "Signature counts (Absolute)",
    fill = "Combined TMC & TIC"
  ) +
  theme_minimal(base_size = 14) 

print(plot)

dev.off()

### Part C. Oversize values in "df211_mapped" 
###
###
# Define Excel's character limit
max_excel_characters <- 32767

# Identify oversized rows for each column
oversized_rows <- df211_mapped %>%
  mutate(across(where(is.character), ~ nchar(.))) %>%  # Compute character count for all text columns
  rowwise() %>%
  mutate(Row_Number = cur_group_id()) %>%  # Preserve row index
  ungroup() %>%
  pivot_longer(-Row_Number, names_to = "Column", values_to = "Char_Count") %>%
  filter(Char_Count > max_excel_characters)  # Only keep rows exceeding the limit

# Print oversized row report
if (nrow(oversized_rows) > 0) {
  print("⚠️ Warning: The following rows exceed the Excel character limit:")
  print(oversized_rows)
} else {
  print("✅ No rows exceed the Excel character limit.")
}

# Rename the selected columns in df211_mapped
df211_mapped <- df211_mapped %>%
  rename(
    Rank = Ranking,
    Elements = Count_source,
    `RCD form` = RCD_types,
    `Omic feature` = Genotype
  )

# Step 4: Save the plot as .pdf and .tiff with 600 DPI in A4 landscape orientation
ggsave("Figure 8_df1157_Distribution of Mutation-specific signatures by immunotherapic potential.pdf", plot = plot, device = "pdf", width = 11.69, height = 8.27, dpi = 600)
ggsave("Figure 8_df1157_Distribution of Mutation-specific signatures by immunotherapic potential.tiff", plot = plot, device = "tiff", width = 11.69, height = 8.27, dpi = 600, bg = "white")

rio::export(df211_mapped, "df1157_Mutation-specific signatures by immunotherapic potential.xlsx")

####
####
####
#### Making a global distribution table of signatures by genomic feature, clinical outcomes and immunotherapeutic informativeness
# Load necessary libraries
library(dplyr)
library(openxlsx)

# Define the specific order for the "Multi-omic Feature" (Genotype)
multiomic_order <- c("mRNA", "Transcript", "miRNA", "Mutation", "CNV", "Methylation", "Protein")

#########
#########
#########
#########
# Generate Table 4 (initial variables)
table4 <- df101_final %>%
  group_by(Genotype) %>% # Grouping by "Genotype" for "Multi-omic Feature"
  summarize(
    Total = n(),
    `Anti-tumoral` = sum(microenvironment_classification == "anti-tumoral"),
    `Pro-tumoral` = sum(microenvironment_classification == "pro-tumoral"),
    Dual = sum(microenvironment_classification == "dual"),
    Hot = sum(immune_classification == "Hot"),
    Cold = sum(immune_classification == "Cold"),
    Variable = sum(immune_classification == "Variable")
  ) %>%
  ungroup() %>%
  # Arrange by the specified order
  mutate(Genotype = factor(Genotype, levels = multiomic_order)) %>%
  arrange(Genotype)

# Add DSS, DFI, PFI, and OS counts
survival_metrics_counts <- df101_final %>%
  group_by(Genotype) %>%
  summarize(
    DSS_Risky = sum(Type_Cox_DSS == "Risky"),
    DSS_Protective = sum(Type_Cox_DSS == "Protective"),
    DFI_Risky = sum(Type_Cox_DFI == "Risky"),
    DFI_Protective = sum(Type_Cox_DFI == "Protective"),
    PFI_Risky = sum(Type_Cox_PFI == "Risky"),
    PFI_Protective = sum(Type_Cox_PFI == "Protective"),
    OS_Risky = sum(Type_Cox_OS == "Risky"),
    OS_Protective = sum(Type_Cox_OS == "Protective")
  )

# Add "Poor" counts
poor_counts <- df101_final %>%
  group_by(Genotype) %>%
  summarize(
    DSS_Poor = sum(Type_log_rank_DSS %in% c("MT", "High", "Low", "WT", "Duplicated", "Deleted", "Deleted/Duplicated")),
    DFI_Poor = sum(Type_log_rank_DFI %in% c("MT", "High", "Low", "WT", "Duplicated", "Deleted", "Deleted/Duplicated")),
    PFI_Poor = sum(Type_log_rank_PFI %in% c("MT", "High", "Low", "WT", "Duplicated", "Deleted", "Deleted/Duplicated")),
    OS_Poor = sum(Type_log_rank_OS %in% c("MT", "High", "Low", "WT", "Duplicated", "Deleted", "Deleted/Duplicated"))
  )

# Merge all calculated metrics into table4
table4 <- table4 %>%
  left_join(survival_metrics_counts, by = c("Genotype" = "Genotype")) %>%
  left_join(poor_counts, by = c("Genotype" = "Genotype"))

# Recalculate the "Total" row, including all new variables
total_row <- table4 %>%
  summarize(
    Genotype = "Total",
    Total = sum(Total, na.rm = TRUE),
    `Anti-tumoral` = sum(`Anti-tumoral`, na.rm = TRUE),
    `Pro-tumoral` = sum(`Pro-tumoral`, na.rm = TRUE),
    Dual = sum(Dual, na.rm = TRUE),
    Hot = sum(Hot, na.rm = TRUE),
    Cold = sum(Cold, na.rm = TRUE),
    Variable = sum(Variable, na.rm = TRUE),
    DSS_Risky = sum(DSS_Risky, na.rm = TRUE),
    DSS_Protective = sum(DSS_Protective, na.rm = TRUE),
    DFI_Risky = sum(DFI_Risky, na.rm = TRUE),
    DFI_Protective = sum(DFI_Protective, na.rm = TRUE),
    PFI_Risky = sum(PFI_Risky, na.rm = TRUE),
    PFI_Protective = sum(PFI_Protective, na.rm = TRUE),
    OS_Risky = sum(OS_Risky, na.rm = TRUE),
    OS_Protective = sum(OS_Protective, na.rm = TRUE),
    DSS_Poor = sum(DSS_Poor, na.rm = TRUE),
    DFI_Poor = sum(DFI_Poor, na.rm = TRUE),
    PFI_Poor = sum(PFI_Poor, na.rm = TRUE),
    OS_Poor = sum(OS_Poor, na.rm = TRUE)
  )

# Append the recalculated total row to the table
table4_final <- bind_rows(table4 %>% filter(Genotype != "Total"), total_row)

# Rename columns to match the desired output
colnames(table4_final) <- c(
  "Multi-omic Feature", "Total", "Anti-tumoral (TME)", "Pro-tumoral (TME)",
  "Dual (TME)", "Hot (TIC)", "Cold (TIC)", "Variable (TIC)",
  "DSS_Risky", "DSS_Protective", "DFI_Risky", "DFI_Protective",
  "PFI_Risky", "PFI_Protective", "OS_Risky", "OS_Protective",
  "DSS_Poor", "DFI_Poor", "PFI_Poor", "OS_Poor"
)

# rio::export the final table to Excel
write.xlsx(
  table4_final, 
  "Table4_Classification_Distribution_Final.xlsx", 
  rowNames = FALSE
)

# Print confirmation
print("Updated Table 4 with all metrics successfully created and rio::exported as 'Table4_Classification_Distribution_Final.xlsx'.")

# Reorder columns to place TME and TIC variables after OS_Poor
table4_final <- table4_final %>%
  select(
    `Multi-omic Feature`, `Total`, 
    `DSS_Risky`, `DSS_Protective`, `DFI_Risky`, `DFI_Protective`, 
    `PFI_Risky`, `PFI_Protective`, `OS_Risky`, `OS_Protective`, 
    `DSS_Poor`, `DFI_Poor`, `PFI_Poor`, `OS_Poor`, 
    `Anti-tumoral (TME)`, `Pro-tumoral (TME)`, `Dual (TME)`,
    `Hot (TIC)`, `Cold (TIC)`, `Variable (TIC)`
  )

# rio::export the updated table to Excel
write.xlsx(
  table4_final, 
  "Table4_Classification_Distribution_Reordered.xlsx", 
  rowNames = FALSE
)

# Print confirmation
print("Updated Table 4 with reordered variables successfully created and rio::exported as 'Table4_Classification_Distribution_Reordered.xlsx'.")

# Reorder columns to match the provided sequence
table4_final <- table4_final %>%
  select(
    `Multi-omic Feature`, `Total`,
    `DSS_Risky`, `DFI_Risky`, `PFI_Risky`, `OS_Risky`,
    `DSS_Protective`, `DFI_Protective`, `PFI_Protective`, `OS_Protective`,
    `DSS_Poor`, `DFI_Poor`, `PFI_Poor`, `OS_Poor`,
    `Anti-tumoral (TME)`, `Pro-tumoral (TME)`, `Dual (TME)`,
    `Hot (TIC)`, `Cold (TIC)`, `Variable (TIC)`
  )

# rio::export the reordered table to Excel
write.xlsx(
  table4_final, 
  "Table4_Classification_Distribution_Strict_Order.xlsx", 
  rowNames = FALSE
)

# Print confirmation
print("Table 4 reordered successfully and rio::exported as 'Table4_Classification_Distribution_Strict_Order.xlsx'.")

#########
#########
######### 
######### Table 5. Classification and Distribution of Multi-omic Signatures 
######### by Clinical Outcomes, Hazard Risk, and Therapeutic Informativeness
######### Define the specific order for the "Multi-omic Feature" (Genotype)
multiomic_order <- c("mRNA", "Transcript", "miRNA", "Mutation", "CNV", "Methylation", "Protein")

# Generate Table 5_v2 (initial variables)
table5_v2 <- df1011_updated %>%
  group_by(Genotype) %>% # Grouping by "Genotype" for "Multi-omic Feature"
  summarize(
    Total = n(),
    `Anti-tumoral` = sum(microenvironment_classification == "anti-tumoral"),
    `Pro-tumoral` = sum(microenvironment_classification == "pro-tumoral"),
    Dual = sum(microenvironment_classification == "dual"),
    Hot = sum(immune_classification == "Hot"),
    Cold = sum(immune_classification == "Cold"),
    Variable = sum(immune_classification == "Variable")
  ) %>%
  ungroup() %>%
  # Arrange by the specified order
  mutate(Genotype = factor(Genotype, levels = multiomic_order)) %>%
  arrange(Genotype)

# Add DSS, DFI, PFI, and OS counts
survival_metrics_counts <- df1011_updated %>%
  group_by(Genotype) %>%
  summarize(
    DSS_Risky = sum(Type_Cox_DSS == "Risky"),
    DSS_Protective = sum(Type_Cox_DSS == "Protective"),
    DFI_Risky = sum(Type_Cox_DFI == "Risky"),
    DFI_Protective = sum(Type_Cox_DFI == "Protective"),
    PFI_Risky = sum(Type_Cox_PFI == "Risky"),
    PFI_Protective = sum(Type_Cox_PFI == "Protective"),
    OS_Risky = sum(Type_Cox_OS == "Risky"),
    OS_Protective = sum(Type_Cox_OS == "Protective")
  )

# Add "Poor" counts
poor_counts <- df1011_updated %>%
  group_by(Genotype) %>%
  summarize(
    DSS_Poor = sum(Type_log_rank_DSS %in% c("MT", "High", "Low", "WT", "Duplicated", "Deleted", "Deleted/Duplicated")),
    DFI_Poor = sum(Type_log_rank_DFI %in% c("MT", "High", "Low", "WT", "Duplicated", "Deleted", "Deleted/Duplicated")),
    PFI_Poor = sum(Type_log_rank_PFI %in% c("MT", "High", "Low", "WT", "Duplicated", "Deleted", "Deleted/Duplicated")),
    OS_Poor = sum(Type_log_rank_OS %in% c("MT", "High", "Low", "WT", "Duplicated", "Deleted", "Deleted/Duplicated"))
  )

# Merge all calculated metrics into table5_v2
table5_v2 <- table5_v2 %>%
  left_join(survival_metrics_counts, by = c("Genotype" = "Genotype")) %>%
  left_join(poor_counts, by = c("Genotype" = "Genotype"))

# Recalculate the "Total" row, including all new variables
total_row <- table5_v2 %>%
  summarize(
    Genotype = "Total",
    Total = sum(Total, na.rm = TRUE),
    `Anti-tumoral` = sum(`Anti-tumoral`, na.rm = TRUE),
    `Pro-tumoral` = sum(`Pro-tumoral`, na.rm = TRUE),
    Dual = sum(Dual, na.rm = TRUE),
    Hot = sum(Hot, na.rm = TRUE),
    Cold = sum(Cold, na.rm = TRUE),
    Variable = sum(Variable, na.rm = TRUE),
    DSS_Risky = sum(DSS_Risky, na.rm = TRUE),
    DSS_Protective = sum(DSS_Protective, na.rm = TRUE),
    DFI_Risky = sum(DFI_Risky, na.rm = TRUE),
    DFI_Protective = sum(DFI_Protective, na.rm = TRUE),
    PFI_Risky = sum(PFI_Risky, na.rm = TRUE),
    PFI_Protective = sum(PFI_Protective, na.rm = TRUE),
    OS_Risky = sum(OS_Risky, na.rm = TRUE),
    OS_Protective = sum(OS_Protective, na.rm = TRUE),
    DSS_Poor = sum(DSS_Poor, na.rm = TRUE),
    DFI_Poor = sum(DFI_Poor, na.rm = TRUE),
    PFI_Poor = sum(PFI_Poor, na.rm = TRUE),
    OS_Poor = sum(OS_Poor, na.rm = TRUE)
  )

# Append the recalculated total row to the table
table5_v2_final <- bind_rows(table5_v2 %>% filter(Genotype != "Total"), total_row)

# Rename columns to match the desired output
colnames(table5_v2_final) <- c(
  "Multi-omic Feature", "Total", "Anti-tumoral (TME)", "Pro-tumoral (TME)",
  "Dual (TME)", "Hot (TIC)", "Cold (TIC)", "Variable (TIC)",
  "DSS_Risky", "DSS_Protective", "DFI_Risky", "DFI_Protective",
  "PFI_Risky", "PFI_Protective", "OS_Risky", "OS_Protective",
  "DSS_Poor", "DFI_Poor", "PFI_Poor", "OS_Poor"
)

# rio::export the final table to Excel
write.xlsx(
  table5_v2_final, 
  "table5_v2_Classification_Distribution_Final.xlsx", 
  rowNames = FALSE
)

# Print confirmation
print("Updated Table 5 with all metrics successfully created and rio::exported as 'table5_v2_Classification_Distribution_Final.xlsx'.")

# Reorder columns to place TME and TIC variables after OS_Poor
table5_v2_final <- table5_v2_final %>%
  select(
    `Multi-omic Feature`, `Total`, 
    `DSS_Risky`, `DSS_Protective`, `DFI_Risky`, `DFI_Protective`, 
    `PFI_Risky`, `PFI_Protective`, `OS_Risky`, `OS_Protective`, 
    `DSS_Poor`, `DFI_Poor`, `PFI_Poor`, `OS_Poor`, 
    `Anti-tumoral (TME)`, `Pro-tumoral (TME)`, `Dual (TME)`,
    `Hot (TIC)`, `Cold (TIC)`, `Variable (TIC)`
  )

# rio::export the updated table to Excel
write.xlsx(
  table5_v2_final, 
  "table5_v2_Classification_Distribution_Reordered.xlsx", 
  rowNames = FALSE
)

# Print confirmation
print("Updated Table 5 with reordered variables successfully created and rio::exported as 'table5_v2_Classification_Distribution_Reordered.xlsx'.")

# Reorder columns to match the provided sequence
table5_v2_final <- table5_v2_final %>%
  select(
    `Multi-omic Feature`, `Total`,
    `DSS_Risky`, `DFI_Risky`, `PFI_Risky`, `OS_Risky`,
    `DSS_Protective`, `DFI_Protective`, `PFI_Protective`, `OS_Protective`,
    `DSS_Poor`, `DFI_Poor`, `PFI_Poor`, `OS_Poor`,
    `Anti-tumoral (TME)`, `Pro-tumoral (TME)`, `Dual (TME)`,
    `Hot (TIC)`, `Cold (TIC)`, `Variable (TIC)`
  )

# rio::export the reordered table to Excel
rio::export(
  table5_v2_final, 
  "table5_v2_Classification_Distribution_Strict_Order.xlsx", 
  rowNames = FALSE
)

# Print confirmation
print("Table 5 reordered successfully and rio::exported as 'table5_v2_Classification_Distribution_Strict_Order.xlsx'.")

####
####
#### Genes in top signatures manuscript Tables 1, 2, 3 without duplicates for PDFAi
####
####
####
# PART A
df1148 <- df1011_select_updated_top20
df1149 <- selected_columns_df1047
df1150 <- selected_columns_combined_genotype_rows

# Rename columns
df1148 <- df1148 %>%
  rename(Rank = Ranking,
         Elements = Count_source,
         `RCD count` = RCD,
         `RCD forms` = RCD_types,
         `Multi-omic feature` = Genotype)

df1149 <- df1149 %>%
  rename(Rank = Ranking,
         Elements = Count_source,
         `RCD count` = RCD,
         `RCD forms` = RCD_types,
         `Multi-omic feature` = Genotype)

# Reorder columns
desired_order <- c("Rank", "Nomenclature", "Signature", "Elements", 
                   "Multi-omic feature", "RCD count", "RCD forms")

df1148 <- df1148[, desired_order]
df1149 <- df1149[, desired_order]

# Check if all column names and order are identical across the three dataframes
if (identical(names(df1148), names(df1149)) && identical(names(df1148), names(df1150))) {
  
  # Perform rbind operation to combine the dataframes
  df1151 <- rbind(df1148, df1149, df1150)
  
  # Print success message and show the combined dataframe
  cat("Column names and order are identical. Combined dataframe created as df1151:\n")
  print(df1151)
  
} else {
  # Print error message and show differences in column names or order
  cat("Column names or order are not identical among the dataframes.\n")
  cat("Column names in df1148:\n")
  print(names(df1148))
  cat("Column names in df1149:\n")
  print(names(df1149))
  cat("Column names in df1150:\n")
  print(names(df1150))
}

# Check for duplicate rows across all columns
duplicates_exist <- any(duplicated(df1151))

if (duplicates_exist) {
  cat("Duplicate rows found. Removing duplicates...\n")
  
  # Remove duplicate rows while retaining one instance
  df1151 <- df1151 %>%
    distinct()
  
  cat("Duplicates removed. Updated dataframe:\n")
  print(df1151)
} else {
  cat("No duplicate rows found.\n")
}

# Check for duplicate rows across all columns
duplicate_indices <- which(duplicated(df1151) | duplicated(df1151, fromLast = TRUE))  # Identify all duplicate rows

if (length(duplicate_indices) > 0) {
  cat("Duplicate rows found. Printing duplicate rows:\n")
  
  # Print the duplicated rows
  duplicated_rows <- df1151[duplicate_indices, ]
  print(duplicated_rows)
  
  # Remove duplicate rows while retaining one instance
  df1151 <- df1151 %>%
    distinct()
  
  cat("Duplicates removed. Updated dataframe:\n")
  print(df1151)
} else {
  cat("No duplicate rows found.\n")
}

#######################
#######################
##### Analysis of top signatures ####
#######################
#######################
# setwd("D:/Pré-artigo 5-optosis model/Dataframes all metrics")

data <- import("Drafts_2/gene_info_trancritos_excel2.xlsx")

data <- data %>%
  mutate_all(trimws)

miRNA_database <-  import("Drafts_2/miRNA_database.xlsx")

miRNA_database <- miRNA_database %>%
  mutate_all(trimws)

protein_info <- import("Drafts_2/protein_info.xlsx") 

protein_info <- protein_info %>%
  mutate_all(trimws)

# Create the "RCD-related" variable in the protein_info dataframe
protein_info <- protein_info %>%
  mutate(`RCD-related` = ifelse(Gene_Symbol %in% Target_genes$`Gene symbol`, "Yes", "No"))

datasetS1C <- protein_info 

# Rename the column "Protein1" to "TCGA-RPPA-pancan-clean.xena (protein composite element REF)"
datasetS1C <- datasetS1C %>%
  rename(`TCGA-RPPA-pancan-clean.xena (protein composite element REF)` = Protein1)
rio::export(datasetS1C, "Dataset S1C.xlsx")

# View the updated dataframe
head(datasetS1C)

df1151 <- df1151 %>%
  mutate_all(trimws)

c <- df1011_updated %>%
  mutate_all(trimws)

# Step 0: Filter df1011_updated based on matching values in the "Nomenclature" column of df1011_updated and "Nomenclature" column in df1151
df1011_updated_filtered <- df1011_updated %>%
  filter(Nomenclature %in% df1151$Nomenclature)

# Step 1: Copy the filtered columns "Signature" and "Nomenclature" from df1011_updated_filtered to df35
df35 <- df1011_updated_filtered %>%
  select(Signature, Nomenclature)

# Step 2: Create a new column "Gene_conversion" after "Signature"
df35 <- df35 %>%
  mutate(Gene_conversion = NA)  # Initialize the column with NA values

# Step 3: Populate "Gene_conversion" using left_join for better matching reliability

# For values starting with "ENST", join with the "data" dataframe on "Transcript_ID"
df35 <- df35 %>%
  left_join(data %>% select(Transcript_ID, Display_Name), by = c("Signature" = "Transcript_ID")) %>%
  mutate(Gene_conversion = ifelse(!is.na(Display_Name), Display_Name, Gene_conversion)) %>%
  select(-Display_Name)  # Remove the temporary column

# For values starting with "hsa-miR", join with the "miRNA_database" on "Mature1" and "Mature2"
df35 <- df35 %>%
  left_join(miRNA_database %>% select(Mature1, Gene_Symbol), by = c("Signature" = "Mature1")) %>%
  mutate(Gene_conversion = ifelse(!is.na(Gene_Symbol), Gene_Symbol, Gene_conversion)) %>%
  select(-Gene_Symbol)  # Remove the temporary column

df35 <- df35 %>%
  left_join(miRNA_database %>% select(Mature2, Gene_Symbol), by = c("Signature" = "Mature2")) %>%
  mutate(Gene_conversion = ifelse(!is.na(Gene_Symbol), Gene_Symbol, Gene_conversion)) %>%
  select(-Gene_Symbol)  # Remove the temporary column

# For protein values join with the "protein_info" on "Protein1" 
df35 <- df35 %>%
  left_join(protein_info %>% select(Protein1, Gene_Symbol), by = c("Signature" = "Protein1")) %>%
  mutate(Gene_conversion = ifelse(!is.na(Gene_Symbol), Gene_Symbol, Gene_conversion)) %>%
  select(-Gene_Symbol)  # Remove the temporary column

# Step 1: Extract gene symbols, transcript names, and miRNA from the "Signature" column
df36 <- df35 %>%
  # Use str_replace_all to remove parentheses and separate by "+" into individual elements
  mutate(Signature_cleaned = str_replace_all(Signature, "[()]", "")) %>%  # Remove parentheses
  separate_rows(Signature_cleaned, sep = " \\+ ") %>%  # Split by " + " into separate rows
  distinct() %>%  # Remove duplicate rows
  filter(str_detect(Signature_cleaned, "^ENST|^hsa-miR|^[A-Za-z]")) %>%  # Keep only gene symbols, ENST transcripts, or miRNA
  arrange(Signature_cleaned)  # Order alphabetically

# Step 2: Create the df36 with the unique and ordered entries
df36 <- df36 %>%
  select(Signature_cleaned) %>%
  rename(Extracted_Signature = Signature_cleaned)

# Step 1: Create a new column "Gene_Symbol" in df36 and populate it based on the mapping logic
df36 <- df36 %>%
  rowwise() %>%
  mutate(Gene_Symbol = case_when(
    
    # For values with prefix "ENST", look up in the "data" dataframe using "Transcript_ID"
    str_detect(Extracted_Signature, "^ENST") ~ 
      ifelse(Extracted_Signature %in% data$Transcript_ID, 
             data$Display_Name[which(data$Transcript_ID == Extracted_Signature)], 
             NA_character_),
    
    # For values with prefix "hsa-miR", look up in "miRNA_database" for Mature1
    str_detect(Extracted_Signature, "^hsa-miR") & Extracted_Signature %in% miRNA_database$Mature1 ~ 
      ifelse(Extracted_Signature %in% miRNA_database$Mature1, 
             miRNA_database$Gene_Symbol[which(miRNA_database$Mature1 == Extracted_Signature)], 
             NA_character_),
    
    # For values with prefix "hsa-miR", look up in "miRNA_database" for Mature2
    str_detect(Extracted_Signature, "^hsa-miR") & Extracted_Signature %in% miRNA_database$Mature2 ~ 
      ifelse(Extracted_Signature %in% miRNA_database$Mature2, 
             miRNA_database$Gene_Symbol[which(miRNA_database$Mature2 == Extracted_Signature)], 
             NA_character_),
    
    # If no match is found, leave as NA
    TRUE ~ NA_character_
  ))

# Fill NA values in "Gene_Symbol" with corresponding values from "Extracted_Signature"
df36 <- df36 %>%
  mutate(Gene_Symbol = ifelse(is.na(Gene_Symbol), Extracted_Signature, Gene_Symbol))

# Remove duplicates based on "Gene_Symbol"
df36 <- df36 %>%
  distinct(Gene_Symbol, .keep_all = TRUE)

# Order the dataframe by "Gene_Symbol" alphabetically
df36 <- df36 %>%
  arrange(Gene_Symbol)

# Step 1: Identify the common values between "Protein1" and "Gene_symbol"
common_proteins <- intersect(protein_info$Protein1, df36$Gene_Symbol)

# Step 2: Create a lookup table for replacements
protein_lookup <- protein_info %>%
  filter(Protein1 %in% common_proteins) %>%
  select(Protein1) %>%
  distinct()

# Step 3: Replace the "Gene_symbol" values in df36 with the corresponding "Protein1" values from protein_info
df36 <- df36 %>%
  mutate(Gene_Symbol = ifelse(Gene_Symbol %in% protein_info$Protein1,  # Check if Gene_symbol is in Protein1
                              protein_info$Gene_Symbol[match(Gene_Symbol, protein_info$Protein1)],  # Replace with corresponding Gene_symbol in protein_info
                              Gene_Symbol))  # Keep original Gene_symbol if no match

# Exclude the "Extracted_Signature" column
df36 <- df36 %>%
  select(-Extracted_Signature)

# Assign dataframe to Gene_Symbols_top_signatures_df1151
Gene_Symbols_top_signatures_df1151 <- df36

# Order the dataframe alphabetically by the "Gene_Symbol" column
Gene_Symbols_top_signatures_df1151 <- Gene_Symbols_top_signatures_df1151[order(Gene_Symbols_top_signatures_df1151$Gene_Symbol), ]

# Remove leading and trailing spaces globally across all columns
Gene_Symbols_top_signatures_df1151 <- Gene_Symbols_top_signatures_df1151 %>%
  mutate_all(trimws)

#  Remove duplicate rows
duplicates_exist <- any(duplicated(Gene_Symbols_top_signatures_df1151))

if (duplicates_exist) {
  cat("Duplicate rows found. Removing duplicates...\n")
  
  # Remove duplicate rows while retaining one instance
  Gene_Symbols_top_signatures_df1151 <- Gene_Symbols_top_signatures_df1151 %>%
    distinct()
  
  cat("Duplicates removed. Updated dataframe:\n")
  print(Gene_Symbols_top_signatures_df1151)
} else {
  cat("No duplicate rows found in the dataframe.\n")
}

# Step 5: rio::export the cleaned dataframe to a TSV file
rio::export(Gene_Symbols_top_signatures_df1151, "Gene_Symbols_top_signatures_df1151.tsv")

####
####
#### Genes in top signatures tables 1, 2, 3 without duplicates
####
####
####
# PART A
df1148 <- df1011_select_updated_top20
df1149 <- selected_columns_df1047
df1150 <- selected_columns_combined_genotype_rows

# Check if all column names and order are identical across the three dataframes
if (identical(names(df1148), names(df1149)) && identical(names(df1148), names(df1150))) {
  
  # Perform rbind operation to combine the dataframes
  df1151 <- rbind(df1148, df1149, df1150)
  
  # Print success message and show the combined dataframe
  cat("Column names and order are identical. Combined dataframe created as df1151:\n")
  print(df1151)
  
} else {
  # Print error message and show differences in column names or order
  cat("Column names or order are not identical among the dataframes.\n")
  cat("Column names in df1148:\n")
  print(names(df1148))
  cat("Column names in df1149:\n")
  print(names(df1149))
  cat("Column names in df1150:\n")
  print(names(df1150))
}

# Check for duplicate rows across all columns
duplicates_exist <- any(duplicated(df1151))

if (duplicates_exist) {
  cat("Duplicate rows found. Removing duplicates...\n")
  
  # Remove duplicate rows while retaining one instance
  df1151 <- df1151 %>%
    distinct()
  
  cat("Duplicates removed. Updated dataframe:\n")
  print(df1151)
} else {
  cat("No duplicate rows found.\n")
}

# Check for duplicate rows across all columns
duplicate_indices <- which(duplicated(df1151) | duplicated(df1151, fromLast = TRUE))  # Identify all duplicate rows

if (length(duplicate_indices) > 0) {
  cat("Duplicate rows found. Printing duplicate rows:\n")
  
  # Print the duplicated rows
  duplicated_rows <- df1151[duplicate_indices, ]
  print(duplicated_rows)
  
  # Remove duplicate rows while retaining one instance
  df1151 <- df1151 %>%
    distinct()
  
  cat("Duplicates removed. Updated dataframe:\n")
  print(df1151)
} else {
  cat("No duplicate rows found.\n")
}

#######################
#######################
##### Analysis of top signatures ####
#######################
#######################
# setwd("D:/Pré-artigo 5-optosis model/Dataframes all metrics")

data <- import("Drafts_2/gene_info_trancritos_excel2.xlsx")

data <- data %>%
  mutate_all(trimws)

miRNA_database <-  import("Drafts_2/miRNA_database.xlsx")

miRNA_database <- miRNA_database %>%
  mutate_all(trimws)

protein_info <- import("Drafts_2/protein_info.xlsx") 

protein_info <- protein_info %>%
  mutate_all(trimws)

df1151 <- df1151 %>%
  mutate_all(trimws)

df1011_updated <- df1011_updated %>%
  mutate_all(trimws)

# Step 0: Filter df1011_updated based on matching values in the "Nomenclature" column of df1011_updated and "Nomenclature" column in df1151
df1011_updated_filtered <- df1011_updated %>%
  filter(Nomenclature %in% df1151$Nomenclature)

# Step 1: Copy the filtered columns "Signature" and "Nomenclature" from df1011_updated_filtered to df35
df35 <- df1011_updated_filtered %>%
  select(Signature, Nomenclature)

# Step 2: Create a new column "Gene_conversion" after "Signature"
df35 <- df35 %>%
  mutate(Gene_conversion = NA)  # Initialize the column with NA values

# Step 3: Populate "Gene_conversion" using left_join for better matching reliability

# For values starting with "ENST", join with the "data" dataframe on "Transcript_ID"
df35 <- df35 %>%
  left_join(data %>% select(Transcript_ID, Display_Name), by = c("Signature" = "Transcript_ID")) %>%
  mutate(Gene_conversion = ifelse(!is.na(Display_Name), Display_Name, Gene_conversion)) %>%
  select(-Display_Name)  # Remove the temporary column

# For values starting with "hsa-miR", join with the "miRNA_database" on "Mature1" and "Mature2"
df35 <- df35 %>%
  left_join(miRNA_database %>% select(Mature1, Gene_Symbol), by = c("Signature" = "Mature1")) %>%
  mutate(Gene_conversion = ifelse(!is.na(Gene_Symbol), Gene_Symbol, Gene_conversion)) %>%
  select(-Gene_Symbol)  # Remove the temporary column

df35 <- df35 %>%
  left_join(miRNA_database %>% select(Mature2, Gene_Symbol), by = c("Signature" = "Mature2")) %>%
  mutate(Gene_conversion = ifelse(!is.na(Gene_Symbol), Gene_Symbol, Gene_conversion)) %>%
  select(-Gene_Symbol)  # Remove the temporary column

# For protein values join with the "protein_info" on "Protein1" 
df35 <- df35 %>%
  left_join(protein_info %>% select(Protein1, Gene_Symbol), by = c("Signature" = "Protein1")) %>%
  mutate(Gene_conversion = ifelse(!is.na(Gene_Symbol), Gene_Symbol, Gene_conversion)) %>%
  select(-Gene_Symbol)  # Remove the temporary column

# Step 1: Extract gene symbols, transcript names, and miRNA from the "Signature" column
df36 <- df35 %>%
  # Use str_replace_all to remove parentheses and separate by "+" into individual elements
  mutate(Signature_cleaned = str_replace_all(Signature, "[()]", "")) %>%  # Remove parentheses
  separate_rows(Signature_cleaned, sep = " \\+ ") %>%  # Split by " + " into separate rows
  distinct() %>%  # Remove duplicate rows
  filter(str_detect(Signature_cleaned, "^ENST|^hsa-miR|^[A-Za-z]")) %>%  # Keep only gene symbols, ENST transcripts, or miRNA
  arrange(Signature_cleaned)  # Order alphabetically

# Step 2: Create the df36 with the unique and ordered entries
df36 <- df36 %>%
  select(Signature_cleaned) %>%
  rename(Extracted_Signature = Signature_cleaned)

# Step 1: Create a new column "Gene_Symbol" in df36 and populate it based on the mapping logic
df36 <- df36 %>%
  rowwise() %>%
  mutate(Gene_Symbol = case_when(
    
    # For values with prefix "ENST", look up in the "data" dataframe using "Transcript_ID"
    str_detect(Extracted_Signature, "^ENST") ~ 
      ifelse(Extracted_Signature %in% data$Transcript_ID, 
             data$Display_Name[which(data$Transcript_ID == Extracted_Signature)], 
             NA_character_),
    
    # For values with prefix "hsa-miR", look up in "miRNA_database" for Mature1
    str_detect(Extracted_Signature, "^hsa-miR") & Extracted_Signature %in% miRNA_database$Mature1 ~ 
      ifelse(Extracted_Signature %in% miRNA_database$Mature1, 
             miRNA_database$Gene_Symbol[which(miRNA_database$Mature1 == Extracted_Signature)], 
             NA_character_),
    
    # For values with prefix "hsa-miR", look up in "miRNA_database" for Mature2
    str_detect(Extracted_Signature, "^hsa-miR") & Extracted_Signature %in% miRNA_database$Mature2 ~ 
      ifelse(Extracted_Signature %in% miRNA_database$Mature2, 
             miRNA_database$Gene_Symbol[which(miRNA_database$Mature2 == Extracted_Signature)], 
             NA_character_),
    
    # If no match is found, leave as NA
    TRUE ~ NA_character_
  ))

# Fill NA values in "Gene_Symbol" with corresponding values from "Extracted_Signature"
df36 <- df36 %>%
  mutate(Gene_Symbol = ifelse(is.na(Gene_Symbol), Extracted_Signature, Gene_Symbol))

# Remove duplicates based on "Gene_Symbol"
df36 <- df36 %>%
  distinct(Gene_Symbol, .keep_all = TRUE)

# Order the dataframe by "Gene_Symbol" alphabetically
df36 <- df36 %>%
  arrange(Gene_Symbol)


# Step 1: Identify the common values between "Protein1" and "Gene_symbol"
common_proteins <- intersect(protein_info$Protein1, df36$Gene_Symbol)

# Step 2: Create a look up table for replacements
protein_lookup <- protein_info %>%
  filter(Protein1 %in% common_proteins) %>%
  select(Protein1) %>%
  distinct()

# Step 3: Replace the "Gene_symbol" values in df36 with the corresponding "Protein1" values from protein_info
df36 <- df36 %>%
  mutate(Gene_Symbol = ifelse(Gene_Symbol %in% protein_info$Protein1,  # Check if Gene_symbol is in Protein1
                              protein_info$Gene_Symbol[match(Gene_Symbol, protein_info$Protein1)],  # Replace with corresponding Gene_symbol in protein_info
                              Gene_Symbol))  # Keep original Gene_symbol if no match

# Exclude the "Extracted_Signature" column
df36 <- df36 %>%
  select(-Extracted_Signature)

# Assign dataframe to Gene_Symbols_top_signatures_df1151
Gene_Symbols_top_signatures_df1151 <- df36

# Order the dataframe alphabetically by the "Gene_Symbol" column
Gene_Symbols_top_signatures_df1151 <- Gene_Symbols_top_signatures_df1151[order(Gene_Symbols_top_signatures_df1151$Gene_Symbol), ]

# Remove leading and trailing spaces globally across all columns
Gene_Symbols_top_signatures_df1151 <- Gene_Symbols_top_signatures_df1151 %>%
  mutate_all(trimws)

#  Remove duplicate rows
duplicates_exist <- any(duplicated(Gene_Symbols_top_signatures_df1151))

if (duplicates_exist) {
  cat("Duplicate rows found. Removing duplicates...\n")
  
  # Remove duplicate rows while retaining one instance
  Gene_Symbols_top_signatures_df1151 <- Gene_Symbols_top_signatures_df1151 %>%
    distinct()
  
  cat("Duplicates removed. Updated dataframe:\n")
  print(Gene_Symbols_top_signatures_df1151)
} else {
  cat("No duplicate rows found in the dataframe.\n")
}

#### Handling "PRKN", "PARK2" ambiguity (02/24/2025) 
#### This is so because: PARK2 is the gene Symbol that codes the Parkin protein, also known as PRKN 
# Trim whitespace, replace "PRKN" with "PARK2", and remove duplicates
df36 <- df36 %>%
  mutate(Gene_Symbol = trimws(Gene_Symbol)) %>%  # Trim whitespace
  mutate(Gene_Symbol = ifelse(Gene_Symbol == "PRKN", "PARK2", Gene_Symbol)) %>%  # Replace PRKN
  distinct()  # Remove duplicate rows

# Print confirmation message
print("Replaced 'PRKN' with 'PARK2', removed duplicates, and trimmed whitespace in 'Gene_Symbol' column.")

# Ensure the dataframes have been loaded or created: Target_genes and df36

# Filter Target_genes based on matching values in Gene_Symbol from df36
filtered_Target_genes <- Target_genes[Target_genes$`Gene symbol` %in% df36$Gene_Symbol, ]

# Identify rows in df36 where Gene_Symbol is NOT in Target_genes$`Gene symbol`
non_matching_rows <- df36[!df36$Gene_Symbol %in% Target_genes$`Gene symbol`, ]

# View the non-matching rows
head(non_matching_rows)

# Step 5: rio::export the cleaned dataframe to a TSV file
rio::export(Gene_Symbols_top_signatures_df1151, "Gene_Symbols_top_signatures_df1151.tsv")
rio::export(Gene_Symbols_top_signatures_df1151, "Gene_Symbols_top_signatures_df1151.xlsx")

cat("Cleaned dataframe rio::exported to 'Gene_Symbols_top_signatures_df1151.tsv'.\n")

######
######
###### MOST_CLINICALLY_meaningful_filtered_signatures_df1011_updated
######
######
# None_effect value_signatures: Select all rows where the value in "HRC_series" is "1A2A3A4A"
filtered_rows_HRC_df1011_updated <- df1011_updated %>%
  filter(HRC_series == "1A2A3A4A") %>%  # Filter rows where HRC_series is "1A2A3A4A"
  filter(TNC == "0") %>%  # Exclude rows where TNC is "0"
  filter(SMC != "0") 


filtered_rows_SMC_df1011_updated <- df1011_updated %>%
  filter(SMC_series == "1A2A3A4A")%>%  # Filter rows where HRC_series is "1A2A3A4A"
  filter(TNC == "0") %>%  # Exclude rows where TNC is "0"
  filter(HRC != "0") 

#### No effect value signatures: OR
filtered_rows_HRC_df1011_updated_SMC <- df1011_updated %>%
  filter(HRC_series == "1A2A3A4A" | SMC_series == "1A2A3A4A")

#### No effect value signatures: AND
filtered_rows_HRC_df1011_updated_SMC <- df1011_updated %>%
  filter(HRC_series == "1A2A3A4A" & SMC_series == "1A2A3A4A") # Logical AND: Both conditions must be true

#### Select all rows that are meaningful by "HRC_series" and "SMC_series"
# Select all rows where the value in "HRC_series" an  ""SMC_series) is different that "1A2A3A4A" 
filtered_rows_meaningful_df1011_updated <- df1011_updated %>%
  filter(HRC_series != "1A2A3A4A", SMC_series != "1A2A3A4A") %>%  # Exclude rows where TNC is "0"
  filter(TNC != "0") 

# Select all rows that are meaningful by the value in "HRC_series", "SMC_series", "immune_classification", and "microenvironment_classification") is different that "1A2A3A4A"
meaningful_filtered_signatures_df1011_updated <- df1011_updated %>%
  filter(HRC_series != "1A2A3A4A", SMC_series != "1A2A3A4A", immune_classification != "NS", microenvironment_classification != "NS")  %>%  # Exclude rows where TNC is "0"
  filter(TNC != "0") 

# Select all rows that are meaningful by the value in "HRC_series", "SMC_series", "immune_classification", and "microenvironment_classification") is different that "1A2A3A4A"
meaningful_filtered_signatures_df1011_updated <- df1011_updated %>%
  filter(HRC_series != "1A2A3A4A", SMC_series != "1A2A3A4A", immune_classification != "NS", microenvironment_classification != "NS",
         Type_log_rank_DSS != "NS", Type_log_rank_DFI != "NS", Type_log_rank_PFI != "NS", Type_log_rank_OS != "NS",
         Type_Cox_DSS != "NS", Type_Cox_DFI != "NS", Type_Cox_PFI != "NS", Type_Cox_OS != "NS" ) %>%  # Exclude rows where TNC is "0"
  filter(TNC != "0") 

#####
#####
##### TOP MOST CLINICALLY MEANINGFUL SIGNATURES
##### Table 6. Top most clinically meaningful signatures 
##### 
##### 

# Load required libraries
library(dplyr)


# Select rows that are meaningful based on specified conditions
MOST_meaningful_filtered_signatures_df1011_updated <- df1011_updated %>%
  filter(
    # Exclude rows where specific series values are "1A2A3A4A"
    HRC_series != "1A2A3A4A",
    SMC_series != "1A2A3A4A",
    
    # Exclude rows where classifications are "NS"
    immune_classification != "NS",
    microenvironment_classification != "NS",
    
    # Exclude rows where specific log rank types are "NS"
    Type_log_rank_DSS != "NS",
    Type_log_rank_DFI != "NS",
    Type_log_rank_PFI != "NS",
    Type_log_rank_OS != "NS",
    
    # Exclude rows where specific Cox types are "NS"
    Type_Cox_DSS != "NS",
    Type_Cox_DFI != "NS",
    Type_Cox_PFI != "NS",
    Type_Cox_OS != "NS",
    
    # Exclude rows where TNC is "0"
    TNC != "0"
  ) %>%
  arrange(desc(Ranking))  # Order by "Ranking" in descending order

df11052 <- MOST_meaningful_filtered_signatures_df1011_updated

# Select rows where HRC_Rank and SMC_Rank are 11, and TMC_Rank and TIC_Rank are 7
top_MOST_ranked_signatures_df11052 <- df11052 %>%
  filter(SMC_ranking == 11, TMC_ranking == 7, TIC_ranking == 7) %>%
  arrange(desc(Ranking))  # Ensure descending order

# Select specified column variables
filtered_columns_top_MOST_ranked_signatures_df11052 <- top_MOST_ranked_signatures_df11052[, c(1, 2, 4, 16, 17, 30, 33, 36, 39, 41, 43, 45, 47, 48, 50)]

# Reorder specified columns and retain the remaining columns in their original order
reordered_filtered_columns_top_MOST_ranked_signatures_df11052 <- filtered_columns_top_MOST_ranked_signatures_df11052 %>%
  select(1, 3, 2, 5, 4, everything()) %>%
  
  # Rename specified columns
  rename(
    Rank = Ranking,
    `Genomic feature` = Genotype,
    `RCD form` = RCD_types,
    TMC = microenvironment_classification,
    TIC = immune_classification
  )

# Reorder the columns in reordered_filtered_columns_top_MOST_ranked_signatures_df11052
reordered_filtered_columns_top_MOST_ranked_signatures_df11052 <- reordered_filtered_columns_top_MOST_ranked_signatures_df11052 %>%
  
  relocate(Type_log_rank_OS, .after = Type_log_rank_PFI) %>%
  relocate(Type_Cox_OS, .after = Type_Cox_PFI)


# rio::export the updated table to Excel
write.xlsx(
  reordered_filtered_columns_top_MOST_ranked_signatures_df11052, 
  "Table 6_TOP_MOST_meaningful_filtered_signatures_df11052_updated_reordered.xlsx", 
  rowNames = FALSE
)
# rio::export ordered and renamed dataframes
rio::export(MOST_meaningful_filtered_signatures_df1011_updated, "MOST_meaningful_filtered_signatures_df11052_updated.tsv") 
rio::export(MOST_meaningful_filtered_signatures_df1011_updated, "MOST_meaningful_filtered_signatures_df1052_updated.xlsx") 

# rio::export the reordered dataframe
rio::export(reordered_filtered_columns_top_MOST_ranked_signatures_df11052, "TOP_MOST_meaningful_filtered_signatures_df11052_updated.tsv") 
rio::export(reordered_filtered_columns_top_MOST_ranked_signatures_df11052, "TOP_MOST_meaningful_filtered_signatures_df11052_updated.xlsx") 

##### 
##### 
##### Table 4. Examples of Transcript-Specific Correlations of MAPK10 with Cancer Types, Phenotypic Features, and Prognostic Outcomes
##### Updated  on 16/2/2025
##### For manuscript 

# Load required library
library(openxlsx)

# Define the data with corrected variable names
data <- data.frame(
  `Signature Identifier` = c(
    "LGG-1814.5.3.P.3.93.72.2.3.2",
    "STAD-1718.5.3.N.1.44.0.3.4.2",
    "LUSC-1549.5.2.P.1.4.0.4.4.2",
    "LUAD-1824.5.1.N.1.0.0.3.4.2"
  ),
  `Transcript ID` = c(
    "ENST00000395169",
    "ENST00000395160",
    "ENST00000486985",
    "ENST00000502302"
  ),
  `Cancer Type` = c(
    "LGG (Lower-Grade Glioma)",
    "STAD (Stomach Adenocarcinoma)",
    "LUSC (Lung Squamous Cell Carcinoma)",
    "LUAD (Lung Adenocarcinoma)"
  ),
  `Phenotypic Correlation` = c(
    "Favorable outcomes",
    "Poor prognosis",
    "Microsatellite Instability (MSI)",
    "Tumor Mutational Burden (TMB)"
  ),
  `Correlation Direction` = c(
    "Protective",
    "Risk",
    "Positive",
    "Negative"
  ),
  `Comment` = c(
    "Correlated with better survival outcomes.",
    "Linked to worse survival outcomes.",
    "Transcript positively correlated with MSI phenotype.",
    "Transcript negatively correlated with high TMB, a hallmark of poor prognosis in LUAD."
  ),
  check.names = FALSE # This prevents R from altering the column names
)

# Save the data to an Excel file
write.xlsx(data, "MAPK10_Transcript_Correlations_Corrected.xlsx", rowNames = FALSE)

# Print confirmation
print("Excel file 'MAPK10_Transcript_Correlations_Corrected.xlsx' has been created successfully.")

##### 
##### 
##### Package for creating publication-ready plots with additional annotations 
##### and customization options.
##### Combining plot with Consistent Fonts and Sizes Across Plots
##### 
##### 
# Load libraries
library(patchwork)

# Shared theme
theme_shared <- theme_minimal(base_size = 14) +
  theme(plot.title = element_text(size = 16, face = "bold"))

# Create plots
p1 <- ggplot(mtcars, aes(x = wt, y = mpg)) +
  geom_point(color = "blue") +
  labs(title = "Scatter Plot: Weight vs MPG") +
  theme_shared

p2 <- ggplot(mtcars, aes(x = factor(cyl), fill = factor(cyl))) +
  geom_bar() +
  labs(title = "Bar Plot: Cylinder Counts") +
  theme_shared

p3 <- ggplot(mtcars, aes(x = mpg)) +
  geom_histogram(binwidth = 2, fill = "orange", color = "black") +
  labs(title = "Histogram: MPG Distribution") +
  theme_shared

# Combine plots
combined_plot <- (p1 | p2) / p3 +
  plot_annotation(
    title = "Combined Plots with Consistent Fonts",
    subtitle = "Example Using Patchwork",
    theme = theme_shared
  )

# Save the combined plot as PDF and TIFF in landscape orientation
ggsave(
  filename = "combined_plot_landscape.pdf",
  plot = combined_plot,
  device = "pdf",
  width = 11.69,  # A4 height in inches (landscape width)
  height = 8.27,  # A4 width in inches (landscape height)
  dpi = 600
)

ggsave(
  filename = "combined_plot_landscape.tiff",
  plot = combined_plot,
  device = "tiff",
  width = 11.69,  # A4 height in inches (landscape width)
  height = 8.27,  # A4 width in inches (landscape height)
  dpi = 600,
  compression = "lzw"  # Compression for TIFF files
)

# Print the plot for visualization
print(combined_plot)

# Print confirmation
print("Combined plot saved as PDF and TIFF in landscape orientation with A4 size at 600 DPI.")

##### 
##### IT DID NOT WORK
#### Data Preparation for Predictive Nomograms: 
#### Comprehensive Mapping and Subsetting of Genomic Signatures Dataset
##### 
##### 
##### 
# Load necessary libraries
library(dplyr)

# Load your dataset
data <- final_signatures_full_ranking

# Function to apply mappings to a column
apply_mapping <- function(column, labels) {
  mapped <- as.character(labels[as.character(column)])
  unmapped <- unique(column[is.na(mapped)])
  if (length(unmapped) > 0) {
    warning(paste("Unmapped values in column:", paste(unmapped, collapse = ", ")))
  }
  return(mapped)
}

# Define all mappings
TNC_labels <- c("1" = "Unchanged", "2" = "Underexpression", "3" = "Overexpression")
PFC_labels <- c("1" = "TMB", "2" = "MSI", "3" = "TSM")
SCS_labels <- c("P" = "positive", "N" = "negative")
GFC_labels <- c("1" = "Protein", "2" = "Mutation", "3" = "CNV", "4" = "miRNA", 
                "5" = "Transcript", "6" = "mRNA", "7" = "Methylation")
TMC_labels <- c("1" = "Anti-tumoral", "2" = "Dual", "3" = "Pro-tumoral", "4" = "NS")
TIC_labels <- c("1" = "Hot", "2" = "Variable", "3" = "Cold", "4" = "NS")

# Subset all 58 variables into `selected_data` (retain all initially for debugging)
selected_data <- data.frame(
  Ranking = data$Ranking,
  Signature = data$Signature,
  Nomenclature = data$Nomenclature,
  CTAB = data$CTAB,
  GFC = apply_mapping(data$GFC, GFC_labels),
  Genotype = factor(
    data$GFC,
    levels = c("1", "2", "3", "4", "5", "6", "7"),
    labels = c("Protein", "Mutation", "CNV", "miRNA", "Transcript", "mRNA", "Methylation")
  ),
  PFC = apply_mapping(data$PFC, PFC_labels),
  TMC = apply_mapping(data$TMC, TMC_labels),
  TIC = apply_mapping(data$TIC, TIC_labels),
  RCD = data$RCD,
  RCD_types = data$RCD_types,
  HRC = data$HRC,
  SMC = data$SMC,
  SMC_series = data$SMC_series,
  HRC_series = data$HRC_series,
  Microenvironment_Classification = data$microenvironment_classification,
  Immune_Classification_Detailed = data$immune_classification,
  SCS = apply_mapping(data$SCS, SCS_labels),
  TNC = apply_mapping(data$TNC, TNC_labels),
  Type_Log_Rank_OS = data$Type_log_rank_OS,
  Type_Log_Rank_DSS = data$Type_log_rank_DSS,
  Type_Log_Rank_DFI = data$Type_log_rank_DFI,
  Type_Log_Rank_PFI = data$Type_log_rank_PFI,
  Type_Cox_OS = data$Type_Cox_OS,
  Type_Cox_DSS = data$Type_Cox_DSS,
  Type_Cox_DFI = data$Type_Cox_DFI,
  Type_Cox_PFI = data$Type_Cox_PFI,
  LogRank_OS_ID = data$logrank_OS_id,
  LogRank_DSS_ID = data$logrank_DSS_id,
  LogRank_DFI_ID = data$logrank_DFI_id,
  LogRank_PFI_ID = data$logrank_PFI_id,
  Cox_OS_ID = data$Cox_OS_id,
  Cox_DSS_ID = data$Cox_DSS_id,
  Cox_DFI_ID = data$Cox_DFI_id,
  Cox_PFI_ID = data$Cox_PFI_id,
  CoxI_OS = data$CoxI_OS,
  CoxI_DSS = data$CoxI_DSS,
  CoxI_DFI = data$CoxI_DFI,
  CoxI_PFI = data$CoxI_PFI
)

# Validate the mappings and transformations
variables_to_validate <- c("GFC", "PFC", "TMC", "TIC", "SCS", "TNC")
for (var in variables_to_validate) {
  cat("\nUnique values in", var, ":\n")
  print(table(selected_data[[var]]))
}

# Define the desired order of the first 15 variables
desired_order <- c(
  "Ranking",
  "Signature",
  "Nomenclature",
  "CTAB",
  "GFC",
  "PFC",
  "SCS",
  "TNC",
  "HRC",
  "SMC",
  "TMC",
  "TIC",
  "RCD",
  "RCD_types",
  "Genotype"
)

# Identify the remaining variables not in the desired order
remaining_vars <- setdiff(names(selected_data), desired_order)

# Create the final order of variables: desired_order first, followed by remaining_vars
final_order <- c(desired_order, remaining_vars)

# Reorder the columns in `selected_data`
selected_data <- selected_data[, final_order]

# Verify the new column order
print(names(selected_data))

# Check for column consistency
cat("\nColumn consistency check:\n")
print(names(selected_data))

# Save the updated dataset
write.csv(selected_data, "debugged_selected_data.csv", row.names = FALSE)

cat("Dataset saved successfully: debugged_selected_data.csv\n")

#####
#####
#####
#### Data Preparation for Predictive Nomograms: 
#### Comprehensive Mapping and Subsetting of Genomic Signatures Dataset

# Load necessary libraries
library(dplyr)

# Load your dataset
data <- final_signatures_full_ranking

outcome_columns <- c(
  "Type_Cox_DSS", "Type_Cox_DFI", "Type_Cox_PFI", "Type_Cox_OS"
)
cat("Validating unique values in clinical outcome columns:\n")
for (col in outcome_columns) {
  cat("\n", col, ":\n")
  print(unique(data[[col]]))
}

# Summarize clinical outcomes
summary_outcomes <- data %>%
  select(all_of(outcome_columns)) %>%
  summarise(across(everything(), ~ paste(unique(.), collapse = ", ")))

print(summary_outcomes)

# Introduce a HRC hierarchy for categorization
data <- data %>%
  rowwise() %>%
  mutate(
    Protective_Count = sum(c_across(all_of(outcome_columns)) == "Protective"),
    Risky_Count = sum(c_across(all_of(outcome_columns)) == "Risky"),
    NS_Count = sum(c_across(all_of(outcome_columns)) == "NS"),
    Combined_Outcome = case_when(
      # Rows with four "Protective" outcomes
      Protective_Count == 4 ~ "Meaningful Protective",
      
      # Rows with 1 to 3 "Protective" outcomes and no "Risky" outcomes
      Protective_Count > 0 & Protective_Count < 4 & Risky_Count == 0 ~ "Protective",
      
      # Rows with four "Risky" outcomes
      Risky_Count == 4 ~ "Meaningful Risky",
      
      # Rows with 1 to 3 "Risky" outcomes and no "Protective" outcomes
      Risky_Count > 0 & Risky_Count < 4 & Protective_Count == 0 ~ "Risky",
      
      # Rows with four "NS" outcomes
      NS_Count == 4 ~ "Meaningless",
      
      # Default to "Mixed" for other combinations
      TRUE ~ "Mixed"
    )
  ) %>%
  ungroup()

# Relocate Combined_Outcome immediately after RCD
data <- data %>%
  relocate(Combined_Outcome, .after = RCD)

# Validate the derived categories
cat("Combined Outcome Categories:\n")
print(table(data$Combined_Outcome))

# Filter rows where Combined_Outcome is "Meaningful_Protective"
Meaningful_Protective_rows <- data %>%
  filter(Combined_Outcome == "Meaningful Protective")

# View the filtered data
print(Meaningful_Protective_rows)

# Filter rows where Combined_Outcome is "Protective"
Protective_rows <- data %>%
  filter(Combined_Outcome == "Protective")

# View the filtered data
print(Protective_rows)

# Filter rows where Combined_Outcome is "Meaningful Risky"
Meaningful_Risky_rows <- data %>%
  filter(Combined_Outcome == "Meaningful Risky")

# View the filtered data
print(Meaningful_Risky_rows)

# Filter rows where Combined_Outcome is "Risky"
Risky_rows <- data %>%
  filter(Combined_Outcome == "Risky")

# View the filtered data
print(Risky_rows)

# Filter rows where Combined_Outcome is "Mixed"
mixed_rows <- data %>%
  filter(Combined_Outcome == "Mixed")

# View the filtered data
print(mixed_rows)

# Filter rows where Combined_Outcome is "Meaningless"
Meaningless_rows <- data %>%
  filter(Combined_Outcome == "Meaningless")

# View the filtered data
print(Meaningless_rows)

# Select only the desired columns
selected_columns <- data %>%
  select(Protective_Count, Risky_Count, NS_Count)

# View the resulting dataframe
print(head(selected_columns))

# Relocate the variables immediately after "Combined_Outcome"
data <- data %>%
  relocate(Protective_Count, Risky_Count, NS_Count, .after = Combined_Outcome)

# Validate the column order
cat("Updated column order:\n")
print(names(data))

# Define the Okabe-Ito color-blind friendly palette
okabe_ito_palette <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#999999"
)

# Prepare data for plotting
plot_data <- data %>%
  count(Combined_Outcome, Genotype) %>%
  arrange(Combined_Outcome, Genotype)

# Create a stacked bar plot
ggplot(plot_data, aes(x = Combined_Outcome, y = n, fill = Genotype)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = okabe_ito_palette) +
  theme_minimal() +
  labs(
    title = "Stacked Distribution by Combined Outcome and Genotype",
    x = "Combined Outcome",
    y = "Number of Rows",
    fill = "Genotype"
  ) +
  theme(
    text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


# Exclude rows with "Meaningless" in Combined_Outcome
filtered_data <- data %>%
  filter(Combined_Outcome != "Meaningless")

# Define Okabe-Ito color palette
okabe_ito_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                       "#0072B2", "#D55E00", "#CC79A7")

# Create the plot
stacked_plot <- ggplot(filtered_data, aes(x = Genotype, fill = Combined_Outcome)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = okabe_ito_palette) +
  labs(
    title = "Stacked Distribution of Combined Outcome by Genotype",
    x = "Genotype",
    y = "Proportion",
    fill = "Combined Outcome"
  ) +
  theme_minimal(base_size = 14)

# Save the plot as TIFF
ggsave(
  filename = "stacked_distribution_plot.tiff",
  plot = stacked_plot,
  device = "tiff",
  dpi = 600,
  width = 8.27, height = 11.69 # A4 dimensions in inches
)

# Save the plot as PDF
ggsave(
  filename = "stacked_distribution_plot.pdf",
  plot = stacked_plot,
  device = "pdf",
  width = 8.27, height = 11.69 # A4 dimensions in inches
)

cat("Plot saved as 'stacked_distribution_plot.tiff' and 'stacked_distribution_plot.pdf' successfully!\n")

dev.off()

# Load necessary libraries
library(ggplot2)

# Exclude "Meaningless" from "Combined_Outcome"
filtered_data <- data %>%
  filter(Combined_Outcome != "Meaningless")

# Create the proportional distribution plot
plot <- ggplot(filtered_data, aes(x = Combined_Outcome, fill = Genotype)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    x = "Combined Outcome",
    y = "Proportional Distribution",
    title = "Proportional Distribution of Genotype by Combined Outcome",
    fill = "Genotype"
  ) +
  scale_fill_manual(values = c(
    "#E69F00", "#56B4E9", "#009E73", "#F0E442",
    "#0072B2", "#D55E00", "#CC79A7"  # Okabe-Ito palette
  )) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
  )

# Save the plot as PNG
ggsave("Proportional_Distribution_by_Combined_Outcome.png", plot = plot, 
       dpi = 600, width = 8.3, height = 11.7, units = "in")  # A4 dimensions

# Save the plot as TIFF
ggsave("Proportional_Distribution_by_Combined_Outcome.tiff", plot = plot, 
       dpi = 600, width = 8.3, height = 11.7, units = "in")  # A4 dimensions

# Print the plot to the console
print(plot)

dev.off()
#####
#####
##### snippet to create an interactive proportional distribution plot of "Genotype" by "Combined_Outcome," excluding the "Meaningless" category.
#####
##### 

# Exclude "Meaningless" from "Combined_Outcome"
filtered_data <- data %>%
  filter(Combined_Outcome != "Meaningless")

# Create the proportional distribution plot
ggplot_obj <- ggplot(filtered_data, aes(x = Combined_Outcome, fill = Genotype)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    x = "Combined Outcome",
    y = "Proportional Distribution",
    title = "Proportional Distribution of Genotype by Combined Outcome",
    fill = "Genotype"
  ) +
  scale_fill_manual(values = c(
    "#E69F00", "#56B4E9", "#009E73", "#F0E442",
    "#0072B2", "#D55E00", "#CC79A7"  # Okabe-Ito palette
  )) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
  )

# Convert ggplot to an interactive plotly object
interactive_plot <- ggplotly(ggplot_obj)

# Display the interactive plot
interactive_plot

#####
#####
#####
##### IT DID NOT WORKED
#### Data Preparation for Predictive Nomograms: 
#### Comprehensive Mapping and Subsetting of Genomic Signatures Dataset
#####
#####
# Load your dataset
data <- final_signatures_full_ranking

# Function to apply mappings to a column
apply_mapping <- function(column, labels) {
  mapped <- as.character(labels[as.character(column)])
  unmapped <- unique(column[is.na(mapped)])
  if (length(unmapped) > 0) {
    warning(paste("Unmapped values in column:", paste(unmapped, collapse = ", ")))
  }
  return(mapped)
}

# Define all mappings
TNC_labels <- c("1" = "Unchanged", "2" = "Underexpression", "3" = "Overexpression")
PFC_labels <- c("1" = "TMB", "2" = "MSI", "3" = "TSM")
SCS_labels <- c("P" = "positive", "N" = "negative")
GFC_labels <- c("1" = "Protein", "2" = "Mutation", "3" = "CNV", "4" = "miRNA", 
                "5" = "Transcript", "6" = "mRNA", "7" = "Methylation")
TMC_labels <- c("1" = "Anti-tumoral", "2" = "Dual", "3" = "Pro-tumoral", "4" = "NS")
TIC_labels <- c("1" = "Hot", "2" = "Variable", "3" = "Cold", "4" = "NS")

# Subset all 58 variables into `selected_data` (retain all initially for debugging)
selected_data <- data.frame(
  Ranking = data$Ranking,
  Signature = data$Signature,
  Nomenclature = data$Nomenclature,
  CTAB = data$CTAB,
  GFC = apply_mapping(data$GFC, GFC_labels),
  Genotype = factor(
    data$GFC,
    levels = c("1", "2", "3", "4", "5", "6", "7"),
    labels = c("Protein", "Mutation", "CNV", "miRNA", "Transcript", "mRNA", "Methylation")
  ),
  PFC = apply_mapping(data$PFC, PFC_labels),
  TMC = apply_mapping(data$TMC, TMC_labels),
  TIC = apply_mapping(data$TIC, TIC_labels),
  RCD = data$RCD,
  RCD_types = data$RCD_types,
  HRC = data$HRC,
  SMC = data$SMC,
  SMC_series = data$SMC_series,
  HRC_series = data$HRC_series,
  Microenvironment_Classification = data$microenvironment_classification,
  Immune_Classification_Detailed = data$immune_classification,
  SCS = apply_mapping(data$SCS, SCS_labels),
  TNC = apply_mapping(data$TNC, TNC_labels),
  Type_Log_Rank_OS = data$Type_log_rank_OS,
  Type_Log_Rank_DSS = data$Type_log_rank_DSS,
  Type_Log_Rank_DFI = data$Type_log_rank_DFI,
  Type_Log_Rank_PFI = data$Type_log_rank_PFI,
  Type_Cox_OS = data$Type_Cox_OS,
  Type_Cox_DSS = data$Type_Cox_DSS,
  Type_Cox_DFI = data$Type_Cox_DFI,
  Type_Cox_PFI = data$Type_Cox_PFI,
  LogRank_OS_ID = data$logrank_OS_id,
  LogRank_DSS_ID = data$logrank_DSS_id,
  LogRank_DFI_ID = data$logrank_DFI_id,
  LogRank_PFI_ID = data$logrank_PFI_id,
  Cox_OS_ID = data$Cox_OS_id,
  Cox_DSS_ID = data$Cox_DSS_id,
  Cox_DFI_ID = data$Cox_DFI_id,
  Cox_PFI_ID = data$Cox_PFI_id,
  CoxI_OS = data$CoxI_OS,
  CoxI_DSS = data$CoxI_DSS,
  CoxI_DFI = data$CoxI_DFI,
  CoxI_PFI = data$CoxI_PFI
)

# Validate the mappings and transformations
variables_to_validate <- c("GFC", "PFC", "TMC", "TIC", "SCS", "TNC")
for (var in variables_to_validate) {
  cat("\nUnique values in", var, ":\n")
  print(table(selected_data[[var]]))
}

# Check for column consistency
cat("\nColumn consistency check:\n")
print(names(selected_data))

# Save the updated dataset
write.csv(selected_data, "debugged_selected_data.csv", row.names = FALSE)

cat("Dataset saved successfully: debugged_selected_data.csv\n")

#### Table 7 Making a global distribution table of signatures by genomic feature, clinical outcomes and immunotherapeutic informativeness
# Load necessary libraries
library(dplyr)
library(openxlsx)

# Define the specific order for the "Multi-omic Feature" (Genotype)
multiomic_order <- c("mRNA", "Transcript", "miRNA", "Mutation", "CNV", "Methylation", "Protein")

# Generate Table 7(initial variables)
table7 <- df1011_updated %>%
  group_by(Genotype) %>% # Grouping by "Genotype" for "Multi-omic Feature"
  summarize(
    Total = n(),
    `Anti-tumoral` = sum(microenvironment_classification == "anti-tumoral"),
    `Pro-tumoral` = sum(microenvironment_classification == "pro-tumoral"),
    Dual = sum(microenvironment_classification == "dual"),
    Hot = sum(immune_classification == "Hot"),
    Cold = sum(immune_classification == "Cold"),
    Variable = sum(immune_classification == "Variable")
  ) %>%
  ungroup() %>%
  # Arrange by the specified order
  mutate(Genotype = factor(Genotype, levels = multiomic_order)) %>%
  arrange(Genotype)

# Add DSS, DFI, PFI, and OS counts
survival_metrics_counts <- df1011_updated %>%
  group_by(Genotype) %>%
  summarize(
    DSS_Risky = sum(Type_Cox_DSS == "Risky"),
    DSS_Protective = sum(Type_Cox_DSS == "Protective"),
    DFI_Risky = sum(Type_Cox_DFI == "Risky"),
    DFI_Protective = sum(Type_Cox_DFI == "Protective"),
    PFI_Risky = sum(Type_Cox_PFI == "Risky"),
    PFI_Protective = sum(Type_Cox_PFI == "Protective"),
    OS_Risky = sum(Type_Cox_OS == "Risky"),
    OS_Protective = sum(Type_Cox_OS == "Protective")
  )

# Add "Poor" counts
poor_counts <- df1011_updated %>%
  group_by(Genotype) %>%
  summarize(
    DSS_Poor = sum(Type_log_rank_DSS %in% c("MT", "High", "Low", "WT", "Duplicated", "Deleted", "Deleted/Duplicated")),
    DFI_Poor = sum(Type_log_rank_DFI %in% c("MT", "High", "Low", "WT", "Duplicated", "Deleted", "Deleted/Duplicated")),
    PFI_Poor = sum(Type_log_rank_PFI %in% c("MT", "High", "Low", "WT", "Duplicated", "Deleted", "Deleted/Duplicated")),
    OS_Poor = sum(Type_log_rank_OS %in% c("MT", "High", "Low", "WT", "Duplicated", "Deleted", "Deleted/Duplicated"))
  )

# Merge all calculated metrics into table7
table7 <- table7 %>%
  left_join(survival_metrics_counts, by = c("Genotype" = "Genotype")) %>%
  left_join(poor_counts, by = c("Genotype" = "Genotype"))

# Recalculate the "Total" row, including all new variables
total_row <- table7 %>%
  summarize(
    Genotype = "Total",
    Total = sum(Total, na.rm = TRUE),
    `Anti-tumoral` = sum(`Anti-tumoral`, na.rm = TRUE),
    `Pro-tumoral` = sum(`Pro-tumoral`, na.rm = TRUE),
    Dual = sum(Dual, na.rm = TRUE),
    Hot = sum(Hot, na.rm = TRUE),
    Cold = sum(Cold, na.rm = TRUE),
    Variable = sum(Variable, na.rm = TRUE),
    DSS_Risky = sum(DSS_Risky, na.rm = TRUE),
    DSS_Protective = sum(DSS_Protective, na.rm = TRUE),
    DFI_Risky = sum(DFI_Risky, na.rm = TRUE),
    DFI_Protective = sum(DFI_Protective, na.rm = TRUE),
    PFI_Risky = sum(PFI_Risky, na.rm = TRUE),
    PFI_Protective = sum(PFI_Protective, na.rm = TRUE),
    OS_Risky = sum(OS_Risky, na.rm = TRUE),
    OS_Protective = sum(OS_Protective, na.rm = TRUE),
    DSS_Poor = sum(DSS_Poor, na.rm = TRUE),
    DFI_Poor = sum(DFI_Poor, na.rm = TRUE),
    PFI_Poor = sum(PFI_Poor, na.rm = TRUE),
    OS_Poor = sum(OS_Poor, na.rm = TRUE)
  )

# Append the recalculated total row to the table
table7_final <- bind_rows(table7 %>% filter(Genotype != "Total"), total_row)

# Rename columns to match the desired output
colnames(table7_final) <- c(
  "Multi-omic Feature", "Total", "Anti-tumoral (TME)", "Pro-tumoral (TME)",
  "Dual (TME)", "Hot (TIC)", "Cold (TIC)", "Variable (TIC)",
  "DSS_Risky", "DSS_Protective", "DFI_Risky", "DFI_Protective",
  "PFI_Risky", "PFI_Protective", "OS_Risky", "OS_Protective",
  "DSS_Poor", "DFI_Poor", "PFI_Poor", "OS_Poor"
)

# rio::export the final table to Excel
write.xlsx(
  table7_final, 
  "table7_Classification_Distribution_Final.xlsx", 
  rowNames = FALSE
)

# Print confirmation
print("Updated Table 7with all metrics successfully created and rio::exported as 'table7_Classification_Distribution_Final.xlsx'.")

# Reorder columns to place TME and TIC variables after OS_Poor
table7_final <- table7_final %>%
  select(
    `Multi-omic Feature`, `Total`, 
    `DSS_Risky`, `DSS_Protective`, `DFI_Risky`, `DFI_Protective`, 
    `PFI_Risky`, `PFI_Protective`, `OS_Risky`, `OS_Protective`, 
    `DSS_Poor`, `DFI_Poor`, `PFI_Poor`, `OS_Poor`, 
    `Anti-tumoral (TME)`, `Pro-tumoral (TME)`, `Dual (TME)`,
    `Hot (TIC)`, `Cold (TIC)`, `Variable (TIC)`
  )

# rio::export the updated table to Excel
write.xlsx(
  table7_final, 
  "table7_Classification_Distribution_Reordered.xlsx", 
  rowNames = FALSE
)

# Print confirmation
print("Updated Table 7 with reordered variables successfully created and rio::exported as 'table7_Classification_Distribution_Reordered.xlsx'.")

# Reorder columns to match the provided sequence
table7_final <- table7_final %>%
  select(
    `Multi-omic Feature`, `Total`,
    `DSS_Risky`, `DFI_Risky`, `PFI_Risky`, `OS_Risky`,
    `DSS_Protective`, `DFI_Protective`, `PFI_Protective`, `OS_Protective`,
    `DSS_Poor`, `DFI_Poor`, `PFI_Poor`, `OS_Poor`,
    `Anti-tumoral (TME)`, `Pro-tumoral (TME)`, `Dual (TME)`,
    `Hot (TIC)`, `Cold (TIC)`, `Variable (TIC)`
  )

# rio::export the reordered table to Excel
write.xlsx(
  table7_final, 
  "table7_Classification_Distribution_Strict_Order.xlsx", 
  rowNames = FALSE
)

# Transpose the dataframe
table7_transposed <- as.data.frame(t(table7_final))

# Create a new column for the original row names
table7_transposed <- cbind(Row_Names = rownames(table7_transposed), table7_transposed)

# Assign the first row as column names
colnames(table7_transposed) <- table7_transposed[1, ]

# Remove the first row as it is now redundant
table7_transposed <- table7_transposed[-1, , drop = FALSE]

# Reset row names to default numeric indices
rownames(table7_transposed) <- NULL  

# rio::export the reordered table to Excel
write.xlsx(
  table7_transposed, 
  "table7_transposed_Classification_Distribution_Strict_Order.xlsx", 
  rowNames = FALSE
)

# Display the transposed dataframe
print(table7_transposed)
# Load necessary libraries
library(ggplot2)

#####
#####
#####
#### data_df1011 Preparation for Predictive Nomograms: 
#### Comprehensive Mapping and Subsetting of Genomic Signatures data_df1011 set
#####
#####
# Load necessary libraries
library(dplyr)

# Load your data_df1011set
data_df1011 <- df1011_updated

outcome_columns <- c(
  "Type_Cox_DSS", "Type_Cox_DFI", "Type_Cox_PFI", "Type_Cox_OS"
)
cat("Validating unique values in clinical outcome columns:\n")
for (col in outcome_columns) {
  cat("\n", col, ":\n")
  print(unique(data_df1011[[col]]))
}

# Summarize clinical outcomes
summary_outcomes <- data_df1011 %>%
  select(all_of(outcome_columns)) %>%
  summarise(across(everything(), ~ paste(unique(.), collapse = ", ")))

print(summary_outcomes)

# Introduce a HRC hierarchy for categorization
data_df1011 <- data_df1011 %>%
  rowwise() %>%
  mutate(
    Protective_Count = sum(c_across(all_of(outcome_columns)) == "Protective"),
    Risky_Count = sum(c_across(all_of(outcome_columns)) == "Risky"),
    NS_Count = sum(c_across(all_of(outcome_columns)) == "NS"),
    Combined_Outcome = case_when(
      # Rows with four "Protective" outcomes
      Protective_Count == 4 ~ "Meaningful Protective",
      
      # Rows with 1 to 3 "Protective" outcomes and no "Risky" outcomes
      Protective_Count > 0 & Protective_Count < 4 & Risky_Count == 0 ~ "Protective",
      
      # Rows with four "Risky" outcomes
      Risky_Count == 4 ~ "Meaningful Risky",
      
      # Rows with 1 to 3 "Risky" outcomes and no "Protective" outcomes
      Risky_Count > 0 & Risky_Count < 4 & Protective_Count == 0 ~ "Risky",
      
      # Rows with four "NS" outcomes
      NS_Count == 4 ~ "Meaningless",
      
      # Default to "Mixed" for other combinations
      TRUE ~ "Mixed"
    )
  ) %>%
  ungroup()

# Relocate Combined_Outcome immediately after RCD
data_df1011 <- data_df1011 %>%
  relocate(Combined_Outcome, .after = RCD)

# Validate the derived categories
cat("Combined Outcome Categories:\n")
print(table(data_df1011$Combined_Outcome))

# Filter rows where Combined_Outcome is "Meaningful_Protective"
Meaningful_Protective_rows <- data_df1011 %>%
  filter(Combined_Outcome == "Meaningful Protective")

# View the filtered data_df1011
print(Meaningful_Protective_rows)

# Filter rows where Combined_Outcome is "Protective"
Protective_rows <- data_df1011 %>%
  filter(Combined_Outcome == "Protective")

# View the filtered data_df1011
print(Protective_rows)

# Filter rows where Combined_Outcome is "Meaningful Risky"
Meaningful_Risky_rows <- data_df1011 %>%
  filter(Combined_Outcome == "Meaningful Risky")

# View the filtered data_df1011
print(Meaningful_Risky_rows)

# Filter rows where Combined_Outcome is "Risky"
Risky_rows <- data_df1011 %>%
  filter(Combined_Outcome == "Risky")

# View the filtered data_df1011
print(Risky_rows)

# Filter rows where Combined_Outcome is "Mixed"
mixed_rows <- data_df1011 %>%
  filter(Combined_Outcome == "Mixed")

# View the filtered data_df1011
print(mixed_rows)

# Filter rows where Combined_Outcome is "Meaningless"
Meaningless_rows <- data_df1011 %>%
  filter(Combined_Outcome == "Meaningless")

# View the filtered data_df1011
print(Meaningless_rows)

# Select only the desired columns
selected_columns <- data_df1011 %>%
  select(Protective_Count, Risky_Count, NS_Count)

# View the resulting data_df1011frame
print(head(selected_columns))

# Relocate the variables immediately after "Combined_Outcome"
data_df1011 <- data_df1011 %>%
  relocate(Protective_Count, Risky_Count, NS_Count, .after = Combined_Outcome)

# Define the correspondence table
PFC_to_Phenotype <- data.frame(
  Phenotype = c("TMB", "MSI", "TSM"),
  PFC = c(1, 2, 3)
)

# Ensure PFC in data_df1011 is numeric for proper joining
data_df1011$PFC <- as.numeric(data_df1011$PFC)

# Create the new "Phenotype" variable in data_df1011 based on the values in "PFC"
data_df1011 <- data_df1011 %>%
  left_join(PFC_to_Phenotype, by = "PFC") %>%  # Add the "Phenotype" column based on "PFC" values
  relocate(Phenotype, .after = Genotype)       # Move "Phenotype" right after "Genotype"

# rio::export the modified dataframe
rio::export(data_df1011, "data_df1011.tsv")


# Validate the column order
cat("Updated column order:\n")
print(names(data_df1011))


# Ensure the variables are converted to numeric
data_df1011$Ranking <- as.numeric(as.character(data_df1011$Ranking))
data_df1011$Count_source <- as.numeric(as.character(data_df1011$Count_source))

# Check for any NA values that might result from non-numeric conversion
if(any(is.na(data_df1011$Ranking) | is.na(data_df1011$Count_source))) {
  warning("Some values could not be converted to numeric. Check for non-numeric entries.")
}

# Display summary to confirm the conversion
summary(data_df1011[, c("Ranking", "Count_source")])

#### Comprehensive Descriptive Statistics with Centiles for 'Count_source' in 'data_df1011' (Transposed Format
# Load necessary package
library(dplyr)
library(tidyr)

# Generate full descriptive statistics including centiles for "Count_source"
count_source_stats <- data_df1011 %>%
  summarise(
    Min = min(Count_source, na.rm = TRUE),
    Max = max(Count_source, na.rm = TRUE),
    Mean = mean(Count_source, na.rm = TRUE),
    Median = median(Count_source, na.rm = TRUE),
    Std_Dev = sd(Count_source, na.rm = TRUE),
    Variance = var(Count_source, na.rm = TRUE),
    Q1 = quantile(Count_source, 0.25, na.rm = TRUE),  # 25th percentile
    Q3 = quantile(Count_source, 0.75, na.rm = TRUE),  # 75th percentile
    IQR = IQR(Count_source, na.rm = TRUE),  # Interquartile Range
    P1 = quantile(Count_source, 0.01, na.rm = TRUE),  # 1st percentile
    P5 = quantile(Count_source, 0.05, na.rm = TRUE),  # 5th percentile
    P10 = quantile(Count_source, 0.10, na.rm = TRUE), # 10th percentile
    P90 = quantile(Count_source, 0.90, na.rm = TRUE), # 90th percentile
    P95 = quantile(Count_source, 0.95, na.rm = TRUE), # 95th percentile
    P99 = quantile(Count_source, 0.99, na.rm = TRUE), # 99th percentile
    Missing_Values = sum(is.na(Count_source))
  ) 

# Transpose the summary statistics
count_source_stats_t <- count_source_stats %>%
  pivot_longer(cols = everything(), names_to = "Statistic", values_to = "Value")

# Display results in RStudio
View(count_source_stats_t)  # Open in Viewer
print(count_source_stats_t) # Print in console

#####
#####
#####
#####
##### Figure 4.Accumulated histogram illustrating the distribution of multi-omic signatures 
##### by multi-omic feature across various cancer types
##### Visualizing Cancer Type and Genomic Feature Distributions with Okabe-Ito Palette for 
##### Enhanced Accessibility in Data Interpretation

# Load required packages
library(ggplot2)
library(dplyr)

# Define Okabe-Ito palette (color-blind friendly)
okabe_ito_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                       "#0072B2", "#D55E00", "#CC79A7")

# Data Preprocessing with Renamed Columns
df_plot <- data_df1011 %>%
  group_by(`Cancer Type Abbreviation` = CTAB, `Genomic Feature` = Genotype) %>%
  summarise(Absolute_Count = n()) %>%
  ungroup()

# Plotting Accumulated Histogram with Absolute Counts and color-blind friendly palette
signature_distribution_plot <- ggplot(df_plot, aes(x = `Cancer Type Abbreviation`, y = Absolute_Count, fill = `Genomic Feature`)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = okabe_ito_palette) +  # Using Okabe-Ito color-blind palette
  theme_classic() + # Ensures a white background
  theme(
    panel.background = element_rect(fill = "white", color = "white"), # Set background to white
    plot.background = element_rect(fill = "white", color = "white"),   # Set plot background to white
    axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5), 
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 14)
  ) +
  labs(
    x = "Cancer Type Abbreviation",
    y = "Absolute Count (Accumulated)",
    fill = "Genomic Feature"
  )
signature_distribution_plot

# Corrected ggsave function usage
ggsave("Signature_Distribution.tif", plot = signature_distribution_plot, 
       width = 29, height = 21, units = "cm", dpi = 600)


#####
#####
#####
#####
##### Figure 4B.Accumulated histogram illustrating the distribution of multi-omic signatures 
##### by multi-omic feature across various cancer types
##### Visualizing Cancer Type and Genomic Feature Distributions with Okabe-Ito Palette for 
##### Enhanced Accessibility in Data Interpretation

# Load required packages
library(ggplot2)
library(dplyr)

# Define Okabe-Ito palette (color-blind friendly)
okabe_ito_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                       "#0072B2", "#D55E00", "#CC79A7")

# Data Preprocessing with Renamed Columns
df_plot <- data_df1011 %>%
  group_by(`Cancer Type Abbreviation` = CTAB, `Genomic Feature` = Genotype) %>%
  summarise(Absolute_Count = n()) %>%
  ungroup()

# Plotting Accumulated Histogram with Absolute Counts and color-blind friendly palette
signature_distribution_plot <- ggplot(df_plot, aes(x = `Cancer Type Abbreviation`, y = Absolute_Count, fill = `Genomic Feature`)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = okabe_ito_palette) +  # Using Okabe-Ito color-blind palette
  theme_classic() + # Ensures a white background
  theme(
    panel.background = element_rect(fill = "white", color = "white"), # Set background to white
    plot.background = element_rect(fill = "white", color = "white"),   # Set plot background to white
    axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5), 
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 14)
  ) +
  labs(
    x = "Cancer Type Abbreviation",
    y = "Absolute Count (Accumulated)",
    fill = "Genomic Feature"
  )
signature_distribution_plot

# Corrected ggsave function usage
ggsave("Figire_4_B_Signature_Distribution.tif", plot = signature_distribution_plot, 
       width = 29, height = 21, units = "cm", dpi = 600)

##### Figure 4C.Accumulated histogram illustrating the distribution of multi-omic signatures 
##### by multi-omic feature across various cancer types
##### Visualizing Cancer Type and Genomic Feature Distributions with Okabe-Ito Palette for 
##### Enhanced Accessibility in Data Interpretation

# Load required packages
library(ggplot2)
library(dplyr)
library(viridis)  # Load the viridis package for color-blind friendly palettes

# Data Preprocessing with Renamed Columns
df_plot <- data_df1011 %>%
  group_by(`Cancer Type Abbreviation` = CTAB, `Genomic Feature` = Genotype) %>%
  summarise(Absolute_Count = n()) %>%
  ungroup()

# Plotting Accumulated Histogram with Absolute Counts and viridis color-blind friendly palette
signature_distribution_plot <- ggplot(df_plot, aes(x = `Cancer Type Abbreviation`, y = Absolute_Count, fill = `Genomic Feature`)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_viridis(discrete = TRUE) +  # Using viridis color-blind friendly palette
  theme_classic() + # Ensures a white background
  theme(
    panel.background = element_rect(fill = "white", color = "white"), # Set background to white
    plot.background = element_rect(fill = "white", color = "white"),   # Set plot background to white
    axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5), 
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 14)
  ) +
  labs(
    x = "Cancer Type Abbreviation",
    y = "Absolute Count (Accumulated)",
    fill = "Genomic Feature"
  )
signature_distribution_plot

# Save the plot
ggsave("Signature_Distribution.tif", plot = signature_distribution_plot, 
       width = 29, height = 21, units = "cm", dpi = 600)

##### Figure 4D.Accumulated histogram illustrating the distribution of multi-omic signatures 
##### by multi-omic feature across various cancer types
##### Visualizing Cancer Type and Genomic Feature Distributions with Okabe-Ito Palette for 
##### Enhanced Accessibility in Data Interpretation
# Load required packages
library(ggplot2)
library(dplyr)
library(RColorBrewer)  # Load the RColorBrewer package for color-blind friendly palettes

# Data Preprocessing with Renamed Columns
df_plot <- data_df1011 %>%
  group_by(`Cancer Type Abbreviation` = CTAB, `Genomic Feature` = Genotype) %>%
  summarise(Absolute_Count = n()) %>%
  ungroup()

# Define the ColorBrewer "Paired" palette
paired_palette <- brewer.pal(12, "Paired")  # 12 colors in the "Paired" palette

# Plotting Accumulated Histogram with Absolute Counts and ColorBrewer "Paired" palette
signature_distribution_plot <- ggplot(df_plot, aes(x = `Cancer Type Abbreviation`, y = Absolute_Count, fill = `Genomic Feature`)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = paired_palette) +  # Using ColorBrewer "Paired" palette
  theme_classic() + # Ensures a white background
  theme(
    panel.background = element_rect(fill = "white", color = "white"), # Set background to white
    plot.background = element_rect(fill = "white", color = "white"),   # Set plot background to white
    axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5), 
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 14)
  ) +
  labs(
    x = "Cancer Type Abbreviation",
    y = "Absolute Count (Accumulated)",
    fill = "Genomic Feature"
  )
signature_distribution_plot

# Save the plot
ggsave("Signature_Distribution.tif", plot = signature_distribution_plot, 
       width = 29, height = 21, units = "cm", dpi = 600)

##### Figure 4E.Accumulated histogram illustrating the distribution of multi-omic signatures 
##### by multi-omic feature across various cancer types
##### Visualizing Cancer Type and Genomic Feature Distributions with Okabe-Ito Palette for 
##### Enhanced Accessibility in Data Interpretation
# Load required packages
library(ggplot2)
library(dplyr)

# Data Preprocessing with Renamed Columns
df_plot <- data_df1011 %>%
  group_by(`Cancer Type Abbreviation` = CTAB, `Genomic Feature` = Genotype) %>%
  summarise(Absolute_Count = n()) %>%
  ungroup()

# Plotting Accumulated Histogram with Absolute Counts and ggplot2 default palette
signature_distribution_plot <- ggplot(df_plot, aes(x = `Cancer Type Abbreviation`, y = Absolute_Count, fill = `Genomic Feature`)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_hue() +  # Use ggplot2's default color palette
  theme_classic() + # Ensures a white background
  theme(
    panel.background = element_rect(fill = "white", color = "white"), # Set background to white
    plot.background = element_rect(fill = "white", color = "white"),   # Set plot background to white
    axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5), 
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 14)
  ) +
  labs(
    x = "Cancer Type Abbreviation",
    y = "Absolute Count (Accumulated)",
    fill = "Genomic Feature"
  )
signature_distribution_plot

# Save the plot
ggsave("Figure_4_Signature_Distribution.tif", plot = signature_distribution_plot, 
       width = 29, height = 21, units = "cm", dpi = 600)

#####
#####


#####
#####
#####
#####
####Top-Ranked Cancer Types by Genomic Feature in Multi-Omic Signatures
# Load necessary package
library(dplyr)

# Aggregate data: Count the number of signatures per Cancer Type (CTAB) and Genotype
top_cancer_types <- data_df1011 %>%
  group_by(CTAB, Genotype) %>%
  summarise(Signature_Count = n(), .groups = "drop") %>%
  arrange(Genotype, desc(Signature_Count))

# Identify the top-ranked cancer types per Genotype
top_ranked <- top_cancer_types %>%
  group_by(Genotype) %>%
  slice_max(order_by = Signature_Count, n = 1) %>%
  ungroup()

# Display results in RStudio
View(top_ranked)  # Opens in RStudio Viewer
# OR print in console
print(top_ranked)



# Define the Okabe-Ito color-blind friendly palette
okabe_ito_palette <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#999999"
)

# Prepare data_df1011 for plotting
plot_data_df1011 <- data_df1011 %>%
  count(Combined_Outcome, Genotype) %>%
  arrange(Combined_Outcome, Genotype)

# Create a stacked bar plot
ggplot(plot_data_df1011, aes(x = Combined_Outcome, y = n, fill = Genotype)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = okabe_ito_palette) +
  theme_minimal() +
  labs(
    title = "Stacked Distribution by Combined Outcome and Genotype",
    x = "Combined Outcome",
    y = "Number of Rows",
    fill = "Genotype"
  ) +
  theme(
    text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


# Exclude rows with "Meaningless" in Combined_Outcome
filtered_data_df1011 <- data_df1011 %>%
  filter(Combined_Outcome != "Meaningless")

# Define Okabe-Ito color palette
okabe_ito_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                       "#0072B2", "#D55E00", "#CC79A7")

# Create the plot
stacked_plot <- ggplot(filtered_data_df1011, aes(x = Genotype, fill = Combined_Outcome)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = okabe_ito_palette) +
  labs(
    title = "Stacked Distribution of Combined Outcome by Genotype",
    x = "Genotype",
    y = "Proportion",
    fill = "Combined Outcome"
  ) +
  theme_minimal(base_size = 14)

# Save the plot as TIFF
ggsave(
  filename = "stacked_distribution_plot.tiff",
  plot = stacked_plot,
  device = "tiff",
  dpi = 600,
  width = 8.27, height = 11.69 # A4 dimensions in inches
)

# Save the plot as PDF
ggsave(
  filename = "stacked_distribution_plot.pdf",
  plot = stacked_plot,
  device = "pdf",
  width = 8.27, height = 11.69 # A4 dimensions in inches
)

cat("Plot saved as 'stacked_distribution_plot.tiff' and 'stacked_distribution_plot.pdf' successfully!\n")

# Load necessary libraries
library(ggplot2)

# Exclude "Meaningless" from "Combined_Outcome"
filtered_data_df1011 <- data_df1011 %>%
  filter(Combined_Outcome != "Meaningless")

# Create the proportional distribution plot
plot <- ggplot(filtered_data_df1011, aes(x = Combined_Outcome, fill = Genotype)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    x = "Combined Outcome",
    y = "Proportional Distribution",
    title = "Proportional Distribution of Genotype by Combined Outcome",
    fill = "Genotype"
  ) +
  scale_fill_manual(values = c(
    "#E69F00", "#56B4E9", "#009E73", "#F0E442",
    "#0072B2", "#D55E00", "#CC79A7"  # Okabe-Ito palette
  )) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
  )

# Save the plot as PNG
ggsave("Proportional_Distribution_by_Combined_Outcome.png", plot = plot, 
       dpi = 600, width = 8.3, height = 11.7, units = "in")  # A4 dimensions

# Save the plot as TIFF
ggsave("Proportional_Distribution_by_Combined_Outcome.tiff", plot = plot, 
       dpi = 600, width = 8.3, height = 11.7, units = "in")  # A4 dimensions

# Print the plot to the console
print(plot)

dev.off()
#####
#####
##### snippet to create an interactive proportional distribution plot of "Genotype" by "Combined_Outcome," excluding the "Meaningless" category.
#####
##### 

# Exclude "Meaningless" from "Combined_Outcome"
filtered_data_df1011 <- data_df1011 %>%
  filter(Combined_Outcome != "Meaningless")

# Create the proportional distribution plot
ggplot_obj <- ggplot(filtered_data_df1011, aes(x = Combined_Outcome, fill = Genotype)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    x = "Combined Outcome",
    y = "Proportional Distribution",
    title = "Proportional Distribution of Genotype by Combined Outcome",
    fill = "Genotype"
  ) +
  scale_fill_manual(values = c(
    "#E69F00", "#56B4E9", "#009E73", "#F0E442",
    "#0072B2", "#D55E00", "#CC79A7"  # Okabe-Ito palette
  )) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
  )

# Convert ggplot to an interactive plotly object
interactive_plot <- ggplotly(ggplot_obj)

# Display the interactive plot
interactive_plot

######
######
######
###### Autophagy RPPA Victor genes
######
######
######
setwd("D:/Pré-artigo 5-optosis model/Dataframes all metrics/Victor_genes")

df500 <- import("lista_final_genes_validados.txt")

setwd("D:/Pré-artigo 5-optosis model/Dataframes all metrics")

df501 <- unique_RPPA_array_mapped

# Assuming df500 and df501 are your dataframes
filtered_df500 <- df500[df500$Genes %in% df501$Gene_Symbol, ]

# Ensure it is a dataframe and explicitly name the variable as "Genes"
filtered_df500 <- as.data.frame(filtered_df500)
colnames(filtered_df500)[colnames(filtered_df500) == "filtered_df500"] <- "Genes"

df502 <- filtered_df500

# Assuming df500 and df501 are your dataframes
filtered_df501<- df501[df501$Gene_Symbol %in% df502$Genes, ]

protein_info_autophagy <- filtered_df501


rio::export(protein_info_autophagy, "protein_info_autophagy.xlsx")

# Save protein_info_autophagy as an Excel file
write.xlsx(protein_info_autophagy, file = "protein_info_autophagy.xlsx")


######
###### Table 7. Selected top of the top for final figures
######
# Define the specific values to filter
nomenclature_values <- c("KIRC-169.2.1.P.2.71.45.1.1.2",
                         "BRCA-1496.1.3.P.3.71.71.1.1.2",
                         "HNSC-1855.4.3.P.3.71.64.1.1.4",
                         "HNSC-156.7.3.N.3.71.55.1.1.1",
                         "KIRP-107.3.2.N.1.44.44.1.1.2",
                         "CESC-283.6.3.N.2.44.44.1.1.3",
                         "CESC-215.5.3.N.2.44.44.1.1.3")

# Filter the rows and arrange by Rank in descending order
filtered_data_top <- data_df1011 %>%
  filter(Nomenclature %in% nomenclature_values) %>%
  arrange(desc(Ranking))

df11053 <- filtered_data_top

rio::export(df11053, "Table 7 Selected top of the top for final figures.xlsx")

# Save top top  as an Excel file
# write.xlsx(df11053, file = "Table 7 Selected top of the top for final figures.xlsx")

####
####
####
#### Revised! Genes in top signatures tables 1, 2, 3, 6 without duplicates: New compilation analysis
####
####
####
# PART A
df1152 <-final_renamed_ordered_df1011
df1153 <- selected_columns_df1047_reordered
df1154 <- selected_columns_combined_genotype_rows
df1155 <- reordered_filtered_columns_top_MOST_ranked_signatures_df11052

# Rename variables in df155
df1155 <- df1155 %>%
  rename(
    `Multi-omic feature` = `Genomic feature`,
    `RCD forms` = `RCD form`
  )

library(dplyr)

# Select the specified columns in each dataframe and ensure "Rank" is numeric
df1152_selected <- df1152 %>%
  select(Rank, Nomenclature, Signature, `Multi-omic feature`, `RCD forms`) %>%
  mutate(Rank = as.numeric(Rank))

df1153_selected <- df1153 %>%
  select(Rank, Nomenclature, Signature, `Multi-omic feature`, `RCD forms`) %>%
  mutate(Rank = as.numeric(Rank))

df1154_selected <- df1154 %>%
  select(Rank, Nomenclature, Signature, `Multi-omic feature`, `RCD forms`) %>%
  mutate(Rank = as.numeric(Rank))

df1155_selected <- df1155 %>%
  select(Rank, Nomenclature, Signature, `Multi-omic feature`, `RCD forms`) %>%
  mutate(Rank = as.numeric(Rank))

# Check if variable names and order are identical in all four dataframes
identical_names <- identical(names(df1152_selected), names(df1153_selected)) &&
  identical(names(df1153_selected), names(df1154_selected)) &&
  identical(names(df1154_selected), names(df1155_selected))

# If names are identical, bind the dataframes
if (identical_names) {
  df1156 <- bind_rows(df1152_selected, df1153_selected, df1154_selected, df1155_selected)
} else {
  stop("The column names or order are not identical across the dataframes.")
}

# Check for duplicate rows across all columns
duplicates_exist <- any(duplicated(df1156))

if (duplicates_exist) {
  cat("Duplicate rows found. Removing duplicates...\n")
  
  # Remove duplicate rows while retaining one instance
  df1156 <- df1156 %>%
    distinct()
  
  cat("Duplicates removed. Updated dataframe:\n")
  print(df1156)
} else {
  cat("No duplicate rows found.\n")
}

# Check for duplicate rows across all columns
duplicate_indices <- which(duplicated(df1156) | duplicated(df1156, fromLast = TRUE))  # Identify all duplicate rows

if (length(duplicate_indices) > 0) {
  cat("Duplicate rows found. Printing duplicate rows:\n")
  
  # Print the duplicated rows
  duplicated_rows <- df1156[duplicate_indices, ]
  print(duplicated_rows)
  
  # Remove duplicate rows while retaining one instance
  df1156 <- df1156 %>%
    distinct()
  
  cat("Duplicates removed. Updated dataframe:\n")
  print(df1156)
} else {
  cat("No duplicate rows found.\n")
}

# Rename, reorder, and sort the columns in descending order by Rank
df1156 <- df1156 %>%
  # Reorder columns
  arrange(desc(Rank))  # Sort in descending order by Rank

#######################
#######################
##### Re-Analysis of top signatures ####
#######################
#######################
# setwd("D:/Pré-artigo 5-optosis model/Dataframes all metrics")

data <- import("Drafts_2/gene_info_trancritos_excel2.xlsx")

data <- data %>%
  mutate_all(trimws)

miRNA_database <-  import("Drafts_2/miRNA_database.xlsx")

miRNA_database <- miRNA_database %>%
  mutate_all(trimws)

protein_info <- import("Drafts_2/protein_info.xlsx") 

protein_info <- protein_info %>%
  mutate_all(trimws)

df1156 <-df1156 %>%
  mutate_all(trimws)

data_df1011 <- data_df1011 %>%
  mutate_all(trimws)

# Step 0: Filter data_df1011 based on matching values in the "Nomenclature" column of data_df1011 and "Nomenclature" column in df1156
data_df1011_filtered <- data_df1011 %>%
  filter(Nomenclature %in%df1156$Nomenclature)

# Step 1: Copy the filtered columns "Signature" and "Nomenclature" from data_df1011_filtered to df38
df38 <- data_df1011_filtered %>%
  select(Signature, Nomenclature)

# Step 2: Create a new column "Gene_conversion" after "Signature"
df38 <- df38 %>%
  mutate(Gene_conversion = NA)  # Initialize the column with NA values

# Step 3: Populate "Gene_conversion" using left_join for better matching reliability

# For values starting with "ENST", join with the "data" dataframe on "Transcript_ID"
df38 <- df38 %>%
  left_join(data %>% select(Transcript_ID, Display_Name), by = c("Signature" = "Transcript_ID")) %>%
  mutate(Gene_conversion = ifelse(!is.na(Display_Name), Display_Name, Gene_conversion)) %>%
  select(-Display_Name)  # Remove the temporary column

# For values starting with "hsa-miR", join with the "miRNA_database" on "Mature1" and "Mature2"
df38 <- df38 %>%
  left_join(miRNA_database %>% select(Mature1, Gene_Symbol), by = c("Signature" = "Mature1")) %>%
  mutate(Gene_conversion = ifelse(!is.na(Gene_Symbol), Gene_Symbol, Gene_conversion)) %>%
  select(-Gene_Symbol)  # Remove the temporary column

df38 <- df38 %>%
  left_join(miRNA_database %>% select(Mature2, Gene_Symbol), by = c("Signature" = "Mature2")) %>%
  mutate(Gene_conversion = ifelse(!is.na(Gene_Symbol), Gene_Symbol, Gene_conversion)) %>%
  select(-Gene_Symbol)  # Remove the temporary column

# For protein values join with the "protein_info" on "Protein1" 
df38 <- df38 %>%
  left_join(protein_info %>% select(Protein1, Gene_Symbol), by = c("Signature" = "Protein1")) %>%
  mutate(Gene_conversion = ifelse(!is.na(Gene_Symbol), Gene_Symbol, Gene_conversion)) %>%
  select(-Gene_Symbol)  # Remove the temporary column

# Step 1: Extract gene symbols, transcript names, and miRNA from the "Signature" column
df39 <- df38 %>%
  # Use str_replace_all to remove parentheses and separate by "+" into individual elements
  mutate(Signature_cleaned = str_replace_all(Signature, "[()]", "")) %>%  # Remove parentheses
  separate_rows(Signature_cleaned, sep = " \\+ ") %>%  # Split by " + " into separate rows
  distinct() %>%  # Remove duplicate rows
  filter(str_detect(Signature_cleaned, "^ENST|^hsa-miR|^[A-Za-z]")) %>%  # Keep only gene symbols, ENST transcripts, or miRNA
  arrange(Signature_cleaned)  # Order alphabetically

# Step 2: Create the df39 with the unique and ordered entries
df39 <- df39 %>%
  select(Signature_cleaned) %>%
  rename(Extracted_Signature = Signature_cleaned)

# Step 1: Create a new column "Gene_Symbol" in df39 and populate it based on the mapping logic
df39 <- df39 %>%
  rowwise() %>%
  mutate(Gene_Symbol = case_when(
    
    # For values with prefix "ENST", look up in the "data" dataframe using "Transcript_ID"
    str_detect(Extracted_Signature, "^ENST") ~ 
      ifelse(Extracted_Signature %in% data$Transcript_ID, 
             data$Display_Name[which(data$Transcript_ID == Extracted_Signature)], 
             NA_character_),
    
    # For values with prefix "hsa-miR", look up in "miRNA_database" for Mature1
    str_detect(Extracted_Signature, "^hsa-miR") & Extracted_Signature %in% miRNA_database$Mature1 ~ 
      ifelse(Extracted_Signature %in% miRNA_database$Mature1, 
             miRNA_database$Gene_Symbol[which(miRNA_database$Mature1 == Extracted_Signature)], 
             NA_character_),
    
    # For values with prefix "hsa-miR", look up in "miRNA_database" for Mature2
    str_detect(Extracted_Signature, "^hsa-miR") & Extracted_Signature %in% miRNA_database$Mature2 ~ 
      ifelse(Extracted_Signature %in% miRNA_database$Mature2, 
             miRNA_database$Gene_Symbol[which(miRNA_database$Mature2 == Extracted_Signature)], 
             NA_character_),
    
    # If no match is found, leave as NA
    TRUE ~ NA_character_
  ))

# Fill NA values in "Gene_Symbol" with corresponding values from "Extracted_Signature"
df39 <- df39 %>%
  mutate(Gene_Symbol = ifelse(is.na(Gene_Symbol), Extracted_Signature, Gene_Symbol))

# Remove duplicates based on "Gene_Symbol"
df39 <- df39 %>%
  distinct(Gene_Symbol, .keep_all = TRUE)

# Order the dataframe by "Gene_Symbol" alphabetically
df39 <- df39 %>%
  arrange(Gene_Symbol)

# Step 1: Identify the common values between "Protein1" and "Gene_symbol"
common_proteins <- intersect(protein_info$Protein1, df39$Gene_Symbol)

# Step 2: Create a look up table for replacements
protein_lookup <- protein_info %>%
  filter(Protein1 %in% common_proteins) %>%
  select(Protein1) %>%
  distinct()

# Step 3: Replace the "Gene_symbol" values in df39 with the corresponding "Protein1" values from protein_info
df39 <- df39 %>%
  mutate(Gene_Symbol = ifelse(Gene_Symbol %in% protein_info$Protein1,  # Check if Gene_symbol is in Protein1
                              protein_info$Gene_Symbol[match(Gene_Symbol, protein_info$Protein1)],  # Replace with corresponding Gene_symbol in protein_info
                              Gene_Symbol))  # Keep original Gene_symbol if no match

# Exclude the "Extracted_Signature" column
df39 <- df39 %>%
  select(-Extracted_Signature)

#### Handling "PRKN", "PARK2" ambiguity (02/24/2025) 
#### This is so because: PARK2 is the gene Symbol that codes the Parkin protein, also known as PRKN 
# Trim whitespace, replace "PRKN" with "PARK2", and remove duplicates
df39 <- df39 %>%
  mutate(Gene_Symbol = trimws(Gene_Symbol)) %>%  # Trim whitespace
  mutate(Gene_Symbol = ifelse(Gene_Symbol == "PRKN", "PARK2", Gene_Symbol)) %>%  # Replace PRKN
  distinct()  # Remove duplicate rows

# Print confirmation message
print("Replaced 'PRKN' with 'PARK2', removed duplicates, and trimmed whitespace in 'Gene_Symbol' column.")

# Ensure the dataframes have been loaded or created: Target_genes and df39

# Filter Target_genes based on matching values in Gene_Symbol from df39
filtered_Target_genes <- Target_genes[Target_genes$`Gene symbol` %in% df39$Gene_Symbol, ]

# Identify rows in df39 where Gene_Symbol is NOT in Target_genes$Gene_Symbol
non_matching_rows <- df39[!df39$Gene_Symbol %in% Target_genes$`Gene symbol`, ]

# View the non-matching rows
head(non_matching_rows)

rio::export(filtered_Target_genes, "Dataset S1Q_NOVO_filtered_Target_genes_Gene_Symbols_top_signatures_df1156.xlsx")

# Assign dataframe to Gene_Symbols_top_signatures_df1156
Gene_Symbols_top_signatures_df1156 <- df39

# Order the dataframe alphabetically by the "Gene_Symbol" column
Gene_Symbols_top_signatures_df1156 <- Gene_Symbols_top_signatures_df1156[order(Gene_Symbols_top_signatures_df1156$Gene_Symbol), ]

# Remove leading and trailing spaces globally across all columns
Gene_Symbols_top_signatures_df1156 <- Gene_Symbols_top_signatures_df1156 %>%
  mutate_all(trimws)

#  Remove duplicate rows
duplicates_exist <- any(duplicated(Gene_Symbols_top_signatures_df1156))

if (duplicates_exist) {
  cat("Duplicate rows found. Removing duplicates...\n")
  
  # Remove duplicate rows while retaining one instance
  Gene_Symbols_top_signatures_df1156 <- Gene_Symbols_top_signatures_df1156 %>%
    distinct()
  
  cat("Duplicates removed. Updated dataframe:\n")
  print(Gene_Symbols_top_signatures_df1156)
} else {
  cat("No duplicate rows found in the dataframe.\n")
}

# Ensure the dataframes have been loaded or created: Target_genes and df39

# Filter Target_genes based on matching values in Gene_Symbol from df39
filtered_Target_genes <- Target_genes[Target_genes$"Gene symbol" %in% df39$Gene_Symbol, ]

# Identify rows in df39 where Gene_Symbol is NOT in Target_genes$Gene_Symbol
non_matching_rows <- df39[!df39$Gene_Symbol %in% Target_genes$"Gene symbol", ]

# View the non-matching rows
head(non_matching_rows)

# Load required libraries
library(readr)    # For writing TSV files
library(openxlsx) # For writing Excel files

# Step 5: rio::export the cleaned dataframe to a TSV file
rio::export(Gene_Symbols_top_signatures_df1156, "Gene_Symbols_top_signatures_df1156.tsv")
write_tsv(Gene_Symbols_top_signatures_df1156, "Gene_Symbols_top_signatures_df1156.tsv")

# Step 6: rio::export the cleaned dataframe to an Excel file
write.xlsx(Gene_Symbols_top_signatures_df1156, "Gene_Symbols_top_signatures_df1156.xlsx", overwrite = TRUE)

# Step 5: rio::export the cleaned dataframe to a TSV file
rio::export(Gene_Symbols_top_signatures_df1156, "Gene_Symbols_top_signatures_df1156.tsv")
rio::export(Gene_Symbols_top_signatures_df1156, "Gene_Symbols_top_signatures_df1156.xlsx")

cat("Cleaned dataframe rio::exported to 'Gene_Symbols_top_signatures_df1156.tsv'.\n")

##### 
##### 
##### 
##### DGIdb 5.0: drug-gene interaction 
##### Making Dataset XXX. filtered_GDIdb_meaningful
#####
#####
setwd("D:/Pré-artigo 5-optosis model/DGIdb_05022025")

GDIdb <- import("interactions.tsv")

setwd("D:/Pré-artigo 5-optosis model/Dataframes all metrics")

GDIdb_genes <- as.data.frame(unique(GDIdb$gene_name))

# Load necessary package
library(dplyr)

# Trim leading and trailing spaces in the "gene_name" column of GDIdb
GDIdb <- GDIdb %>%
  mutate(gene_name = trimws(gene_name))

# Trim leading and trailing spaces in the "Gene_Symbol" column of df39
df39 <- df39 %>%
  mutate(Gene_Symbol = trimws(Gene_Symbol))

# Filter GDIdb for rows where "gene_name" matches "Gene_Symbol"
filtered_GDIdb <- GDIdb %>%
  filter(gene_name %in% df39$Gene_Symbol)

filtered_GDIdb_meaningful <- filtered_GDIdb %>%
  filter(interaction_type != "NULL" & interaction_type != "other/unknown")

GDIdb_gene_distribution <- filtered_GDIdb_meaningful %>%
  group_by(interaction_type, gene_name) %>%
  summarise(count = n(), .groups = 'drop')

gene_counts <- filtered_GDIdb_meaningful %>%
  group_by(gene_name) %>%
  summarise(total_count = n(), .groups = 'drop')

# Filter genes with at least 5 counts
genes_with_min_counts <- gene_counts[gene_counts$total_count >= 5, ]

# Order by total_count in descending order
genes_with_min_counts <- genes_with_min_counts[order(-genes_with_min_counts$total_count), ]

# Extract the ordered gene names
ordered_gene_list <- genes_with_min_counts$gene_name

# View the result
print(ordered_gene_list)

set_genes_with_min_counts <- filtered_GDIdb_meaningful %>%
  group_by(gene_name) %>%
  filter(n() >= 5) %>%
  ungroup()

# Load necessary library
library(ggplot2)

# Create the stacked histogram
ggplot(filtered_GDIdb_meaningful, aes(x = gene_name, fill = interaction_type)) +
  geom_histogram(stat = "count", position = "stack", color = "black") +
  labs(
    title = "Distribution of Gene Names by Interaction Type",
    x = "Gene Name",
    y = "Count",
    fill = "Interaction Type"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    plot.title = element_text(hjust = 0.5)
  )

# Load necessary libraries
library(ggplot2)
library(forcats)

# Define a color palette with sufficient unique colors
interaction_types <- unique(filtered_GDIdb_meaningful$interaction_type)
num_interaction_types <- length(interaction_types)

# Generate a color palette with the required number of unique colors
# You can use the 'hue_pal' function from the 'scales' package to generate distinct colors
library(scales)
color_palette <- hue_pal()(num_interaction_types)

# Create a named vector for the color palette
names(color_palette) <- interaction_types

# Create the stacked bar plot with the specified color palette
ggplot(filtered_GDIdb_meaningful, aes(x = fct_infreq(gene_name), fill = interaction_type)) +
  geom_bar(position = "stack", color = "black") +
  labs(
    title = "Distribution of Gene Names by Interaction Type",
    x = "Gene Name",
    y = "Count",
    fill = "Interaction Type"
  ) +
  scale_fill_manual(values = color_palette) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    plot.title = element_text(hjust = 0.5)
  )

# Install ggpattern if not already installed
if (!requireNamespace("ggpattern", quietly = TRUE)) {
  install.packages("ggpattern", repos = "https://cinc.rud.is")
}

#### Creating bar histogram using different colors and different patterns
# Load required libraries
library(ggplot2)
library(ggpattern)
library(forcats)

# Extract unique interaction types
interaction_types <- unique(filtered_GDIdb_meaningful$interaction_type)
num_interaction_types <- length(interaction_types)

# Define a color palette with sufficient unique colors
color_palette <- scales::hue_pal()(num_interaction_types)

# Define a set of patterns
available_patterns <- c("stripe", "crosshatch", "wave", "weave", "circle", "none")
pattern_palette <- rep(available_patterns, length.out = num_interaction_types)

# Create named vectors for colors and patterns
names(color_palette) <- interaction_types
names(pattern_palette) <- interaction_types

# Create the plot
plot_GDIdb_genes <- ggplot(filtered_GDIdb_meaningful, aes(x = fct_infreq(gene_name), fill = interaction_type, pattern = interaction_type)) +
  geom_bar_pattern(
    position = "stack",
    color = "black",
    pattern_density = 0.3,  # Adjusted for better clarity
    pattern_fill = "grey80",  # Adjusted for better contrast
    pattern_spacing = 0.02
  ) +
  scale_fill_manual(values = color_palette) +
  scale_pattern_manual(values = pattern_palette) +
  labs(
    title = "Distribution of Gene Names by Interaction Type",
    x = "Gene Name",
    y = "Count",
    fill = "Interaction Type",
    pattern = "Interaction Type"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  )

# Save the plot as PDF
ggsave(
  filename = "plot_GDIdb_genes.pdf",
  plot = plot_GDIdb_genes,
  device = "pdf",
  width = 11.69,  # A4 width in inches (landscape orientation)
  height = 8.27,  # A4 height in inches
  dpi = 600
)

# Save the plot as TIFF
ggsave(
  filename = "plot_GDIdb_genes.tiff",
  plot = plot_GDIdb_genes,
  device = "tiff",
  width = 11.69,
  height = 8.27,
  dpi = 600
)

rio::export(filtered_GDIdb_meaningful, "filtered_GDIdb_meaningful.xlsx")

dev.off()

#### 
#### 
#### 
#### DGI Network analysis plotting, Part A. using the meaningful set of interactions
#### 
#### 
#### 
#### # Load required packages
library(igraph)
library(ggraph)
library(tidyverse)

# Filter meaningful interactions (excluding unknown interactions)
filtered_GDIdb_network_meaningful <- filtered_GDIdb %>%
  filter(interaction_type != "NULL" & interaction_type != "other/unknown") %>%
  select(gene_name, drug_name, interaction_type) %>%
  distinct()

# Create an edge list (interactions between genes and drugs)
edges <- filtered_GDIdb_network_meaningful %>%
  select(gene_name, drug_name) %>%
  distinct()

# Create a graph object
g <- graph_from_data_frame(edges, directed = FALSE)

# Assign node types (Gene or Drug)
V(g)$type <- ifelse(V(g)$name %in% filtered_GDIdb_network_meaningful$gene_name, "Gene", "Drug")

# Assign colors based on node type
node_colors <- ifelse(V(g)$type == "Gene", "steelblue", "tomato")

# Plot using ggraph
ggraph(g, layout = "fr") + 
  geom_edge_link(aes(edge_alpha = 0.5), color = "gray50") + 
  geom_node_point(aes(color = V(g)$type), size = 5) + 
  geom_node_text(aes(label = name), repel = TRUE, size = 4) + 
  scale_color_manual(values = c("Gene" = "steelblue", "Drug" = "tomato")) + 
  theme_void() + 
  ggtitle("Drug-Gene Interaction Network")


######
######
######
###### Part B network coloring by drug effect type
###### 
###### 
# Load required packages
library(igraph)
library(ggraph)
library(tidyverse)

# Filter meaningful interactions (excluding unknown interactions)
filtered_GDIdb_network_meaningful <- filtered_GDIdb %>%
  filter(interaction_type != "NULL" & interaction_type != "other/unknown") %>%
  select(gene_name, drug_name, interaction_type) %>%
  distinct()

# Create an edge list (interactions between genes and drugs)
edges <- filtered_GDIdb_network_meaningful %>%
  select(gene_name, drug_name) %>%
  distinct()

# Create a graph object
g <- graph_from_data_frame(edges, directed = FALSE)

# Assign node types (Gene or Drug)
V(g)$type <- ifelse(V(g)$name %in% filtered_GDIdb_network_meaningful$gene_name, "Gene", "Drug")

# Assign interaction types to drugs
drug_interaction_mapping <- filtered_GDIdb_network_meaningful %>%
  distinct(drug_name, interaction_type) %>%
  group_by(drug_name) %>%
  summarise(interaction_type = first(interaction_type), .groups = "drop") # Keep first type if multiple

# Merge drug interaction types into graph nodes
V(g)$interaction_type <- ifelse(V(g)$name %in% drug_interaction_mapping$drug_name,
                                drug_interaction_mapping$interaction_type[match(V(g)$name, drug_interaction_mapping$drug_name)],
                                NA)

# Generate a color palette for different interaction types
interaction_types <- unique(na.omit(V(g)$interaction_type))
interaction_colors <- setNames(scales::hue_pal()(length(interaction_types)), interaction_types)

# Assign colors: Genes in "steelblue", Drugs colored by interaction type
node_colors <- ifelse(V(g)$type == "Gene", "steelblue", interaction_colors[V(g)$interaction_type])

# Plot using ggraph
Network_GDIdb_genes <- ggraph(g, layout = "fr") + 
  geom_edge_link(aes(edge_alpha = 0.5), color = "gray50") + 
  geom_node_point(aes(color = factor(V(g)$interaction_type, levels = interaction_types)), size = 5) + 
  geom_node_text(aes(label = name), repel = TRUE, size = 4) + 
  scale_color_manual(name = "Drug Effects", values = c("Gene" = "steelblue", interaction_colors)) +  # Updated legend title
  theme_void() + 
  ggtitle("Drug-Gene Interaction Network: Colored by Interaction Type") + 
  theme(legend.position = "right")

# Save the Network plot as PDF
ggsave(
  filename = "Network_plot_GDIdb_genes.pdf",
  plot = Network_GDIdb_genes,
  device = "pdf",
  width = 11.69,  # A4 width in inches (landscape orientation)
  height = 8.27,  # A4 height in inches
  dpi = 600
)

# Save the Network plot as TIFF
ggsave(
  filename = "Network_plot_GDIdb_genes.tiff",
  plot = Network_GDIdb_genes,
  device = "tiff",
  width = 11.69,  # A4 width in inches (landscape orientation)
  height = 8.27,  # A4 height in inches
  dpi = 600
)

######
######
###### Part C improving display
###### 
# Load required packages
library(igraph)
library(ggraph)
library(tidyverse)

# Filter meaningful interactions (excluding unknown interactions)
filtered_GDIdb_network_meaningful <- filtered_GDIdb %>%
  filter(interaction_type != "NULL" & interaction_type != "other/unknown") %>%
  select(gene_name, drug_name, interaction_type) %>%
  distinct()

# Create an edge list (interactions between genes and drugs)
edges <- filtered_GDIdb_network_meaningful %>%
  select(gene_name, drug_name) %>%
  distinct()

# Create a graph object
g <- graph_from_data_frame(edges, directed = FALSE)

# Assign node types (Gene or Drug)
V(g)$type <- ifelse(V(g)$name %in% filtered_GDIdb_network_meaningful$gene_name, "Gene", "Drug")

# Assign interaction types to drugs
drug_interaction_mapping <- filtered_GDIdb_network_meaningful %>%
  distinct(drug_name, interaction_type) %>%
  group_by(drug_name) %>%
  summarise(interaction_type = first(interaction_type), .groups = "drop") # Keep first type if multiple

# Merge drug interaction types into graph nodes
V(g)$interaction_type <- ifelse(V(g)$name %in% drug_interaction_mapping$drug_name,
                                drug_interaction_mapping$interaction_type[match(V(g)$name, drug_interaction_mapping$drug_name)],
                                NA)

# Generate a color palette for different interaction types
interaction_types <- unique(na.omit(V(g)$interaction_type))
interaction_colors <- setNames(scales::hue_pal()(length(interaction_types)), interaction_types)

# Assign colors: Genes in "steelblue", Drugs colored by interaction type
node_colors <- ifelse(V(g)$type == "Gene", "steelblue", interaction_colors[V(g)$interaction_type])

# Compute node degree (number of connections per node)
node_degree <- degree(g)

# Define a threshold for labeling (only highly connected nodes)
label_threshold <- quantile(node_degree, 0.75)  # Top 25% most connected nodes

# Filter labels: Show only nodes with high degree
V(g)$label <- ifelse(node_degree >= label_threshold, V(g)$name, NA)

# Plot using ggraph with improved labeling
ggraph(g, layout = "fr") + 
  geom_edge_link(aes(edge_alpha = 0.5), color = "gray50") + 
  geom_node_point(aes(color = factor(V(g)$interaction_type, levels = interaction_types)), size = 5) + 
  geom_node_text(aes(label = label), repel = TRUE, size = 4, max.overlaps = 50) +  # Increased max.overlaps
  scale_color_manual(name = "Drug Effects", values = c("Gene" = "steelblue", interaction_colors)) +  
  theme_void() + 
  ggtitle("Drug-Gene Interaction Network: Colored by Interaction Type") + 
  theme(legend.position = "right")

#### 
#### 
#### Part J. Improving display
#### 
#### 
# Load required packages
library(igraph)
library(ggraph)
library(tidyverse)

# Filter meaningful interactions (excluding unknown interactions)
filtered_GDIdb_network_meaningful <- filtered_GDIdb %>%
  filter(interaction_type != "NULL" & interaction_type != "other/unknown") %>%
  select(gene_name, drug_name, interaction_type) %>%
  distinct()

# Create an edge list (interactions between genes and drugs)
edges <- filtered_GDIdb_network_meaningful %>%
  select(gene_name, drug_name) %>%
  distinct()

# Create a graph object
g <- graph_from_data_frame(edges, directed = FALSE)

# Assign node types (Gene or Drug)
V(g)$type <- ifelse(V(g)$name %in% filtered_GDIdb_network_meaningful$gene_name, "Gene", "Drug")

# Assign interaction types to drugs
drug_interaction_mapping <- filtered_GDIdb_network_meaningful %>%
  distinct(drug_name, interaction_type) %>%
  group_by(drug_name) %>%
  summarise(interaction_type = first(interaction_type), .groups = "drop") # Keep first type if multiple

# Merge drug interaction types into graph nodes
V(g)$interaction_type <- ifelse(V(g)$name %in% drug_interaction_mapping$drug_name,
                                drug_interaction_mapping$interaction_type[match(V(g)$name, drug_interaction_mapping$drug_name)],
                                NA)

# Generate a color palette for different interaction types
interaction_types <- unique(na.omit(V(g)$interaction_type))
interaction_colors <- setNames(scales::hue_pal()(length(interaction_types)), interaction_types)

# Assign colors: Genes in "steelblue", Drugs colored by interaction type
node_colors <- ifelse(V(g)$type == "Gene", "steelblue", interaction_colors[V(g)$interaction_type])

# Compute node degree (number of connections per node)
node_degree <- degree(g)

# Define a threshold for labeling (only the top 20% most connected nodes)
label_threshold <- quantile(node_degree, 0.80)  # Top 20% most connected nodes

# Filter labels: Show only nodes with high degree
V(g)$label <- ifelse(node_degree >= label_threshold, V(g)$name, NA)

# Plot using ggraph with improved labeling
ggraph(g, layout = "fr") + 
  geom_edge_link(aes(edge_alpha = 0.4), color = "gray50") + 
  geom_node_point(aes(color = factor(V(g)$interaction_type, levels = interaction_types)), size = 4) + 
  geom_node_label(aes(label = label), fill = "white", color = "black", size = 3.5, repel = TRUE, max.overlaps = 20, label.padding = unit(0.15, "lines")) + 
  scale_color_manual(name = "Drug Effects", values = c("Gene" = "steelblue", interaction_colors)) +  
  theme_void() + 
  ggtitle("Drug-Gene Interaction Network: Colored by Interaction Type") + 
  theme(legend.position = "right")

####  
####  
####  Interactive Network
####  
####  
# Load required packages
library(igraph)
library(visNetwork)
library(tidyverse)
library(htmlwidgets)  # Required to save HTML output

# Filter meaningful interactions (excluding unknown interactions)
filtered_GDIdb_network_meaningful <- filtered_GDIdb %>%
  filter(interaction_type != "NULL" & interaction_type != "other/unknown") %>%
  select(gene_name, drug_name, interaction_type) %>%
  distinct()

# Create an edge list (interactions between genes and drugs)
edges <- filtered_GDIdb_network_meaningful %>%
  select(gene_name, drug_name) %>%
  distinct()

# Create a graph object
g <- graph_from_data_frame(edges, directed = FALSE)

# Assign node types (Gene or Drug)
V(g)$type <- ifelse(V(g)$name %in% filtered_GDIdb_network_meaningful$gene_name, "Gene", "Drug")

# Assign interaction types to drugs
drug_interaction_mapping <- filtered_GDIdb_network_meaningful %>%
  distinct(drug_name, interaction_type) %>%
  group_by(drug_name) %>%
  summarise(interaction_type = first(interaction_type), .groups = "drop") # Keep first type if multiple

# Merge drug interaction types into graph nodes
V(g)$interaction_type <- ifelse(V(g)$name %in% drug_interaction_mapping$drug_name,
                                drug_interaction_mapping$interaction_type[match(V(g)$name, drug_interaction_mapping$drug_name)],
                                NA)

# Generate a color palette for different interaction types
interaction_types <- unique(na.omit(V(g)$interaction_type))
interaction_colors <- setNames(scales::hue_pal()(length(interaction_types)), interaction_types)

# Assign colors: Genes in "steelblue", Drugs colored by interaction type
V(g)$color <- ifelse(V(g)$type == "Gene", "steelblue", interaction_colors[V(g)$interaction_type])

# Convert graph into a data frame format for visNetwork
nodes <- data.frame(
  id = V(g)$name,
  label = V(g)$name,
  color = V(g)$color,
  group = V(g)$type,
  title = ifelse(V(g)$type == "Gene",
                 paste0("<b>Gene:</b> ", V(g)$name),
                 paste0("<b>Drug:</b> ", V(g)$name, "<br><b>Effect:</b> ", V(g)$interaction_type))
)

edges <- data.frame(
  from = edges$gene_name,
  to = edges$drug_name
)

# Define legend for Gene and Drug categories
legend_nodes <- data.frame(
  label = c("Gene", "Drug"),
  shape = c("dot", "dot"),
  color = c("steelblue", "gray"),
  id = c("Gene", "Drug")
)

# Create interactive network visualization
network <- visNetwork(nodes, edges) %>%
  visEdges(smooth = FALSE) %>%
  visNodes(size = 10, font = list(size = 20)) %>%
  visInteraction(navigationButtons = TRUE, tooltipDelay = 50) %>%
  visGroups(groupname = "Gene", color = "steelblue") %>%
  visGroups(groupname = "Drug", color = "gray") %>%
  visLegend(addNodes = legend_nodes) %>%  # Corrected legend format
  visPhysics(enabled = TRUE) %>% 
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) 

# Save as HTML
saveWidget(network, "Drug_Gene_Interaction_Network.html", selfcontained = TRUE)

# Display the network in RStudio Viewer
network

#####
#####
#####
# Load required packages
library(igraph)
library(visNetwork)
library(tidyverse)
library(htmlwidgets)  # Required to save HTML output

# Filter meaningful interactions (excluding unknown interactions)
filtered_GDIdb_network_meaningful <- filtered_GDIdb %>%
  filter(interaction_type != "NULL" & interaction_type != "other/unknown") %>%
  select(gene_name, drug_name, interaction_type) %>%
  distinct()

# Create an edge list (interactions between genes and drugs)
edges <- filtered_GDIdb_network_meaningful %>%
  select(gene_name, drug_name) %>%
  distinct()

# Create a graph object
g <- graph_from_data_frame(edges, directed = FALSE)

# Assign node types (Gene or Drug)
V(g)$type <- ifelse(V(g)$name %in% filtered_GDIdb_network_meaningful$gene_name, "Gene", "Drug")

# Assign interaction types to drugs
drug_interaction_mapping <- filtered_GDIdb_network_meaningful %>%
  distinct(drug_name, interaction_type) %>%
  group_by(drug_name) %>%
  summarise(interaction_type = first(interaction_type), .groups = "drop") # Keep first type if multiple

# Merge drug interaction types into graph nodes
V(g)$interaction_type <- ifelse(V(g)$name %in% drug_interaction_mapping$drug_name,
                                drug_interaction_mapping$interaction_type[match(V(g)$name, drug_interaction_mapping$drug_name)],
                                NA)

# Generate a distinct color palette for different interaction types
interaction_types <- unique(na.omit(V(g)$interaction_type))
interaction_colors <- setNames(scales::hue_pal()(length(interaction_types)), interaction_types)

# Assign colors: Genes in "steelblue", Drugs colored by interaction type
V(g)$color <- ifelse(V(g)$type == "Gene", "steelblue", interaction_colors[V(g)$interaction_type])

# Convert graph into a data frame format for visNetwork
nodes <- data.frame(
  id = V(g)$name,
  label = V(g)$name,
  color = V(g)$color,
  group = ifelse(V(g)$type == "Gene", "Gene", V(g)$interaction_type), # Use drug effect type as group
  title = ifelse(V(g)$type == "Gene",
                 paste0("<b>Gene:</b> ", V(g)$name),
                 paste0("<b>Drug:</b> ", V(g)$name, "<br><b>Effect:</b> ", V(g)$interaction_type))
)

edges <- data.frame(
  from = edges$gene_name,
  to = edges$drug_name
)

# Initialize the network
network <- visNetwork(nodes, edges) %>%
  visEdges(smooth = FALSE) %>%
  visNodes(size = 10, font = list(size = 20)) %>%
  visInteraction(navigationButtons = TRUE, tooltipDelay = 50)

# Add the "Gene" group
network <- network %>% visGroups(groupname = "Gene", color = "steelblue")

# Add each drug effect type separately in a loop
for (interaction in interaction_types) {
  network <- network %>% visGroups(groupname = interaction, color = interaction_colors[[interaction]])
}

# Add the correctly formatted legend
network <- network %>%
  visLegend(useGroups = TRUE, main = "Drug Effect Categorization") %>%
  visPhysics(enabled = TRUE) %>% 
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) 

# Save as HTML
saveWidget(network, "Drug_Gene_Interaction_Network.html", selfcontained = TRUE)

# Display the network in RStudio Viewer
network

####
####
#### Coverage of Target genes in signature database without duplicates
#### How many Target genes are components of the signatures database?
####
####
##### PART A
df1157 <- data_df1011

library(dplyr)

# List of variables to convert to numeric
numeric_vars <- c("Ranking", "Count_source", "GSI", "GFC", "PFC", "TNC", "HRC", "SMC", "TMC", "TIC", 
                  "RCD", "Protective_Count", "Risky_Count", "NS_Count", "GFC_ranking", 
                  "PFC_ranking", "SCS_ranking", "TNC_ranking", "HRC_ranking", "SMC_ranking", "TMC_ranking", 
                  "TIC_ranking", "RCD_ranking", "Final_Rank", "Cox_DSS_id", "Cox_DFI_id", "Cox_PFI_id", "Cox_OS_id",
                  "log_Count_source")

# Convert specified variables to numeric
df1157 <- df1157 %>%
  mutate(across(all_of(numeric_vars), as.numeric))

# Display the structure to verify changes
str(df1157)

# Check for duplicate rows across all columns
duplicates_exist <- any(duplicated(df1157))

if (duplicates_exist) {
  cat("Duplicate rows found. Removing duplicates...\n")
  
  # Remove duplicate rows while retaining one instance
  df1157 <-df1157 %>%
    distinct()
  
  cat("Duplicates removed. Updated dataframe:\n")
  print(df1157)
} else {
  cat("No duplicate rows found.\n")
}

# Check for duplicate rows across all columns
duplicate_indices <- which(duplicated(df1157) | duplicated(df1157, fromLast = TRUE))  # Identify all duplicate rows

if (length(duplicate_indices) > 0) {
  cat("Duplicate rows found. Printing duplicate rows:\n")
  
  # Print the duplicated rows
  duplicated_rows <-df1157[duplicate_indices, ]
  print(duplicated_rows)
  
  # Remove duplicate rows while retaining one instance
  df1157 <-df1157 %>%
    distinct()
  
  cat("Duplicates removed. Updated dataframe:\n")
  print(df1157)
} else {
  cat("No duplicate rows found.\n")
}

# Rename, reorder, and sort the columns in descending order by Rank
df1157 <-df1157 %>%
  # Reorder columns
  arrange(desc(Ranking))  # Sort in descending order by Rank

# Ensure "Count_source" is numeric before computing statistics
df1157 <- df1157 %>%
  mutate(Count_source = as.numeric(Count_source))

# Compute full statistics for "Count_source"
stats_Count_source <- df1157 %>%
  summarise(
    Mean = mean(Count_source, na.rm = TRUE),
    Median = median(Count_source, na.rm = TRUE),
    Q1 = quantile(Count_source, 0.25, na.rm = TRUE),  # First quartile (25th percentile)
    Q3 = quantile(Count_source, 0.75, na.rm = TRUE),  # Third quartile (75th percentile)
    P10 = quantile(Count_source, 0.10, na.rm = TRUE), # 10th percentile
    P90 = quantile(Count_source, 0.90, na.rm = TRUE), # 90th percentile
    Min = min(Count_source, na.rm = TRUE),
    Max = max(Count_source, na.rm = TRUE),
    Std_Dev = sd(Count_source, na.rm = TRUE),
    Variance = var(Count_source, na.rm = TRUE)
  )

# Display results
print(stats_Count_source)


df1157 <-df1157 %>%
  mutate_all(trimws)

saveRDS(df1157, "Dataset_S2_full.rds")

Dataset_S2_full <-  df1157

# Rename variables and relocate "Omic feature" after "Elements"
Dataset_S2_full <- Dataset_S2_full %>%
  rename(
    Rank = Ranking,          # Rename "Ranking" to "Ranlk"
    Elements = Count_source,  # Rename "Count_source" to "Elements"
    `RCD form` = RCD_types,   # Rename "RCD_types" to "RCD form"
    `Omic feature` = Genotype # Rename "Genotype" to "Omic feature"
  ) %>%
  relocate(`Omic feature`, .after = Elements)  # Relocate "Omic feature" after "Elements"

# Step 1: Convert "Elements" to numeric
Dataset_S2_full <- Dataset_S2_full %>%
  mutate(Elements = as.numeric(Elements))

# Optional: Summary of basic statistics (includes min, 1st quartile, median, mean, 3rd quartile, and max)
summary_stats <- summary(Dataset_S2_full$Elements)
print(summary_stats)

# Generate full descriptive statistics including centiles for "Elements"
Elements_stats <- Dataset_S2_full %>%
  summarise(
    Min = min(Elements, na.rm = TRUE),
    Max = max(Elements, na.rm = TRUE),
    Mean = mean(Elements, na.rm = TRUE),
    Median = median(Elements, na.rm = TRUE),
    Std_Dev = sd(Elements, na.rm = TRUE),
    Variance = var(Elements, na.rm = TRUE),
    Q1 = quantile(Elements, 0.25, na.rm = TRUE),  # 25th percentile
    Q3 = quantile(Elements, 0.75, na.rm = TRUE),  # 75th percentile
    IQR = IQR(Elements, na.rm = TRUE),  # Interquartile Range
    P1 = quantile(Elements, 0.01, na.rm = TRUE),  # 1st percentile
    P5 = quantile(Elements, 0.05, na.rm = TRUE),  # 5th percentile
    P10 = quantile(Elements, 0.10, na.rm = TRUE), # 10th percentile
    P90 = quantile(Elements, 0.90, na.rm = TRUE), # 90th percentile
    P95 = quantile(Elements, 0.95, na.rm = TRUE), # 95th percentile
    P99 = quantile(Elements, 0.99, na.rm = TRUE), # 99th percentile
    Missing_Values = sum(is.na(Elements))
  )

rio::export(df1157, "Dataset_S2_full.tsv")

#### Renaming variables in df1157 ####
df1157 <- df1157 %>%
  rename(
    Rank = Ranking,
    Elements = Count_source,
    `RCD form` = RCD_types,
    `Omic feature` = Genotype
  )

# Print confirmation message
message("Variables renamed successfully in df1157!")

####
#### sub-setting df1157 by `Omic feature` values to created Omic feature-specifis dataset1
####

# Step 1: Get unique "Omic feature" values
unique_genotypes_df1157 <- unique(df1157$`Omic feature`)

# Step 2: Loop through each unique Omic feature to create and store separate datasets
for (genotype in unique_genotypes_df1157) {
  # Create a filtered dataframe for the current Omic feature
  genotype_df <- df1157 %>% filter(`Omic feature` == genotype)
  
  # Ensure the dataframe is properly stored in the global environment
  assign(paste0("Dataset_S1_", genotype), genotype_df, envir = .GlobalEnv)
  
  # Define the .tsv filename
  file_name_tsv  <- paste0("Dataset_S1_", genotype, ".tsv")
  
  # Save dataframe as .tsv only
  rio::export(genotype_df, file_name_tsv)
  
  # Print confirmation message
  message("File created: ", file_name_tsv, " | Object available in environment: df1157_", genotype)
}

# Print final message
message("All genotype-specific objects are now in the RStudio environment and saved as .tsv files!")

#####
#####
##### Stats of Dataset_S1_mRNA
#####
##### 

# Rename variables and relocate "Omic feature" after "Elements"
Dataset_S1_mRNA_full <- Dataset_S1_mRNA %>%
  relocate(`Omic feature`, .after = Elements) %>%
  mutate(Elements = as.numeric(Elements)) # Relocate "Omic feature" after "Elements"

# Optional: Summary of basic statistics (includes min, 1st quartile, median, mean, 3rd quartile, and max)
summary_stats <- summary(Dataset_S1_mRNA_full$Elements)
print(summary_stats)

# Generate full descriptive statistics including centiles for "Elements"
Elements_stats <- Dataset_S1_mRNA_full %>%
  summarise(
    Min = min(Elements, na.rm = TRUE),
    Max = max(Elements, na.rm = TRUE),
    Mean = mean(Elements, na.rm = TRUE),
    Median = median(Elements, na.rm = TRUE),
    Std_Dev = sd(Elements, na.rm = TRUE),
    Variance = var(Elements, na.rm = TRUE),
    Q1 = quantile(Elements, 0.25, na.rm = TRUE),  # 25th percentile
    Q3 = quantile(Elements, 0.75, na.rm = TRUE),  # 75th percentile
    IQR = IQR(Elements, na.rm = TRUE),  # Interquartile Range
    P1 = quantile(Elements, 0.01, na.rm = TRUE),  # 1st percentile
    P5 = quantile(Elements, 0.05, na.rm = TRUE),  # 5th percentile
    P10 = quantile(Elements, 0.10, na.rm = TRUE), # 10th percentile
    P90 = quantile(Elements, 0.90, na.rm = TRUE), # 90th percentile
    P95 = quantile(Elements, 0.95, na.rm = TRUE), # 95th percentile
    P99 = quantile(Elements, 0.99, na.rm = TRUE), # 99th percentile
    Missing_Values = sum(is.na(Elements))
  )

print(Elements_stats)

####
####
####
# sub-setting Dataset_S2_full by `Omic feature` values to created `Omic feature`-specific dataset1
# subset for CancerRCDShiby application downloads
library(dplyr)
library(openxlsx)
library(readr)

# Step 1: Get unique "`Omic feature`" values
unique_genotypes_Dataset_S2_full <- unique(Dataset_S2_full$`Omic feature`)

# Step 2: Loop through each unique `Omic feature` to create and save separate datasets
for (genotype in unique_genotypes_Dataset_S2_full) {
  # Create a filtered dataframe for the current `Omic feature`
  genotype_df <- Dataset_S2_full %>% filter(`Omic feature` == genotype)
  
  # Assign dataframe to a dynamically created object in the environment
  assign(paste0("Dataset_S2_full_", genotype), genotype_df, envir = .GlobalEnv)
  
}

rio::export(Dataset_S2_full, "Dataset_S2_full.tsv")
rio::export(Dataset_S2_full_CNV, "Dataset_S2_full_CNV.tsv")
rio::export(Dataset_S2_full_Methylation, "Dataset_S2_full_Methylation.tsv")
rio::export(Dataset_S2_full_miRNA, "Dataset_S2_full_miRNA.tsv")
rio::export(Dataset_S2_full_mRNA, "Dataset_S2_full_mRNA.tsv")
rio::export(Dataset_S2_full_Mutation, "Dataset_S2_full_Mutation.tsv")
rio::export(Dataset_S2_full_Protein, "Dataset_S2_full_Protein.tsv")
rio::export(Dataset_S2_full_Transcript, "Dataset_S2_full_Transcript.tsv")

# Print final message
print("All genotype-specific objects created successfully!")

#######################
#######################
##### Analysis of top signatures ####
#######################
#######################
setwd("D:/Pré-artigo 5-optosis model/Dataframes all metrics")

data <- import("Drafts_2/gene_info_trancritos_excel2.xlsx")

data <- data %>%
  mutate_all(trimws)

miRNA_database <-  import("Drafts_2/miRNA_database.xlsx")

library(dplyr)
library(stringr)

# Rename "hsa-mir-" to "hsa-miR-" in the "Mature1" variable
miRNA_database <- miRNA_database %>%
  mutate(Mature1 = str_replace(Mature1, "^hsa-mir-", "hsa-miR-"))

# Display unique values to verify the change
unique(miRNA_database$Mature1)

# Display unique values to verify the change
unique(miRNA_database$Mature2)

miRNA_database <- miRNA_database %>%
  mutate_all(trimws)

protein_info <- import("Drafts_2/protein_info.xlsx") 

protein_info <- protein_info %>%
  mutate_all(trimws)

df1157 <-df1157 %>%
  mutate_all(trimws)

data_df1011 <- data_df1011 %>%
  mutate_all(trimws)

# Step 0: Filter data_df1011 based on matching values in the "Nomenclature" column of data_df1011 and "Nomenclature" column indf1157
data_df1011_filtered <- data_df1011 %>%
  filter(Nomenclature %in%df1157$Nomenclature)

# Step 1: Copy the filtered columns "Signature" and "Nomenclature" from data_df1011_filtered to df40
df40 <- data_df1011_filtered %>%
  select(Signature, Nomenclature)

# Step 2: Create a new column "Gene_conversion" after "Signature"
df40 <- df40 %>%
  mutate(Gene_conversion = NA)  # Initialize the column with NA values

# Step 3: Populate "Gene_conversion" using left_join for better matching reliability

# For values starting with "ENST", join with the "data" dataframe on "Transcript_ID"
df40 <- df40 %>%
  left_join(data %>% select(Transcript_ID, Display_Name), by = c("Signature" = "Transcript_ID")) %>%
  mutate(Gene_conversion = ifelse(!is.na(Display_Name), Display_Name, Gene_conversion)) %>%
  select(-Display_Name)  # Remove the temporary column

# For values starting with "hsa-miR", join with the "miRNA_database" on "Mature1" and "Mature2"
df40 <- df40 %>%
  left_join(miRNA_database %>% select(Mature1, Gene_Symbol), by = c("Signature" = "Mature1")) %>%
  mutate(Gene_conversion = ifelse(!is.na(Gene_Symbol), Gene_Symbol, Gene_conversion)) %>%
  select(-Gene_Symbol)  # Remove the temporary column

df40 <- df40 %>%
  left_join(miRNA_database %>% select(Mature2, Gene_Symbol), by = c("Signature" = "Mature2")) %>%
  mutate(Gene_conversion = ifelse(!is.na(Gene_Symbol), Gene_Symbol, Gene_conversion)) %>%
  select(-Gene_Symbol)  # Remove the temporary column

# For protein values join with the "protein_info" on "Protein1" 
df40 <- df40 %>%
  left_join(protein_info %>% select(Protein1, Gene_Symbol), by = c("Signature" = "Protein1")) %>%
  mutate(Gene_conversion = ifelse(!is.na(Gene_Symbol), Gene_Symbol, Gene_conversion)) %>%
  select(-Gene_Symbol)  # Remove the temporary column

# Step 1: Extract gene symbols, transcript names, and miRNA from the "Signature" column
df41 <- df40 %>%
  # Use str_replace_all to remove parentheses and separate by "+" into individual elements
  mutate(Signature_cleaned = str_replace_all(Signature, "[()]", "")) %>%  # Remove parentheses
  separate_rows(Signature_cleaned, sep = " \\+ ") %>%  # Split by " + " into separate rows
  distinct() %>%  # Remove duplicate rows
  filter(str_detect(Signature_cleaned, "^ENST|^hsa-miR|^[A-Za-z]")) %>%  # Keep only gene symbols, ENST transcripts, or miRNA
  arrange(Signature_cleaned)  # Order alphabetically

# Step 2: Create the df41 with the unique and ordered entries
df41 <- df41 %>%
  select(Signature_cleaned) %>%
  rename(Extracted_Signature = Signature_cleaned)

# Step 1: Create a new column "Gene_Symbol" in df41 and populate it based on the mapping logic
df41 <- df41 %>%
  rowwise() %>%
  mutate(Gene_Symbol = case_when(
    
    # For values with prefix "ENST", look up in the "data" dataframe using "Transcript_ID"
    str_detect(Extracted_Signature, "^ENST") ~ 
      ifelse(Extracted_Signature %in% data$Transcript_ID, 
             data$Display_Name[which(data$Transcript_ID == Extracted_Signature)], 
             NA_character_),
    
    # For values with prefix "hsa-miR", look up in "miRNA_database" for Mature1
    str_detect(Extracted_Signature, "^hsa-miR") & Extracted_Signature %in% miRNA_database$Mature1 ~ 
      ifelse(Extracted_Signature %in% miRNA_database$Mature1, 
             miRNA_database$Gene_Symbol[which(miRNA_database$Mature1 == Extracted_Signature)], 
             NA_character_),
    
    # For values with prefix "hsa-miR", look up in "miRNA_database" for Mature2
    str_detect(Extracted_Signature, "^hsa-miR") & Extracted_Signature %in% miRNA_database$Mature2 ~ 
      ifelse(Extracted_Signature %in% miRNA_database$Mature2, 
             miRNA_database$Gene_Symbol[which(miRNA_database$Mature2 == Extracted_Signature)], 
             NA_character_),
    
    # For values with prefix "hsa-let", look up in "miRNA_database" for Mature1
    str_detect(Extracted_Signature, "^hsa-let") & Extracted_Signature %in% miRNA_database$Mature1 ~ 
      ifelse(Extracted_Signature %in% miRNA_database$Mature1, 
             miRNA_database$Gene_Symbol[which(miRNA_database$Mature1 == Extracted_Signature)], 
             NA_character_),
    
    # For values with prefix "hsa-let", look up in "miRNA_database" for Mature2
    str_detect(Extracted_Signature, "^hsa-let") & Extracted_Signature %in% miRNA_database$Mature2 ~ 
      ifelse(Extracted_Signature %in% miRNA_database$Mature2, 
             miRNA_database$Gene_Symbol[which(miRNA_database$Mature2 == Extracted_Signature)], 
             NA_character_),
    
    # If no match is found, leave as NA
    TRUE ~ NA_character_
  ))

# Fill NA values in "Gene_Symbol" with corresponding values from "Extracted_Signature"
df41 <- df41 %>%
  mutate(Gene_Symbol = ifelse(is.na(Gene_Symbol), Extracted_Signature, Gene_Symbol))

# Remove duplicates based on "Gene_Symbol"
df41 <- df41 %>%
  distinct(Gene_Symbol, .keep_all = TRUE)

# Order the dataframe by "Gene_Symbol" alphabetically
df41 <- df41 %>%
  arrange(Gene_Symbol)

# Step 1: Identify the common values between "Protein1" and "Gene_symbol"
common_proteins <- intersect(protein_info$Protein1, df41$Gene_Symbol)

# Step 2: Create a lookup table for replacements
protein_lookup <- protein_info %>%
  filter(Protein1 %in% common_proteins) %>%
  select(Protein1) %>%
  distinct()

# Step 3: Replace the "Gene_symbol" values in df41 with the corresponding "Protein1" values from protein_info
df41 <- df41 %>%
  mutate(Gene_Symbol = ifelse(Gene_Symbol %in% protein_info$Protein1,  # Check if Gene_symbol is in Protein1
                              protein_info$Gene_Symbol[match(Gene_Symbol, protein_info$Protein1)],  # Replace with corresponding Gene_symbol in protein_info
                              Gene_Symbol))  # Keep original Gene_symbol if no match

# Exclude the "Extracted_Signature" column
df41 <- df41 %>%
  select(-Extracted_Signature)

# Remove leading and trailing spaces globally across all columns
df41 <- df41 %>%
  mutate_all(trimws)

#  Remove duplicate rows
duplicates_exist <- any(duplicated(df41))

if (duplicates_exist) {
  cat("Duplicate rows found. Removing duplicates...\n")
  
  # Remove duplicate rows while retaining one instance
  df41 <- df41 %>%
    distinct()
  
  cat("Duplicates removed. Updated dataframe:\n")
  print(df41)
} else {
  cat("No duplicate rows found in the dataframe.\n")
}

# Assign dataframe to Gene_Symbols_signatures_df1157
Gene_Symbols_signatures_df1157 <- df41

# Order the dataframe alphabetically by the "Gene_Symbol" column
Gene_Symbols_signatures_df1157 <- Gene_Symbols_signatures_df1157[order(Gene_Symbols_signatures_df1157$Gene_Symbol), ]

# Remove leading and trailing spaces globally across all columns
Gene_Symbols_signatures_df1157 <- Gene_Symbols_signatures_df1157 %>%
  mutate_all(trimws)

#  Remove duplicate rows
duplicates_exist <- any(duplicated(Gene_Symbols_signatures_df1157))

if (duplicates_exist) {
  cat("Duplicate rows found. Removing duplicates...\n")
  
  # Remove duplicate rows while retaining one instance
  Gene_Symbols_signatures_df1157 <- Gene_Symbols_signatures_df1157 %>%
    distinct()
  
  cat("Duplicates removed. Updated dataframe:\n")
  print(Gene_Symbols_signatures_df1157)
} else {
  cat("No duplicate rows found in the dataframe.\n")
}

# Ensure the dataframes have been loaded or created: Target_genes and df41

# Filter Target_genes based on matching values in Gene_Symbol from df41
filtered_Target_genes <- Target_genes[Target_genes$Gene_Symbol %in% df41$Gene_Symbol, ]

# Identify rows in df41 where Gene_Symbol is NOT in Target_genes$Gene_Symbol
non_matching_rows <- df41[!df41$Gene_Symbol %in% Target_genes$Gene_Symbol, ]

# View the non-matching rows
head(non_matching_rows)

#### Handling "PRKN", "PARK2" ambiguity (02/23/2025) 
#### This is so because: PARK2 is the gene Symbol that codes the Parkin protein, also known as PRKN 
# Trim whitespace, replace "PRKN" with "PARK2", and remove duplicates
Gene_Symbols_signatures_df1157 <- Gene_Symbols_signatures_df1157 %>%
  mutate(Gene_Symbol = trimws(Gene_Symbol)) %>%  # Trim whitespace
  mutate(Gene_Symbol = ifelse(Gene_Symbol == "PRKN", "PARK2", Gene_Symbol)) %>%  # Replace PRKN
  distinct()  # Remove duplicate rows

# Print confirmation message
print("Replaced 'PRKN' with 'PARK2', removed duplicates, and trimmed whitespace in 'Gene_Symbol' column.")

# Load required libraries
library(readr)    # For writing TSV files
library(openxlsx) # For writing Excel files

# Step 5: rio::export the cleaned dataframe to a TSV file
write_tsv(Gene_Symbols_signatures_df1157, "Gene_Symbols_signatures_df1157.tsv")

# Step 6: rio::export the cleaned dataframe to an Excel file
write.xlsx(Gene_Symbols_signatures_df1157, "Gene_Symbols_signatures_df1157.xlsx", overwrite = TRUE)


# Step 5: rio::export the cleaned dataframe to a TSV file
rio::export(Gene_Symbols_signatures_df1157, "Gene_Symbols_signatures_df1157.tsv")
rio::export(Gene_Symbols_signatures_df1157, "Gene_Symbols_signatures_df1157.xlsx")

cat("Cleaned dataframe rio::exported to 'Gene_Symbols_signatures_df1157.tsv'.\n")

####
####
#### Estimating signatures frequency by specific features and needs
####
####
####
library(dplyr)

# Define the correspondence table
TNC_to_Expression <- data.frame(
  Expression = c("No_data", "Unchanged", "Underexpression", "Overexpression"),  # Original TNC values
  TNC = c(0, 1, 2, 3)  # The corresponding IDs for each TNC value
)
# Ensure TNC in data_df1011 is numeric for proper joining
df1157$TNC <- as.numeric(df1157$TNC)

# Create the new "Expression" variable in data_df1011 based on the values in "TNC"
df1157 <- df1157 %>%
  left_join(TNC_to_Expression, by = "TNC") %>%  # Add the "Expression" column based on "TNC" values
  relocate(Expression, .after = Phenotype)       # Move "Expression" right after "Genotype"

# rio::export the modified dataframe
rio::export(data_df1011, "df1157.tsv")

# Filter rows where "RCD" is not equal to "1"
df1157_filtered_multimodular <- df1157 %>%
  filter(RCD != 1)

library(dplyr)

# Filter rows where "RCD" is not equal to "P"
df1157_filtered_multimodular_SCS_N <- df1157_filtered_multimodular %>%
  filter(SCS != "P")

# Filter rows where "RCD" is not equal to "N"
df1157_filtered_multimodular_SCS_P <- df1157_filtered_multimodular %>%
  filter(SCS != "N")

# Filter rows where "RCD" is not equal to "Overexpression"
df1157_filtered_multimodular_overexpression <- df1157_filtered_multimodular %>%
  filter(Expression == "Overexpression")

# Filter rows where "RCD" is equal to "1"
df1157_filtered_unimodular <- df1157 %>%
  filter(RCD == "1")

library(dplyr)

# Filter rows where "RCD" is not equal to "P"
df1157_filtered_multimodular_SCS_N <- df1157_filtered_multimodular %>%
  filter(SCS != "P")

# Filter rows where "RCD" is not equal to "N"
df1157_filtered_multimodular_SCS_P <- df1157_filtered_multimodular %>%
  filter(SCS != "N")

# Filter rows where "RCD" is not equal to "Overexpression"
df1157_filtered_multimodular_overexpression <- df1157_filtered_multimodular %>%
  filter(Expression == "Overexpression")

######
######
###### Global analysis of dataframes by Omic feature ######
###### 
###### 
library(dplyr)

# Define total row count of the original dataset dynamically
total_rows <- nrow(df1157_mRNA)

# Step 1: Filter rows where "SCS" is equal to "N"
df1157_mRNA_N <- df1157_mRNA %>%
  filter(SCS == "N")
n_N <- nrow(df1157_mRNA_N)  # Store row count
percent_N <- round((n_N / total_rows) * 100, 2)  # Percentage relative to df1157_mRNA_N

# Step 2: Further filter the subset where "Phenotype" is equal to "TSM"
df1157_mRNA_N_TSM <- df1157_mRNA_N %>%
  filter(Phenotype == "TSM")
n_TSM <- nrow(df1157_mRNA_N_TSM)  # Store row count
percent_TSM <- round((n_TSM / n_N) * 100, 2)  # Percentage relative to df1157_mRNA_N

# Step 3: Further filter rows that contain "Risky" in any column
df1157_mRNA_N_TSM_Risky <- df1157_mRNA_N_TSM %>%
  filter(if_any(everything(), ~ grepl("Risky", .)))
n_Risky <- nrow(df1157_mRNA_N_TSM_Risky)  # Store row count
percent_Risky <- round((n_Risky / n_TSM) * 100, 2)  # Percentage relative to df1157_mRNA_N_TSM

# Step 4: Further filter rows that contain "Protective" in any column
df1157_mRNA_N_TSM_Protective <- df1157_mRNA_N_TSM %>%
  filter(if_any(everything(), ~ grepl("Protective", .)))
n_Protective <- nrow(df1157_mRNA_N_TSM_Protective)  # Store row count
percent_Protective <- round((n_Protective / n_TSM) * 100, 2)  # Percentage relative to df1157_mRNA_N_TSM

# Display summaries of the new dataframes with percentages
# Display summaries of the new dataframes with percentages
print(paste("Rows in df1157_mRNA:", total_rows))
print(paste("Rows in df1157_mRNA_N:", n_N, "(", percent_N, "% of df1157_mRNA)"))
print(paste("Rows in df1157_mRNA_N_TSM:", n_TSM, "(", percent_TSM, "% of df1157_mRNA_N)"))
print(paste("Rows in df1157_mRNA_N_TSM_Risky:", n_Risky, "(", percent_Risky, "% of df1157_mRNA_N_TSM)"))
print(paste("Rows in df1157_mRNA_N_TSM_Protective:", n_Protective, "(", percent_Protective, "% of df1157_mRNA_N_TSM)"))


#######
#######
####### PART A - df1157_mRNA

# Define total row count of the original dataset dynamically
total_rows <- nrow(df1157_mRNA)

# Function to filter by microenvironment_classification value and compute percentages
filter_microenvironment <- function(df, classification) {
  filtered_df <- df %>%
    filter(if_any(c(microenvironment_classification), ~ grepl(classification, ., ignore.case = TRUE)))
  
  n_filtered <- nrow(filtered_df)
  percent_filtered <- round((n_filtered / total_rows) * 100, 2)  # Percentage relative to df1157_mRNA
  
  # Print summary
  print(paste("Rows in df1157_mRNA_", classification, ":", n_filtered, "(", percent_filtered, "% of df1157_mRNA)", sep = ""))
  
  return(filtered_df)  # Return the filtered dataframe
}

# Apply filtering separately for each classification
df1157_mRNA_AntiTumoral <- filter_microenvironment(df1157_mRNA, "anti-tumoral")
df1157_mRNA_Dual <- filter_microenvironment(df1157_mRNA, "dual")
df1157_mRNA_ProTumoral <- filter_microenvironment(df1157_mRNA, "pro-tumoral")
df1157_mRNA_NS <- filter_microenvironment(df1157_mRNA, "NS")


library(dplyr)

###### Global analysis of dataframes by immune classification ######

# Define total row count of the original dataset dynamically
total_rows <- nrow(df1157_mRNA)

# Function to filter by immune_classification value and compute percentages
filter_immune_classification <- function(df, classification) {
  filtered_df <- df %>%
    filter(if_any(c(immune_classification), ~ grepl(classification, ., ignore.case = TRUE)))
  
  n_filtered <- nrow(filtered_df)
  percent_filtered <- round((n_filtered / total_rows) * 100, 2)  # Percentage relative to df1157_mRNA
  
  # Print summary
  print(paste("Rows in df1157_mRNA_", classification, ":", n_filtered, "(", percent_filtered, "% of df1157_mRNA)", sep = ""))
  
  return(filtered_df)  # Return the filtered dataframe
}

# Apply filtering separately for each classification
df1157_mRNA_Variable <- filter_immune_classification(df1157_mRNA, "Variable")
df1157_mRNA_Hot <- filter_immune_classification(df1157_mRNA, "Hot")
df1157_mRNA_NS <- filter_immune_classification(df1157_mRNA, "NS")
df1157_mRNA_Cold <- filter_immune_classification(df1157_mRNA, "Cold")

######
######
######
###### PART B - df1157_Transcript
######
######
######

#### Renaming variables in df1157 ####
df1157_Transcript <- df1157_Transcript %>%
  rename(
    Ranking = Rank,
    Count_source = Elements,
    RCD_types = `RCD form`,
    Genotype = `Omic feature`
  )

# Ensure "Count_source" is numeric before computing statistics
df1157_Transcript <- df1157_Transcript %>%
  mutate(Count_source = as.numeric(Count_source))

# Ensure "Risky_Count" is numeric before computing statistics
df1157_Transcript <- df1157_Transcript %>%
  mutate(Risky_Count = as.numeric(Risky_Count))

# Ensure "Count_source" is numeric before computing statistics
df1157_Transcript <- df1157_Transcript %>%
  mutate(Count_source = as.numeric(Count_source))

# Compute full statistics for "Count_source"
stats_Count_source <- df1157_Transcript %>%
  summarise(
    Mean = mean(Count_source, na.rm = TRUE),
    Median = median(Count_source, na.rm = TRUE),
    Q1 = quantile(Count_source, 0.25, na.rm = TRUE),  # First quartile (25th percentile)
    Q3 = quantile(Count_source, 0.75, na.rm = TRUE),  # Third quartile (75th percentile)
    P10 = quantile(Count_source, 0.10, na.rm = TRUE), # 10th percentile
    P90 = quantile(Count_source, 0.90, na.rm = TRUE), # 90th percentile
    Min = min(Count_source, na.rm = TRUE),
    Max = max(Count_source, na.rm = TRUE),
    Std_Dev = sd(Count_source, na.rm = TRUE),
    Variance = var(Count_source, na.rm = TRUE)
  )

# Display results
print(stats_Count_source)

# Define total row count of the original dataset dynamically
total_rows <- nrow(df1157_Transcript)

# Step 1: Further filter rows that contain "Risky" or "Protective" in any column
df1157_Transcript_SMC_R_P<- df1157_Transcript %>%
  filter(if_any(everything(), ~ grepl("Risky|Protective", ., ignore.case = TRUE)))  # Case-insensitive search for both terms

# Store row count
n_SMC_R_P <- nrow(df1157_Transcript_SMC_R_P)  

# Compute percentage relative to df1157_Transcript
percent_SMC_R_P <- round((n_SMC_R_P / nrow(df1157_Transcript)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Transcript_SMC_R_P:", n_SMC_R_P, "(", percent_SMC_R_P, "% of df1157_Transcript)"))

####
####
####
# Step 2: Further filter rows that contain "Risky" in all four metrics
df1157_Transcript_SMC_R_4<- df1157_Transcript_SMC_R_P %>%
  filter(Risky_Count == 4)  

# Store row count
n_SMC_R_4 <- nrow(df1157_Transcript_SMC_R_4)  

# Compute percentage relative to df1157_Transcript_SMC_R_4
percent_SMC_R_4 <- round((n_SMC_R_4 / nrow(df1157_Transcript_SMC_R_P)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Transcript_SMC_R_4:", n_SMC_R_4, "(", percent_SMC_R_4, "% of df1157_Transcript_SMC_R_P)"))

#####
#####
#####
# Step 3: Further filter rows that contain "Protective" in all four metrics
df1157_Transcript_SMC_P_4<- df1157_Transcript_SMC_R_P %>%
  filter(Protective_Count == 4)  

# Store row count
n_SMC_P_4 <- nrow(df1157_Transcript_SMC_P_4)  

# Compute percentage relative to df1157_Transcript_SMC_P_4
percent_SMC_P_4 <- round((n_SMC_P_4 / nrow(df1157_Transcript_SMC_R_P)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Transcript_SMC_P_4:", n_SMC_P_4, "(", percent_SMC_P_4, "% of df1157_Transcript_SMC_R_P)"))

#####
#####
#####
# Step 4: Further filter rows that contain "Risky_Count == 4" and "Phenotype == TSM"
df1157_Transcript_SMC_R_4_TSM <- df1157_Transcript_SMC_R_4 %>%
  filter(Phenotype == "TSM")  

# Store row count
n_R_4_TSM <- nrow(df1157_Transcript_SMC_R_4_TSM)  

# Compute percentage relative to df1157_Transcript_SMC_P_4
percent_R_4_TSM <- round((n_R_4_TSM / nrow(df1157_Transcript_SMC_R_4)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Transcript_SMC_R_4_TSM:", n_R_4_TSM, "(", percent_R_4_TSM, "% of df1157_Transcript_SMC_R_4)"))

#####
#####
#####
#### Step 5: Further filter rows that contain "Protective_Count == 4" and "Phenotype == TSM"
df1157_Transcript_SMC_P_4_TSM <- df1157_Transcript_SMC_P_4 %>%
  filter(Phenotype == "TSM")  

# Store row count
n_P_4_TSM <- nrow(df1157_Transcript_SMC_P_4_TSM)  

# Compute percentage relative to df1157_Transcript_SMC_P_4
percent_P_4_TSM <- round((n_P_4_TSM / nrow(df1157_Transcript_SMC_P_4)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Transcript_SMC_P_4_TSM:", n_P_4_TSM, "(", percent_P_4_TSM, "% of df1157_Transcript_SMC_P_4)"))

#####
#####
#####
# Step 6: Filter rows where "SCS" is equal to "N"
df1157_Transcript_N <- df1157_Transcript %>%
  filter(SCS == "N")
n_N <- nrow(df1157_Transcript_N)  # Store row count
percent_N <- round((n_N / total_rows) * 100, 2)  # Percentage relative to df1157_Transcript_N

# Step 7: Further filter the subset where "Phenotype" is equal to "TSM"
df1157_Transcript_N_TSM <- df1157_Transcript_N %>%
  filter(Phenotype == "TSM")
n_TSM <- nrow(df1157_Transcript_N_TSM)  # Store row count
percent_TSM <- round((n_TSM / n_N) * 100, 2)  # Percentage relative to df1157_Transcript_N

# Step 8: Further filter rows that contain "Risky" in any column
df1157_Transcript_N_TSM_Risky <- df1157_Transcript_N_TSM %>%
  filter(if_any(everything(), ~ grepl("Risky", .)))
n_Risky <- nrow(df1157_Transcript_N_TSM_Risky)  # Store row count
percent_Risky <- round((n_Risky / n_TSM) * 100, 2)  # Percentage relative to df1157_Transcript_N_TSM

# Step 9: Further filter rows that contain "Protective" in any column
df1157_Transcript_N_TSM_Protective <- df1157_Transcript_N_TSM %>%
  filter(if_any(everything(), ~ grepl("Protective", .)))
n_Protective <- nrow(df1157_Transcript_N_TSM_Protective)  # Store row count
percent_Protective <- round((n_Protective / n_TSM) * 100, 2)  # Percentage relative to df1157_Transcript_N_TSM

# Display summaries of the new dataframes with percentages
# Display summaries of the new dataframes with percentages
print(paste("Rows in df1157_Transcript:", total_rows))
print(paste("Rows in df1157_Transcript_N:", n_N, "(", percent_N, "% of df1157_Transcript)"))
print(paste("Rows in df1157_Transcript_N_TSM:", n_TSM, "(", percent_TSM, "% of df1157_Transcript_N)"))
print(paste("Rows in df1157_Transcript_N_TSM_Risky:", n_Risky, "(", percent_Risky, "% of df1157_Transcript_N_TSM)"))
print(paste("Rows in df1157_Transcript_N_TSM_Protective:", n_Protective, "(", percent_Protective, "% of df1157_Transcript_N_TSM)"))

#######
#######
#######

# Define total row count of the original dataset dynamically
total_rows <- nrow(df1157_Transcript)

# Function to filter by microenvironment_classification value and compute percentages
filter_microenvironment <- function(df, classification) {
  filtered_df <- df %>%
    filter(if_any(c(microenvironment_classification), ~ grepl(classification, ., ignore.case = TRUE)))
  
  n_filtered <- nrow(filtered_df)
  percent_filtered <- round((n_filtered / total_rows) * 100, 2)  # Percentage relative to df1157_Transcript
  
  # Print summary
  print(paste("Rows in df1157_Transcript_", classification, ":", n_filtered, "(", percent_filtered, "% of df1157_Transcript)", sep = ""))
  
  return(filtered_df)  # Return the filtered dataframe
}

# Apply filtering separately for each classification
df1157_Transcript_AntiTumoral <- filter_microenvironment(df1157_Transcript, "anti-tumoral")
df1157_Transcript_Dual <- filter_microenvironment(df1157_Transcript, "dual")
df1157_Transcript_ProTumoral <- filter_microenvironment(df1157_Transcript, "pro-tumoral")
df1157_Transcript_NS <- filter_microenvironment(df1157_Transcript, "NS")

###### Global analysis of dataframes by immune classification ######

# Define total row count of the original dataset dynamically
total_rows <- nrow(df1157_Transcript)

# Function to filter by immune_classification value and compute percentages
filter_immune_classification <- function(df, classification) {
  filtered_df <- df %>%
    filter(if_any(c(immune_classification), ~ grepl(classification, ., ignore.case = TRUE)))
  
  n_filtered <- nrow(filtered_df)
  percent_filtered <- round((n_filtered / total_rows) * 100, 2)  # Percentage relative to df1157_Transcript
  
  # Print summary
  print(paste("Rows in df1157_Transcript_", classification, ":", n_filtered, "(", percent_filtered, "% of df1157_Transcript)", sep = ""))
  
  return(filtered_df)  # Return the filtered dataframe
}

# Apply filtering separately for each classification
df1157_Transcript_Variable <- filter_immune_classification(df1157_Transcript, "Variable")
df1157_Transcript_Hot <- filter_immune_classification(df1157_Transcript, "Hot")
df1157_Transcript_NS <- filter_immune_classification(df1157_Transcript, "NS")
df1157_Transcript_Cold <- filter_immune_classification(df1157_Transcript, "Cold")

####
####
###
# Step 10: Further filter rows that contain "Risky" and "Overexpression" in all four metrics
df1157_Transcript_SMC_R_4_Over<- df1157_Transcript_SMC_R_4 %>%
  filter(Expression == "Overexpression")  

# Store row count
n_SMC_R_4_Over <- nrow(df1157_Transcript_SMC_R_4_Over)  

# Compute percentage relative to df1157_Transcript_SMC_R_4
percent_SMC_R_4_Over <- round((n_SMC_R_4_Over / nrow(df1157_Transcript_SMC_R_4)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Transcript_SMC_R_4_over:", n_SMC_R_4_Over, "(", percent_SMC_R_4_Over, "% of df1157_Transcript_SMC_R_4)"))

####
####
####
# Step 11: Further filter rows that contain "Protective" and "Underexpression" in all four metrics
df1157_Transcript_SMC_P_4_Under<- df1157_Transcript_SMC_P_4 %>%
  filter(Expression == "Underexpression")  

# Store row count
n_SMC_P_4_Under <- nrow(df1157_Transcript_SMC_P_4_Under)  

# Compute percentage relative to df1157_Transcript_SMC_P_4
percent_SMC_P_4_Under <- round((n_SMC_P_4_Under / nrow(df1157_Transcript_SMC_P_4)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Transcript_SMC_P_4_Under:", n_SMC_P_4_Under, "(", percent_SMC_P_4_Under, "% of df1157_Transcript_SMC_P_4)"))

#####
#####
#####
# Step 12 Stats
# Ensure "Count_source" is numeric before computing statistics
df1157_Transcript_SMC_R_4 <- df1157_Transcript_SMC_R_4 %>%
  mutate(Count_source = as.numeric(Count_source))

# Compute full statistics for "Count_source"
stats_Count_source <- df1157_Transcript_SMC_R_4 %>%
  summarise(
    Mean = mean(Count_source, na.rm = TRUE),
    Median = median(Count_source, na.rm = TRUE),
    Q1 = quantile(Count_source, 0.25, na.rm = TRUE),  # First quartile (25th percentile)
    Q3 = quantile(Count_source, 0.75, na.rm = TRUE),  # Third quartile (75th percentile)
    P10 = quantile(Count_source, 0.10, na.rm = TRUE), # 10th percentile
    P90 = quantile(Count_source, 0.90, na.rm = TRUE), # 90th percentile
    Min = min(Count_source, na.rm = TRUE),
    Max = max(Count_source, na.rm = TRUE),
    Std_Dev = sd(Count_source, na.rm = TRUE),
    Variance = var(Count_source, na.rm = TRUE)
  )

# Display results
print(stats_Count_source)

#####
#####
#####
# Step 13 Stats
# Ensure "Count_source" is numeric before computing statistics
df1157_Transcript_SMC_P_4 <- df1157_Transcript_SMC_P_4 %>%
  mutate(Count_source = as.numeric(Count_source))

# Compute full statistics for "Count_source"
stats_Count_source <- df1157_Transcript_SMC_P_4 %>%
  summarise(
    Mean = mean(Count_source, na.rm = TRUE),
    Median = median(Count_source, na.rm = TRUE),
    Q1 = quantile(Count_source, 0.25, na.rm = TRUE),  # First quartile (25th percentile)
    Q3 = quantile(Count_source, 0.75, na.rm = TRUE),  # Third quartile (75th percentile)
    P10 = quantile(Count_source, 0.10, na.rm = TRUE), # 10th percentile
    P90 = quantile(Count_source, 0.90, na.rm = TRUE), # 90th percentile
    Min = min(Count_source, na.rm = TRUE),
    Max = max(Count_source, na.rm = TRUE),
    Std_Dev = sd(Count_source, na.rm = TRUE),
    Variance = var(Count_source, na.rm = TRUE)
  )

# Display results
print(stats_Count_source)

####
####
####

###### PART C - df1157_miRNA

library(dplyr)

# Ensure "Count_source" is numeric before computing statistics
df1157_miRNA <- df1157_miRNA %>%
  mutate(Count_source = as.numeric(Count_source))

# Ensure "Risky_Count" is numeric before computing statistics
df1157_miRNA <- df1157_miRNA %>%
  mutate(Risky_Count = as.numeric(Risky_Count))

# Ensure "Count_source" is numeric before computing statistics
df1157_miRNA <- df1157_miRNA %>%
  mutate(Count_source = as.numeric(Count_source))

# Compute full statistics for "Count_source"
stats_Count_source <- df1157_miRNA %>%
  summarise(
    Mean = mean(Count_source, na.rm = TRUE),
    Median = median(Count_source, na.rm = TRUE),
    Q1 = quantile(Count_source, 0.25, na.rm = TRUE),  # First quartile (25th percentile)
    Q3 = quantile(Count_source, 0.75, na.rm = TRUE),  # Third quartile (75th percentile)
    P10 = quantile(Count_source, 0.10, na.rm = TRUE), # 10th percentile
    P90 = quantile(Count_source, 0.90, na.rm = TRUE), # 90th percentile
    Min = min(Count_source, na.rm = TRUE),
    Max = max(Count_source, na.rm = TRUE),
    Std_Dev = sd(Count_source, na.rm = TRUE),
    Variance = var(Count_source, na.rm = TRUE)
  )

# Display results
print(stats_Count_source)

# Define total row count of the original dataset dynamically
total_rows <- nrow(df1157_miRNA)

# Step 1: Further filter rows that contain "Count_source" == 1

df1157_miRNA_Count_1 <- df1157_miRNA %>%
  filter(Count_source == 1)  

# Store row count
n_miRNA_Count_1 <- nrow(df1157_miRNA_Count_1)  

# Compute percentage relative to df1157_miRNA
percent_n_miRNA_Count_1<- round((n_miRNA_Count_1/ nrow(df1157_miRNA)) * 100, 2)

# Display summary
print(paste("Rows in df1157_miRNA_Count_1:", n_miRNA_Count_1, "(", percent_n_miRNA_Count_1, "% of df1157_miRNA)"))

# Step 1: Further filter rows that contain "Risky" or "Protective" in any column
df1157_miRNA_SMC_R_P<- df1157_miRNA %>%
  filter(if_any(everything(), ~ grepl("Risky|Protective", ., ignore.case = TRUE)))  # Case-insensitive search for both terms

# Store row count
n_SMC_R_P <- nrow(df1157_miRNA_SMC_R_P)  

# Compute percentage relative to df1157_miRNA
percent_SMC_R_P <- round((n_SMC_R_P / nrow(df1157_miRNA)) * 100, 2)

# Display summary
print(paste("Rows in df1157_miRNA_SMC_R_P:", n_SMC_R_P, "(", percent_SMC_R_P, "% of df1157_miRNA)"))


# Step 2: Further filter rows that contain "Count_source" == 1

df1157_miRNA_SMC_R_P_1 <- df1157_miRNA_SMC_R_P %>%
  filter(Count_source == 1)  

# Store row count
n_SMC_R_P_1 <- nrow(df1157_miRNA_SMC_R_P_1)  

# Compute percentage relative to df1157_miRNA
percent_SMC_R_P_1<- round((n_SMC_R_P_1 / nrow(df1157_miRNA_SMC_R_P)) * 100, 2)

# Display summary
print(paste("Rows in df1157_miRNA_SMC_R_P_1:", n_SMC_R_P_1, "(", percent_SMC_R_P_1, "% of df1157_miRNA_SMC_R_P)"))

saveRDS(df1157, "Dataset_S1_full.rds")

rio::export(df1157, "Dataset_S1_full.tsv")

####
####
####
# Step 3: Further filter rows that contain "Risky" in all four metrics
df1157_miRNA_SMC_R_4<- df1157_miRNA_SMC_R_P %>%
  filter(Risky_Count == 4)  

# Store row count
n_SMC_R_4 <- nrow(df1157_miRNA_SMC_R_4)  

# Compute percentage relative to df1157_miRNA_SMC_R_4
percent_SMC_R_4 <- round((n_SMC_R_4 / nrow(df1157_miRNA_SMC_R_P)) * 100, 2)

# Display summary
print(paste("Rows in df1157_miRNA_SMC_R_4:", n_SMC_R_4, "(", percent_SMC_R_4, "% of df1157_miRNA_SMC_R_P)"))

#####
#####
#####
# Step 4: Further filter rows that contain "Protective" in all four metrics
df1157_miRNA_SMC_P_4<- df1157_miRNA_SMC_R_P %>%
  filter(Protective_Count == 4)  

# Store row count
n_SMC_P_4 <- nrow(df1157_miRNA_SMC_P_4)  

# Compute percentage relative to df1157_miRNA_SMC_P_4
percent_SMC_P_4 <- round((n_SMC_P_4 / nrow(df1157_miRNA_SMC_R_P)) * 100, 2)

# Display summary
print(paste("Rows in df1157_miRNA_SMC_P_4:", n_SMC_P_4, "(", percent_SMC_P_4, "% of df1157_miRNA_SMC_R_P)"))

##### 
##### 
##### 
###### PART D - df1157_Methylation
##### 
#####
#####

# Ensure "Count_source" is numeric before computing statistics
df1157_Methylation <- df1157_Methylation %>%
  mutate(Count_source = as.numeric(Count_source))

# Ensure "Risky_Count" is numeric before computing statistics
df1157_Methylation <- df1157_Methylation %>%
  mutate(Risky_Count = as.numeric(Risky_Count))

# Ensure "Count_source" is numeric before computing statistics
df1157_Methylation <- df1157_Methylation %>%
  mutate(Count_source = as.numeric(Count_source))

# Step 1: Summary stats
# Compute full statistics for "Count_source"
stats_Count_source <- df1157_Methylation %>%
  summarise(
    Mean = mean(Count_source, na.rm = TRUE),
    Median = median(Count_source, na.rm = TRUE),
    Q1 = quantile(Count_source, 0.25, na.rm = TRUE),  # First quartile (25th percentile)
    Q3 = quantile(Count_source, 0.75, na.rm = TRUE),  # Third quartile (75th percentile)
    P10 = quantile(Count_source, 0.10, na.rm = TRUE), # 10th percentile
    P90 = quantile(Count_source, 0.90, na.rm = TRUE), # 90th percentile
    Min = min(Count_source, na.rm = TRUE),
    Max = max(Count_source, na.rm = TRUE),
    Std_Dev = sd(Count_source, na.rm = TRUE),
    Variance = var(Count_source, na.rm = TRUE)
  )

# Display results
print(stats_Count_source)

# Define total row count of the original dataset dynamically
total_rows <- nrow(df1157_Methylation)

# Step 2: Further filter rows that contain "Count_source" == 1

df1157_Methylation_Count_1 <- df1157_Methylation %>%
  filter(Count_source == 1)  

# Store row count
n_Methylation_Count_1 <- nrow(df1157_Methylation_Count_1)  

# Compute percentage relative to df1157_Methylation
percent_n_Methylation_Count_1<- round((n_Methylation_Count_1/ nrow(df1157_Methylation)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Methylation_Count_1:", n_Methylation_Count_1, "(", percent_n_Methylation_Count_1, "% of df1157_Methylation)"))

# Step 3: Further filter rows that contain "Phenotype" == "TSM"

df1157_Methylation_TSM <- df1157_Methylation %>%
  filter(Phenotype == "TSM")  

# Store row count
n_Methylation_TSM <- nrow(df1157_Methylation_TSM)  

# Compute percentage relative to df1157_Methylation
percent_n_Methylation_TSM<- round((n_Methylation_TSM / nrow(df1157_Methylation)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Methylation_TSM:", n_Methylation_TSM, "(", percent_n_Methylation_TSM, "% of df1157_Methylation_TSM)"))

# Step 4: Further filter rows that contain "Risky_count" == 4

df1157_Methylation_TSM_R_4 <- df1157_Methylation_TSM %>%
  filter(Risky_Count == 4)  

# Store row count
n_Methylation_TSM_R_4 <- nrow(df1157_Methylation_TSM_R_4)  

# Compute percentage relative to df1157_Methylation
percent_n_Methylation_TSM_R_4<- round((n_Methylation_TSM_R_4 / nrow(df1157_Methylation_TSM)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Methylation_TSM_R_4:", n_Methylation_TSM_R_4, "(", percent_n_Methylation_TSM_R_4, "% of df1157_Methylation_TSM)"))

# Step 5: Further filter rows that contain "Protective_count" == 4

df1157_Methylation_TSM_P_4 <- df1157_Methylation_TSM %>%
  filter(Protective_Count == 4)  

# Store row count
n_Methylation_TSM_P_4 <- nrow(df1157_Methylation_TSM_P_4)  

# Compute percentage relative to df1157_Methylation
percent_n_Methylation_TSM_P_4<- round((n_Methylation_TSM_P_4 / nrow(df1157_Methylation_TSM)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Methylation_TSM_P_4:", n_Methylation_TSM_P_4, "(", percent_n_Methylation_TSM_P_4, "% of df1157_Methylation_TSM)"))

# Step 6: Further filter rows that contain "microenvironment_classification" == "anti-tumoral"

df1157_Methylation_TSM_P_4_anti <- df1157_Methylation_TSM_P_4 %>%
  filter(microenvironment_classification == "anti-tumoral")  

# Store row count
n_Methylation_TSM_P_4_anti <- nrow(df1157_Methylation_TSM_P_4_anti)  

# Compute percentage relative to df1157_Methylation
percent_n_Methylation_TSM_P_4_anti<- round((n_Methylation_TSM_P_4_anti / nrow(df1157_Methylation_TSM_P_4)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Methylation_TSM_P_4_anti:", n_Methylation_TSM_P_4_anti, "(", percent_n_Methylation_TSM_P_4_anti, "% of df1157_Methylation_TSM_P_4_anti)"))

# Step 7: Further filter rows that contain "Phenotype" == "TSM" and Risky_Count == 4 or Protective_Count == 4

df1157_Methylation_TSM_R_P_4 <- df1157_Methylation %>%
  filter(Phenotype == "TSM" & (Risky_Count == 4 | Protective_Count == 4))

# Store row count
n_Methylation_TSM_R_P_4 <- nrow(df1157_Methylation_TSM_R_P_4)  

# Compute percentage relative to df1157_Methylation
percent_n_Methylation_TSM_R_P_4<- round((n_Methylation_TSM_R_P_4 / nrow(df1157_Methylation_TSM)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Methylation_TSM_R_P_4:", n_Methylation_TSM_R_P_4, "(", percent_n_Methylation_TSM_R_P_4, "% of df1157_Methylation_TSM_R_P_4)"))

##tep 8: Further filter rows that contain "microenvironment_classification" == "anti-tumoral"

df1157_Methylation_TSM_R_P_4_anti <- df1157_Methylation_TSM_R_P_4 %>%
  filter(microenvironment_classification == "anti-tumoral")

# Store row count
n_Methylation_TSM_R_P_4_anti <- nrow(df1157_Methylation_TSM_R_P_4_anti)  

# Compute percentage relative to df1157_Methylation
percent_n_Methylation_TSM_R_P_4_anti<- round((n_Methylation_TSM_R_P_4_anti / nrow(df1157_Methylation_TSM_R_P_4)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Methylation_TSM_R_P_4_anti:", n_Methylation_TSM_R_P_4_anti, "(", percent_n_Methylation_TSM_R_P_4_anti, "% of df1157_Methylation_TSM_R_P_4)"))

## Step 9: Further filter rows that contain "microenvironment_classification" == "pro-tumoral"

df1157_Methylation_TSM_R_P_4_pro <- df1157_Methylation_TSM_R_P_4 %>%
  filter(microenvironment_classification == "pro-tumoral")

# Store row count
n_Methylation_TSM_R_P_4_pro <- nrow(df1157_Methylation_TSM_R_P_4_pro)  

# Compute percentage relative to df1157_Methylation
percent_n_Methylation_TSM_R_P_4_pro<- round((n_Methylation_TSM_R_P_4_pro / nrow(df1157_Methylation_TSM_R_P_4)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Methylation_TSM_R_P_4_pro:", n_Methylation_TSM_R_P_4_pro, "(", percent_n_Methylation_TSM_R_P_4_pro, "% of df1157_Methylation_TSM_R_P_4)"))

## Step 10: Further filter rows that contain "microenvironment_classification" == "dual"

df1157_Methylation_TSM_R_P_4_dual <- df1157_Methylation_TSM_R_P_4 %>%
  filter(microenvironment_classification == "dual")

# Store row count
n_Methylation_TSM_R_P_4_dual <- nrow(df1157_Methylation_TSM_R_P_4_dual)  

# Compute percentage relative to df1157_Methylation
percent_n_Methylation_TSM_R_P_4_dual<- round((n_Methylation_TSM_R_P_4_dual / nrow(df1157_Methylation_TSM_R_P_4)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Methylation_TSM_R_P_4_dual:", n_Methylation_TSM_R_P_4_dual, "(", percent_n_Methylation_TSM_R_P_4_dual, "% of df1157_Methylation_TSM_R_P_4)"))

## Step 11: Further filter rows that contain "immune_classification" == "Hot"

df1157_Methylation_TSM_R_P_4_Hot <- df1157_Methylation_TSM_R_P_4 %>%
  filter(immune_classification == "Hot")

# Store row count
n_Methylation_TSM_R_P_4_Hot <- nrow(df1157_Methylation_TSM_R_P_4_Hot)  

# Compute percentage relative to df1157_Methylation
percent_n_Methylation_TSM_R_P_4_Hot<- round((n_Methylation_TSM_R_P_4_Hot / nrow(df1157_Methylation_TSM_R_P_4)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Methylation_TSM_R_P_4_Hot:", n_Methylation_TSM_R_P_4_Hot, "(", percent_n_Methylation_TSM_R_P_4_Hot, "% of df1157_Methylation_TSM_R_P_4)"))

## Step 12: Further filter rows that contain "immune_classification" == "Cold"

df1157_Methylation_TSM_R_P_4_Cold <- df1157_Methylation_TSM_R_P_4 %>%
  filter(immune_classification == "Cold")

# Store row count
n_Methylation_TSM_R_P_4_Cold <- nrow(df1157_Methylation_TSM_R_P_4_Cold)  

# Compute percentage relative to df1157_Methylation
percent_n_Methylation_TSM_R_P_4_Cold<- round((n_Methylation_TSM_R_P_4_Cold / nrow(df1157_Methylation_TSM_R_P_4)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Methylation_TSM_R_P_4_Cold:", n_Methylation_TSM_R_P_4_Cold, "(", percent_n_Methylation_TSM_R_P_4_Cold, "% of df1157_Methylation_TSM_R_P_4)"))

## Step 13: Further filter rows that contain "immune_classification" == "Variable"

df1157_Methylation_TSM_R_P_4_Variable <- df1157_Methylation_TSM_R_P_4 %>%
  filter(immune_classification == "Variable")

# Store row count
n_Methylation_TSM_R_P_4_Variable <- nrow(df1157_Methylation_TSM_R_P_4_Variable)  

# Compute percentage relative to df1157_Methylation
percent_n_Methylation_TSM_R_P_4_Variable<- round((n_Methylation_TSM_R_P_4_Variable / nrow(df1157_Methylation_TSM_R_P_4)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Methylation_TSM_R_P_4_Variable:", n_Methylation_TSM_R_P_4_Variable, "(", percent_n_Methylation_TSM_R_P_4_Variable, "% of df1157_Methylation_TSM_R_P_4)"))

##### 
##### 
##### 
###### PART E - df1157_Protein
##### 
#####
#####

# Ensure "Count_source" is numeric before computing statistics
df1157_Protein <- df1157_Protein %>%
  mutate(Count_source = as.numeric(Count_source))

# Ensure "Risky_Count" is numeric before computing statistics
df1157_Protein <- df1157_Protein %>%
  mutate(Risky_Count = as.numeric(Risky_Count))

# Ensure "Count_source" is numeric before computing statistics
df1157_Protein <- df1157_Protein %>%
  mutate(Count_source = as.numeric(Count_source))

# Step 1: Summary stats
# Compute full statistics for "Count_source"
stats_Count_source <- df1157_Protein %>%
  summarise(
    Mean = mean(Count_source, na.rm = TRUE),
    Median = median(Count_source, na.rm = TRUE),
    Q1 = quantile(Count_source, 0.25, na.rm = TRUE),  # First quartile (25th percentile)
    Q3 = quantile(Count_source, 0.75, na.rm = TRUE),  # Third quartile (75th percentile)
    P10 = quantile(Count_source, 0.10, na.rm = TRUE), # 10th percentile
    P90 = quantile(Count_source, 0.90, na.rm = TRUE), # 90th percentile
    Min = min(Count_source, na.rm = TRUE),
    Max = max(Count_source, na.rm = TRUE),
    Std_Dev = sd(Count_source, na.rm = TRUE),
    Variance = var(Count_source, na.rm = TRUE)
  )

# Display results
print(stats_Count_source)

# Define total row count of the original dataset dynamically
total_rows <- nrow(df1157_Protein)

# Step 2: Further filter rows that contain "Count_source" == 1

df1157_Protein_Count_1 <- df1157_Protein %>%
  filter(Count_source == 1)  

# Store row count
n_Protein_Count_1 <- nrow(df1157_Protein_Count_1)  

# Compute percentage relative to df1157_Protein
percent_n_Protein_Count_1<- round((n_Protein_Count_1/ nrow(df1157_Protein)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Protein_Count_1:", n_Protein_Count_1, "(", percent_n_Protein_Count_1, "% of df1157_Protein)"))

# Step 3: Further filter rows that contain "Phenotype" == "TSM"

df1157_Protein_TSM <- df1157_Protein %>%
  filter(Phenotype == "TSM")  

# Store row count
n_Protein_TSM <- nrow(df1157_Protein_TSM)  

# Compute percentage relative to df1157_Protein
percent_n_Protein_TSM<- round((n_Protein_TSM / nrow(df1157_Protein)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Protein_TSM:", n_Protein_TSM, "(", percent_n_Protein_TSM, "% of df1157_Protein_TSM)"))

# Step 4: Further by filter(SCS == "N")

df1157_Protein_TSM_N <- df1157_Protein_TSM %>%
  filter(SCS == "N")  

# Store row count
n_Protein_TSM_N <- nrow(df1157_Protein_TSM_N)  

# Compute percentage relative to df1157_Protein
percent_n_Protein_TSM_N<- round((n_Protein_TSM_N / nrow(df1157_Protein_TSM)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Protein_TSM_N:", n_Protein_TSM_N, "(", percent_n_Protein_TSM_N, "% of df1157_Protein_TSM_N)"))

# Step 5: Further by filter(SCS == "P")

df1157_Protein_TSM_P <- df1157_Protein_TSM %>%
  filter(SCS == "P")  

# Store row count
n_Protein_TSM_P <- nrow(df1157_Protein_TSM_P)  

# Compute percentage relative to df1157_Protein
percent_n_Protein_TSM_P <- round((n_Protein_TSM_P / nrow(df1157_Protein_TSM)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Protein_TSM_P:", n_Protein_TSM_P, "(", percent_n_Protein_TSM_P, "% of df1157_Protein_TSM_N)"))

# Step 6: Further filter rows that contain "Risky_count" == 4

df1157_Protein_TSM_R_4 <- df1157_Protein_TSM %>%
  filter(Risky_Count == 4)  

# Store row count
n_Protein_TSM_R_4 <- nrow(df1157_Protein_TSM_R_4)  

# Compute percentage relative to df1157_Protein
percent_n_Protein_TSM_R_4<- round((n_Protein_TSM_R_4 / nrow(df1157_Protein_TSM)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Protein_TSM_R_4:", n_Protein_TSM_R_4, "(", percent_n_Protein_TSM_R_4, "% of df1157_Protein_TSM)"))

# Step 7: Further filter rows that contain "Protective_count" == 4

df1157_Protein_TSM_P_4 <- df1157_Protein_TSM %>%
  filter(Protective_Count == 4)  

# Store row count
n_Protein_TSM_P_4 <- nrow(df1157_Protein_TSM_P_4)  

# Compute percentage relative to df1157_Protein
percent_n_Protein_TSM_P_4<- round((n_Protein_TSM_P_4 / nrow(df1157_Protein_TSM)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Protein_TSM_P_4:", n_Protein_TSM_P_4, "(", percent_n_Protein_TSM_P_4, "% of df1157_Protein_TSM)"))

# Step 8: Filter df1157_Protein_pro rows that contain "microenvironment_classification" == "pro-tumoral"

df1157_Protein_pro <- df1157_Protein %>%
  filter(microenvironment_classification == "pro-tumoral")  

# Store row count
n_Protein_pro <- nrow(df1157_Protein_pro)  

# Compute percentage relative to df1157_Protein
percent_n_Protein_pro<- round((n_Protein_pro / nrow(df1157_Protein)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Protein_pro:", n_Protein_pro, "(", percent_n_Protein_pro, "% of df1157_Protein_anti)"))

# Step 9: filter df1157_Protein_pro rows that contain "microenvironment_classification" == "anti-tumoral"

df1157_Protein_anti <- df1157_Protein %>%
  filter(microenvironment_classification == "anti-tumoral")  

# Store row count
n_Protein_anti <- nrow(df1157_Protein_anti)  

# Compute percentage relative to df1157_Protein
percent_n_Protein_anti<- round((n_Protein_anti / nrow(df1157_Protein)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Protein_anti:", n_Protein_anti, "(", percent_n_Protein_anti, "% of df1157_Protein_anti)"))

# Step 10: filter df1157_Protein_pro rows that contain "microenvironment_classification" == "dual"

df1157_Protein_dual <- df1157_Protein %>%
  filter(microenvironment_classification == "dual")  

# Store row count
n_Protein_dual <- nrow(df1157_Protein_dual)  

# Compute percentage relative to df1157_Protein
percent_n_Protein_dual<- round((n_Protein_dual / nrow(df1157_Protein)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Protein_dual:", n_Protein_dual, "(", percent_n_Protein_dual, "% of df1157_Protein_anti)"))

# Step 11: filter df1157_Protein_pro rows that contain "immune_classification" == "Hot"

df1157_Protein_Hot <- df1157_Protein %>%
  filter(immune_classification == "Hot")  

# Store row count
n_Protein_Hot <- nrow(df1157_Protein_Hot)  

# Compute percentage relative to df1157_Protein
percent_n_Protein_Hot<- round((n_Protein_Hot / nrow(df1157_Protein)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Protein_Hot:", n_Protein_Hot, "(", percent_n_Protein_Hot, "% of df1157_Protein_anti)"))

# Step 12: filter df1157_Protein_pro rows that contain "immune_classification" == "Cold"

df1157_Protein_Cold <- df1157_Protein %>%
  filter(immune_classification == "Cold")  

# Store row count
n_Protein_Cold <- nrow(df1157_Protein_Cold)  

# Compute percentage relative to df1157_Protein
percent_n_Protein_Cold<- round((n_Protein_Cold / nrow(df1157_Protein)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Protein_Cold:", n_Protein_Cold, "(", percent_n_Protein_Cold, "% of df1157_Protein_anti)"))

# Step 13: filter df1157_Protein_pro rows that contain "immune_classification" == "Variable"

df1157_Protein_Variable <- df1157_Protein %>%
  filter(immune_classification == "Variable")  

# Store row count
n_Protein_Variable <- nrow(df1157_Protein_Variable)  

# Compute percentage relative to df1157_Protein
percent_n_Protein_Variable<- round((n_Protein_Variable / nrow(df1157_Protein)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Protein_Variable:", n_Protein_Variable, "(", percent_n_Protein_Variable, "% of df1157_Protein_anti)"))

# Step 14: Further filter rows that contain "Phenotype" == "TSM" and Risky_Count == 4 or Protective_Count == 4

df1157_Protein_TSM_R_P_4 <- df1157_Protein %>%
  filter(Phenotype == "TSM" & (Risky_Count == 4 | Protective_Count == 4))

# Store row count
n_Protein_TSM_R_P_4 <- nrow(df1157_Protein_TSM_R_P_4)  

# Compute percentage relative to df1157_Protein
percent_n_Protein_TSM_R_P_4<- round((n_Protein_TSM_R_P_4 / nrow(df1157_Protein_TSM)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Protein_TSM_R_P_4:", n_Protein_TSM_R_P_4, "(", percent_n_Protein_TSM_R_P_4, "% of df1157_Protein_TSM_R_P_4)"))

##tep 8: Further filter rows that contain "microenvironment_classification" == "anti-tumoral"

df1157_Protein_TSM_R_P_4_anti <- df1157_Protein_TSM_R_P_4 %>%
  filter(microenvironment_classification == "anti-tumoral")

# Store row count
n_Protein_TSM_R_P_4_anti <- nrow(df1157_Protein_TSM_R_P_4_anti)  

# Compute percentage relative to df1157_Protein
percent_n_Protein_TSM_R_P_4_anti<- round((n_Protein_TSM_R_P_4_anti / nrow(df1157_Protein_TSM_R_P_4)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Protein_TSM_R_P_4_anti:", n_Protein_TSM_R_P_4_anti, "(", percent_n_Protein_TSM_R_P_4_anti, "% of df1157_Protein_TSM_R_P_4)"))

## Step 15: Further filter rows that contain "microenvironment_classification" == "pro-tumoral"

df1157_Protein_TSM_R_P_4_pro <- df1157_Protein_TSM_R_P_4 %>%
  filter(microenvironment_classification == "pro-tumoral")

# Store row count
n_Protein_TSM_R_P_4_pro <- nrow(df1157_Protein_TSM_R_P_4_pro)  

# Compute percentage relative to df1157_Protein
percent_n_Protein_TSM_R_P_4_pro<- round((n_Protein_TSM_R_P_4_pro / nrow(df1157_Protein_TSM_R_P_4)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Protein_TSM_R_P_4_pro:", n_Protein_TSM_R_P_4_pro, "(", percent_n_Protein_TSM_R_P_4_pro, "% of df1157_Protein_TSM_R_P_4)"))

## Step 16: Further filter rows that contain "microenvironment_classification" == "dual"

df1157_Protein_TSM_R_P_4_dual <- df1157_Protein_TSM_R_P_4 %>%
  filter(microenvironment_classification == "dual")

# Store row count
n_Protein_TSM_R_P_4_dual <- nrow(df1157_Protein_TSM_R_P_4_dual)  

# Compute percentage relative to df1157_Protein
percent_n_Protein_TSM_R_P_4_dual<- round((n_Protein_TSM_R_P_4_dual / nrow(df1157_Protein_TSM_R_P_4)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Protein_TSM_R_P_4_dual:", n_Protein_TSM_R_P_4_dual, "(", percent_n_Protein_TSM_R_P_4_dual, "% of df1157_Protein_TSM_R_P_4)"))

## Step 17: Further filter rows that contain "immune_classification" == "Hot"

df1157_Protein_TSM_R_P_4_Hot <- df1157_Protein_TSM_R_P_4 %>%
  filter(immune_classification == "Hot")

# Store row count
n_Protein_TSM_R_P_4_Hot <- nrow(df1157_Protein_TSM_R_P_4_Hot)  

# Compute percentage relative to df1157_Protein
percent_n_Protein_TSM_R_P_4_Hot<- round((n_Protein_TSM_R_P_4_Hot / nrow(df1157_Protein_TSM_R_P_4)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Protein_TSM_R_P_4_Hot:", n_Protein_TSM_R_P_4_Hot, "(", percent_n_Protein_TSM_R_P_4_Hot, "% of df1157_Protein_TSM_R_P_4)"))

## Step 18: Further filter rows that contain "immune_classification" == "Cold"

df1157_Protein_TSM_R_P_4_Cold <- df1157_Protein_TSM_R_P_4 %>%
  filter(immune_classification == "Cold")

# Store row count
n_Protein_TSM_R_P_4_Cold <- nrow(df1157_Protein_TSM_R_P_4_Cold)  

# Compute percentage relative to df1157_Protein
percent_n_Protein_TSM_R_P_4_Cold<- round((n_Protein_TSM_R_P_4_Cold / nrow(df1157_Protein_TSM_R_P_4)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Protein_TSM_R_P_4_Cold:", n_Protein_TSM_R_P_4_Cold, "(", percent_n_Protein_TSM_R_P_4_Cold, "% of df1157_Protein_TSM_R_P_4)"))

## Step 19: Further filter rows that contain "immune_classification" == "Variable"

df1157_Protein_TSM_R_P_4_Variable <- df1157_Protein_TSM_R_P_4 %>%
  filter(immune_classification == "Variable")

# Store row count
n_Protein_TSM_R_P_4_Variable <- nrow(df1157_Protein_TSM_R_P_4_Variable)  

# Compute percentage relative to df1157_Protein
percent_n_Protein_TSM_R_P_4_Variable<- round((n_Protein_TSM_R_P_4_Variable / nrow(df1157_Protein_TSM_R_P_4)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Protein_TSM_R_P_4_Variable:", n_Protein_TSM_R_P_4_Variable, "(", percent_n_Protein_TSM_R_P_4_Variable, "% of df1157_Protein_TSM_R_P_4)"))

##### 
##### 
##### 
###### PART F - df1157_Mutation
##### 
#####
#####

# Ensure "Count_source" is numeric before computing statistics
df1157_Mutation <- df1157_Mutation %>%
  mutate(Count_source = as.numeric(Count_source))

# Ensure "Risky_Count" is numeric before computing statistics
df1157_Mutation <- df1157_Mutation %>%
  mutate(Risky_Count = as.numeric(Risky_Count))

# Ensure "Count_source" is numeric before computing statistics
df1157_Mutation <- df1157_Mutation %>%
  mutate(Count_source = as.numeric(Count_source))

# Step 1: Summary stats
# Compute full statistics for "Count_source"
stats_Count_source <- df1157_Mutation %>%
  summarise(
    Mean = mean(Count_source, na.rm = TRUE),
    Median = median(Count_source, na.rm = TRUE),
    Q1 = quantile(Count_source, 0.25, na.rm = TRUE),  # First quartile (25th percentile)
    Q3 = quantile(Count_source, 0.75, na.rm = TRUE),  # Third quartile (75th percentile)
    P10 = quantile(Count_source, 0.10, na.rm = TRUE), # 10th percentile
    P90 = quantile(Count_source, 0.90, na.rm = TRUE), # 90th percentile
    Min = min(Count_source, na.rm = TRUE),
    Max = max(Count_source, na.rm = TRUE),
    Std_Dev = sd(Count_source, na.rm = TRUE),
    Variance = var(Count_source, na.rm = TRUE)
  )

# Display results
print(stats_Count_source)

# Define total row count of the original dataset dynamically
total_rows <- nrow(df1157_Mutation)

# Step 2: Further filter rows that contain "Count_source" == 1

df1157_Mutation_Count_1 <- df1157_Mutation %>%
  filter(Count_source == 1)  

# Store row count
n_Mutation_Count_1 <- nrow(df1157_Mutation_Count_1)  

# Compute percentage relative to df1157_Mutation
percent_n_Mutation_Count_1<- round((n_Mutation_Count_1/ nrow(df1157_Mutation)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Mutation_Count_1:", n_Mutation_Count_1, "(", percent_n_Mutation_Count_1, "% of df1157_Mutation)"))

# Step 3: Further filter rows that contain "Phenotype" == "TMB"

df1157_Mutation_TMB <- df1157_Mutation %>%
  filter(Phenotype == "TMB")  


df1157_Mutation_TMB_P <- df1157_Mutation_TMB %>%
  filter(SCS == "P") 

# Store row count
n_Mutation_TMB <- nrow(df1157_Mutation_TMB)  

# Compute percentage relative to df1157_Mutation
percent_n_Mutation_TMB<- round((n_Mutation_TMB / nrow(df1157_Mutation)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Mutation_TMB:", n_Mutation_TMB, "(", percent_n_Mutation_TMB, "% of df1157_Mutation_TMB)"))

# Step 4: Further filter rows that contain "Phenotype" == "TSM"

df1157_Mutation_TSM <- df1157_Mutation %>%
  filter(Phenotype == "TSM")  

# Store row count
n_Mutation_TSM <- nrow(df1157_Mutation_TSM)  

# Compute percentage relative to df1157_Mutation
percent_n_Mutation_TSM<- round((n_Mutation_TSM / nrow(df1157_Mutation)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Mutation_TSM:", n_Mutation_TSM, "(", percent_n_Mutation_TSM, "% of df1157_Mutation_TSM)"))

df1157_Mutation_TSM_P <- df1157_Mutation_TSM %>%
  filter(SCS == "P") 

# Store row count
n_Mutation_TSM_P <- nrow(df1157_Mutation_TSM_P)  

# Compute percentage relative to df1157_Mutation
percent_n_Mutation_TSM_P <- round((n_Mutation_TSM_P / nrow(df1157_Mutation)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Mutation_TSM_P:", n_Mutation_TSM_P, "(", percent_n_Mutation_TSM_P, "% of df1157_Mutation_TSM_P)"))

# Step 5: Further filter rows that contain "Phenotype" == "MSI"

df1157_Mutation_MSI <- df1157_Mutation %>%
  filter(Phenotype == "MSI")  

df1157_Mutation_MSI_P <- df1157_Mutation_MSI %>%
  filter(SCS == "P")  

# Store row count
n_Mutation_MSI <- nrow(df1157_Mutation_MSI)  

# Compute percentage relative to df1157_Mutation
percent_n_Mutation_MSI<- round((n_Mutation_MSI / nrow(df1157_Mutation)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Mutation_MSI:", n_Mutation_MSI, "(", percent_n_Mutation_MSI, "% of df1157_Mutation_MSI)"))

# Step 6: Further by filter(SCS == "N")

df1157_Mutation_TMB_N <- df1157_Mutation_TMB %>%
  filter(SCS == "N")  

# Store row count
n_Mutation_TMB_N <- nrow(df1157_Mutation_TMB_N)  

# Compute percentage relative to df1157_Mutation
percent_n_Mutation_TMB_N<- round((n_Mutation_TMB_N / nrow(df1157_Mutation_TMB)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Mutation_TMB_N:", n_Mutation_TMB_N, "(", percent_n_Mutation_TMB_N, "% of df1157_Mutation_TMB_N)"))

# Step 6a: Further by filter(SCS == "P")

df1157_Mutation_TMB_P <- df1157_Mutation_TMB %>%
  filter(SCS == "P")  

# Store row count
n_Mutation_TMB_P <- nrow(df1157_Mutation_TMB_P)  

# Compute percentage relative to df1157_Mutation
percent_n_Mutation_TMB_P <- round((n_Mutation_TMB_P / nrow(df1157_Mutation_TMB)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Mutation_TMB_P:", n_Mutation_TMB_P, "(", percent_n_Mutation_TMB_P, "% of df1157_Mutation_TMB_N)"))

# Step 7: Further filter rows that contain "Risky_count" == 4

df1157_Mutation_TMB_R_4 <- df1157_Mutation_TMB %>%
  filter(Risky_Count == 4)  

# Store row count
n_Mutation_TMB_R_4 <- nrow(df1157_Mutation_TMB_R_4)  

# Compute percentage relative to df1157_Mutation
percent_n_Mutation_TMB_R_4<- round((n_Mutation_TMB_R_4 / nrow(df1157_Mutation_TMB)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Mutation_TMB_R_4:", n_Mutation_TMB_R_4, "(", percent_n_Mutation_TMB_R_4, "% of df1157_Mutation_TMB)"))

# Step 7a: Further filter rows that contain "Protective_count" == 4

df1157_Mutation_TMB_P_4 <- df1157_Mutation_TMB %>%
  filter(Protective_Count == 4)  

# Store row count
n_Mutation_TMB_P_4 <- nrow(df1157_Mutation_TMB_P_4)  

# Compute percentage relative to df1157_Mutation
percent_n_Mutation_TMB_P_4<- round((n_Mutation_TMB_P_4 / nrow(df1157_Mutation_TMB)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Mutation_TMB_P_4:", n_Mutation_TMB_P_4, "(", percent_n_Mutation_TMB_P_4, "% of df1157_Mutation_TMB)"))

# Step 8: Filter df1157_Mutation_TMB rows that contain "microenvironment_classification" == "pro-tumoral"

df1157_Mutation_TMB_pro <- df1157_Mutation_TMB %>%
  filter(microenvironment_classification == "pro-tumoral")  

# Store row count
n_Mutation_TMB_pro <- nrow(df1157_Mutation_TMB_pro)  

# Compute percentage relative to df1157_Mutation
percent_n_Mutation_TMB_pro<- round((n_Mutation_TMB_pro / nrow(df1157_Mutation_TMB)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Mutation_TMB_pro:", n_Mutation_TMB_pro, "(", percent_n_Mutation_TMB_pro, "% of df1157_Mutation_TMB_pro)"))

# Step 9: Filter df1157_Mutation_TMB rows that contain "microenvironment_classification" == "anti-tumoral"

df1157_Mutation_TMB_anti <- df1157_Mutation_TMB %>%
  filter(microenvironment_classification == "anti-tumoral")  

# Store row count
n_Mutation_TMB_anti <- nrow(df1157_Mutation_TMB_anti)  

# Compute percentage relative to df1157_Mutation
percent_n_Mutation_TMB_anti<- round((n_Mutation_TMB_anti / nrow(df1157_Mutation_TMB)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Mutation_TMB_anti:", n_Mutation_TMB_anti, "(", percent_n_Mutation_TMB_anti, "% of df1157_Mutation_TMB_anti)"))

# Step 10: Filter df1157_Mutation_TMB rows that contain "microenvironment_classification" == "dual"

df1157_Mutation_TMB_dual <- df1157_Mutation_TMB %>%
  filter(microenvironment_classification == "dual")  

# Store row count
n_Mutation_TMB_dual <- nrow(df1157_Mutation_TMB_dual)  

# Compute percentage relative to df1157_Mutation
percent_n_Mutation_TMB_dual<- round((n_Mutation_TMB_dual / nrow(df1157_Mutation_TMB)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Mutation_TMB_dual:", n_Mutation_TMB_dual, "(", percent_n_Mutation_TMB_dual, "% of df1157_Mutation_TMB)"))

# Step 11: filter df1157_Mutation_TMB rows that contain "immune_classification" == "Hot"

df1157_Mutation_TMB_Hot <- df1157_Mutation_TMB %>%
  filter(immune_classification == "Hot")  

# Store row count
n_Mutation_TMB_Hot <- nrow(df1157_Mutation_TMB_Hot)  

# Compute percentage relative to df1157_Mutation
percent_n_Mutation_TMB_Hot<- round((n_Mutation_TMB_Hot / nrow(df1157_Mutation_TMB)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Mutation_TMB_Hot:", n_Mutation_TMB_Hot, "(", percent_n_Mutation_TMB_Hot, "% of df1157_Mutation_TMB)"))

# Step 12: filter df1157_Mutation_TMB rows that contain "immune_classification" == "Variable"

df1157_Mutation_TMB_Variable <- df1157_Mutation_TMB %>%
  filter(immune_classification == "Variable")  

# Store row count
n_Mutation_TMB_Variable <- nrow(df1157_Mutation_TMB_Variable)  

# Compute percentage relative to df1157_Mutation
percent_n_Mutation_TMB_Variable<- round((n_Mutation_TMB_Variable / nrow(df1157_Mutation_TMB)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Mutation_TMB_Variable:", n_Mutation_TMB_Variable, "(", percent_n_Mutation_TMB_Variable, "% of df1157_Mutation_TMB)"))

# Step 13: filter df1157_Mutation_pro rows that contain "immune_classification" == "Cold"

df1157_Mutation_Cold <- df1157_Mutation %>%
  filter(immune_classification == "Cold")  

# Store row count
n_Mutation_Cold <- nrow(df1157_Mutation_Cold)  

# Compute percentage relative to df1157_Mutation
percent_n_Mutation_Cold<- round((n_Mutation_Cold / nrow(df1157_Mutation)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Mutation_Cold:", n_Mutation_Cold, "(", percent_n_Mutation_Cold, "% of df1157_Mutation_anti)"))

# Step 14: filter df1157_Mutation_pro rows that contain "immune_classification" == "Variable"

df1157_Mutation_Variable <- df1157_Mutation %>%
  filter(immune_classification == "Variable")  

# Store row count
n_Mutation_Variable <- nrow(df1157_Mutation_Variable)  

# Compute percentage relative to df1157_Mutation
percent_n_Mutation_Variable<- round((n_Mutation_Variable / nrow(df1157_Mutation)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Mutation_Variable:", n_Mutation_Variable, "(", percent_n_Mutation_Variable, "% of df1157_Mutation_anti)"))

# Step 15: filter df1157_Mutation_pro rows that contain "microenvironment_classification" == "anti-tumoral"

df1157_Mutation_anti <- df1157_Mutation %>%
  filter(microenvironment_classification == "anti-tumoral")  

# Store row count
n_Mutation_anti <- nrow(df1157_Mutation_anti)  

# Compute percentage relative to df1157_Mutation
percent_n_Mutation_anti<- round((n_Mutation_anti / nrow(df1157_Mutation)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Mutation_anti:", n_Mutation_anti, "(", percent_n_Mutation_anti, "% of df1157_Mutation_anti)"))

# Step 16: filter df1157_Mutation_pro rows that contain "microenvironment_classification" == "dual"

df1157_Mutation_dual <- df1157_Mutation %>%
  filter(microenvironment_classification == "dual")  

# Store row count
n_Mutation_dual <- nrow(df1157_Mutation_dual)  

# Compute percentage relative to df1157_Mutation
percent_n_Mutation_dual<- round((n_Mutation_dual / nrow(df1157_Mutation)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Mutation_dual:", n_Mutation_dual, "(", percent_n_Mutation_dual, "% of df1157_Mutation_anti)"))

# Step 17: filter df1157_Mutation_pro rows that contain "immune_classification" == "Hot"

df1157_Mutation_Hot <- df1157_Mutation %>%
  filter(immune_classification == "Hot")  

# Store row count
n_Mutation_Hot <- nrow(df1157_Mutation_Hot)  

# Compute percentage relative to df1157_Mutation
percent_n_Mutation_Hot<- round((n_Mutation_Hot / nrow(df1157_Mutation)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Mutation_Hot:", n_Mutation_Hot, "(", percent_n_Mutation_Hot, "% of df1157_Mutation_anti)"))

# Step 18: filter df1157_Mutation_pro rows that contain "immune_classification" == "Cold"

df1157_Mutation_Cold <- df1157_Mutation %>%
  filter(immune_classification == "Cold")  

# Store row count
n_Mutation_Cold <- nrow(df1157_Mutation_Cold)  

# Compute percentage relative to df1157_Mutation
percent_n_Mutation_Cold<- round((n_Mutation_Cold / nrow(df1157_Mutation)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Mutation_Cold:", n_Mutation_Cold, "(", percent_n_Mutation_Cold, "% of df1157_Mutation_anti)"))

# Step 19: filter df1157_Mutation_pro rows that contain "immune_classification" == "Variable"

df1157_Mutation_Variable <- df1157_Mutation %>%
  filter(immune_classification == "Variable")  

# Store row count
n_Mutation_Variable <- nrow(df1157_Mutation_Variable)  

# Compute percentage relative to df1157_Mutation
percent_n_Mutation_Variable<- round((n_Mutation_Variable / nrow(df1157_Mutation)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Mutation_Variable:", n_Mutation_Variable, "(", percent_n_Mutation_Variable, "% of df1157_Mutation_anti)"))

# Step 20: Further filter rows that contain "Phenotype" == "TMB" and Risky_Count == 4 or Protective_Count == 4

df1157_Mutation_TMB_R_P_4 <- df1157_Mutation %>%
  filter(Phenotype == "TMB" & (Risky_Count == 4 | Protective_Count == 4))

# Store row count
n_Mutation_TMB_R_P_4 <- nrow(df1157_Mutation_TMB_R_P_4)  

# Compute percentage relative to df1157_Mutation
percent_n_Mutation_TMB_R_P_4<- round((n_Mutation_TMB_R_P_4 / nrow(df1157_Mutation_TMB)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Mutation_TMB_R_P_4:", n_Mutation_TMB_R_P_4, "(", percent_n_Mutation_TMB_R_P_4, "% of df1157_Mutation_TMB_R_P_4)"))

##tep 8: Further filter rows that contain "microenvironment_classification" == "anti-tumoral"

df1157_Mutation_TMB_R_P_4_anti <- df1157_Mutation_TMB_R_P_4 %>%
  filter(microenvironment_classification == "anti-tumoral")

# Store row count
n_Mutation_TMB_R_P_4_anti <- nrow(df1157_Mutation_TMB_R_P_4_anti)  

# Compute percentage relative to df1157_Mutation
percent_n_Mutation_TMB_R_P_4_anti<- round((n_Mutation_TMB_R_P_4_anti / nrow(df1157_Mutation_TMB_R_P_4)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Mutation_TMB_R_P_4_anti:", n_Mutation_TMB_R_P_4_anti, "(", percent_n_Mutation_TMB_R_P_4_anti, "% of df1157_Mutation_TMB_R_P_4)"))

## Step 21: Further filter rows that contain "microenvironment_classification" == "pro-tumoral"

df1157_Mutation_TMB_R_P_4_pro <- df1157_Mutation_TMB_R_P_4 %>%
  filter(microenvironment_classification == "pro-tumoral")

# Store row count
n_Mutation_TMB_R_P_4_pro <- nrow(df1157_Mutation_TMB_R_P_4_pro)  

# Compute percentage relative to df1157_Mutation
percent_n_Mutation_TMB_R_P_4_pro<- round((n_Mutation_TMB_R_P_4_pro / nrow(df1157_Mutation_TMB_R_P_4)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Mutation_TMB_R_P_4_pro:", n_Mutation_TMB_R_P_4_pro, "(", percent_n_Mutation_TMB_R_P_4_pro, "% of df1157_Mutation_TMB_R_P_4)"))

## Step 22: Further filter rows that contain "microenvironment_classification" == "dual"

df1157_Mutation_TMB_R_P_4_dual <- df1157_Mutation_TMB_R_P_4 %>%
  filter(microenvironment_classification == "dual")

# Store row count
n_Mutation_TMB_R_P_4_dual <- nrow(df1157_Mutation_TMB_R_P_4_dual)  

# Compute percentage relative to df1157_Mutation
percent_n_Mutation_TMB_R_P_4_dual<- round((n_Mutation_TMB_R_P_4_dual / nrow(df1157_Mutation_TMB_R_P_4)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Mutation_TMB_R_P_4_dual:", n_Mutation_TMB_R_P_4_dual, "(", percent_n_Mutation_TMB_R_P_4_dual, "% of df1157_Mutation_TMB_R_P_4)"))

## Step 23: Further filter rows that contain "immune_classification" == "Hot"

df1157_Mutation_TMB_R_P_4_Hot <- df1157_Mutation_TMB_R_P_4 %>%
  filter(immune_classification == "Hot")

# Store row count
n_Mutation_TMB_R_P_4_Hot <- nrow(df1157_Mutation_TMB_R_P_4_Hot)  

# Compute percentage relative to df1157_Mutation
percent_n_Mutation_TMB_R_P_4_Hot<- round((n_Mutation_TMB_R_P_4_Hot / nrow(df1157_Mutation_TMB_R_P_4)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Mutation_TMB_R_P_4_Hot:", n_Mutation_TMB_R_P_4_Hot, "(", percent_n_Mutation_TMB_R_P_4_Hot, "% of df1157_Mutation_TMB_R_P_4)"))

## Step 24: Further filter rows that contain "immune_classification" == "Cold"

df1157_Mutation_TMB_R_P_4_Cold <- df1157_Mutation_TMB_R_P_4 %>%
  filter(immune_classification == "Cold")

# Store row count
n_Mutation_TMB_R_P_4_Cold <- nrow(df1157_Mutation_TMB_R_P_4_Cold)  

# Compute percentage relative to df1157_Mutation
percent_n_Mutation_TMB_R_P_4_Cold<- round((n_Mutation_TMB_R_P_4_Cold / nrow(df1157_Mutation_TMB_R_P_4)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Mutation_TMB_R_P_4_Cold:", n_Mutation_TMB_R_P_4_Cold, "(", percent_n_Mutation_TMB_R_P_4_Cold, "% of df1157_Mutation_TMB_R_P_4)"))

## Step 25: Further filter rows that contain "immune_classification" == "Variable"

df1157_Mutation_TMB_R_P_4_Variable <- df1157_Mutation_TMB_R_P_4 %>%
  filter(immune_classification == "Variable")

# Store row count
n_Mutation_TMB_R_P_4_Variable <- nrow(df1157_Mutation_TMB_R_P_4_Variable)  

# Compute percentage relative to df1157_Mutation
percent_n_Mutation_TMB_R_P_4_Variable<- round((n_Mutation_TMB_R_P_4_Variable / nrow(df1157_Mutation_TMB_R_P_4)) * 100, 2)

# Display summary
print(paste("Rows in df1157_Mutation_TMB_R_P_4_Variable:", n_Mutation_TMB_R_P_4_Variable, "(", percent_n_Mutation_TMB_R_P_4_Variable, "% of df1157_Mutation_TMB_R_P_4)"))

##### 
##### 
##### 
###### PART G - df1157_CNV
##### 
#####
#####

# Ensure "Count_source" is numeric before computing statistics
df1157_CNV <- df1157_CNV %>%
  mutate(Count_source = as.numeric(Count_source))

# Ensure "Risky_Count" is numeric before computing statistics
df1157_CNV <- df1157_CNV %>%
  mutate(Risky_Count = as.numeric(Risky_Count))

# Ensure "Count_source" is numeric before computing statistics
df1157_CNV <- df1157_CNV %>%
  mutate(Count_source = as.numeric(Count_source))

# Step 1: Summary stats
# Compute full statistics for "Count_source"
stats_Count_source <- df1157_CNV %>%
  summarise(
    Mean = mean(Count_source, na.rm = TRUE),
    Median = median(Count_source, na.rm = TRUE),
    Q1 = quantile(Count_source, 0.25, na.rm = TRUE),  # First quartile (25th percentile)
    Q3 = quantile(Count_source, 0.75, na.rm = TRUE),  # Third quartile (75th percentile)
    P10 = quantile(Count_source, 0.10, na.rm = TRUE), # 10th percentile
    P90 = quantile(Count_source, 0.90, na.rm = TRUE), # 90th percentile
    Min = min(Count_source, na.rm = TRUE),
    Max = max(Count_source, na.rm = TRUE),
    Std_Dev = sd(Count_source, na.rm = TRUE),
    Variance = var(Count_source, na.rm = TRUE)
  )

# Display results
print(stats_Count_source)

# Define total row count of the original dataset dynamically
total_rows <- nrow(df1157_CNV)

# Step 2: Further filter rows that contain "Count_source" == 1

df1157_CNV_Count_1 <- df1157_CNV %>%
  filter(Count_source == 1)  

# Store row count
n_CNV_Count_1 <- nrow(df1157_CNV_Count_1)  

# Compute percentage relative to df1157_CNV
percent_n_CNV_Count_1<- round((n_CNV_Count_1/ nrow(df1157_CNV)) * 100, 2)

# Display summary
print(paste("Rows in df1157_CNV_Count_1:", n_CNV_Count_1, "(", percent_n_CNV_Count_1, "% of df1157_CNV)"))

# Step 2a: Further filter rows that contain "Count_source" != 1

df1157_CNV_Count_d1 <- df1157_CNV %>%
  filter(Count_source !=1)  

# Store row count
n_CNV_Count_d1 <- nrow(df1157_CNV_Count_d1)  

# Compute percentage relative to df1157_CNV
percent_n_CNV_Count_d1<- round((n_CNV_Count_d1/ nrow(df1157_CNV)) * 100, 2)

# Display summary
print(paste("Rows in df1157_CNV_Count_d1:", n_CNV_Count_d1, "(", percent_n_CNV_Count_d1, "% of df1157_CNV)"))

# Step 3: Further filter rows that contain "Phenotype" == "TMB"

df1157_CNV_TMB <- df1157_CNV %>%
  filter(Phenotype == "TMB")  


df1157_CNV_TMB_P <- df1157_CNV_TMB %>%
  filter(SCS == "P") 

# Store row count
n_CNV_TMB <- nrow(df1157_CNV_TMB)  

# Compute percentage relative to df1157_CNV
percent_n_CNV_TMB<- round((n_CNV_TMB / nrow(df1157_CNV)) * 100, 2)

# Display summary
print(paste("Rows in df1157_CNV_TMB:", n_CNV_TMB, "(", percent_n_CNV_TMB, "% of df1157_CNV_TMB)"))

# Step 4: Further filter rows that contain "Phenotype" == "TSM"

df1157_CNV_TSM <- df1157_CNV %>%
  filter(Phenotype == "TSM")  

# Store row count
n_CNV_TSM <- nrow(df1157_CNV_TSM)  

# Compute percentage relative to df1157_CNV
percent_n_CNV_TSM<- round((n_CNV_TSM / nrow(df1157_CNV)) * 100, 2)

# Display summary
print(paste("Rows in df1157_CNV_TSM:", n_CNV_TSM, "(", percent_n_CNV_TSM, "% of df1157_CNV_TSM)"))

df1157_CNV_TSM_N <- df1157_CNV_TSM %>%
  filter(SCS == "N") 

# Store row count
n_CNV_TSM_N <- nrow(df1157_CNV_TSM_N)  

# Compute percentage relative to df1157_CNV
percent_n_CNV_TSM_N <- round((n_CNV_TSM_N / nrow(df1157_CNV_TSM)) * 100, 2)

# Display summary
print(paste("Rows in df1157_CNV_TSM_N:", n_CNV_TSM_N, "(", percent_n_CNV_TSM_N, "% of df1157_CNV_TSM_N)"))

df1157_CNV_TSM_P <- df1157_CNV_TSM %>%
  filter(SCS == "P") 

# Store row count
n_CNV_TSM_P <- nrow(df1157_CNV_TSM_P)  

# Compute percentage relative to df1157_CNV
percent_n_CNV_TSM_P <- round((n_CNV_TSM_P / nrow(df1157_CNV_TSM)) * 100, 2)

# Display summary
print(paste("Rows in df1157_CNV_TSM_P:", n_CNV_TSM_P, "(", percent_n_CNV_TSM_P, "% of df1157_CNV_TSM_P)"))

# Step 5: Further filter rows that contain "Phenotype" == "MSI"

df1157_CNV_MSI <- df1157_CNV %>%
  filter(Phenotype == "MSI")  

df1157_CNV_MSI_P <- df1157_CNV_MSI %>%
  filter(SCS == "P")  

# Store row count
n_CNV_MSI <- nrow(df1157_CNV_MSI)  

# Compute percentage relative to df1157_CNV
percent_n_CNV_MSI<- round((n_CNV_MSI / nrow(df1157_CNV)) * 100, 2)

# Display summary
print(paste("Rows in df1157_CNV_MSI:", n_CNV_MSI, "(", percent_n_CNV_MSI, "% of df1157_CNV_MSI)"))

# Step 6: Further by filter(SCS == "N")

df1157_CNV_TMB_N <- df1157_CNV_TMB %>%
  filter(SCS == "N")  

# Store row count
n_CNV_TMB_N <- nrow(df1157_CNV_TMB_N)  

# Compute percentage relative to df1157_CNV
percent_n_CNV_TMB_N<- round((n_CNV_TMB_N / nrow(df1157_CNV_TMB)) * 100, 2)

# Display summary
print(paste("Rows in df1157_CNV_TMB_N:", n_CNV_TMB_N, "(", percent_n_CNV_TMB_N, "% of df1157_CNV_TMB_N)"))

# Step 6a: Further by filter(SCS == "P")

df1157_CNV_TMB_P <- df1157_CNV_TMB %>%
  filter(SCS == "P")  

# Store row count
n_CNV_TMB_P <- nrow(df1157_CNV_TMB_P)  

# Compute percentage relative to df1157_CNV
percent_n_CNV_TMB_P <- round((n_CNV_TMB_P / nrow(df1157_CNV_TMB)) * 100, 2)

# Display summary
print(paste("Rows in df1157_CNV_TMB_P:", n_CNV_TMB_P, "(", percent_n_CNV_TMB_P, "% of df1157_CNV_TMB_N)"))

# Step 7: Further filter rows that contain "Risky_count" == 4

df1157_CNV_TMB_R_4 <- df1157_CNV_TMB %>%
  filter(Risky_Count == 4)  

# Store row count
n_CNV_TMB_R_4 <- nrow(df1157_CNV_TMB_R_4)  

# Compute percentage relative to df1157_CNV
percent_n_CNV_TMB_R_4<- round((n_CNV_TMB_R_4 / nrow(df1157_CNV_TMB)) * 100, 2)

# Display summary
print(paste("Rows in df1157_CNV_TMB_R_4:", n_CNV_TMB_R_4, "(", percent_n_CNV_TMB_R_4, "% of df1157_CNV_TMB)"))

# Step 7a: Further filter rows that contain "Protective_count" == 4

df1157_CNV_TMB_P_4 <- df1157_CNV_TMB %>%
  filter(Protective_Count == 4)  

# Store row count
n_CNV_TMB_P_4 <- nrow(df1157_CNV_TMB_P_4)  

# Compute percentage relative to df1157_CNV
percent_n_CNV_TMB_P_4<- round((n_CNV_TMB_P_4 / nrow(df1157_CNV_TMB)) * 100, 2)

# Display summary
print(paste("Rows in df1157_CNV_TMB_P_4:", n_CNV_TMB_P_4, "(", percent_n_CNV_TMB_P_4, "% of df1157_CNV_TMB)"))

# Step 8: Filter df1157_CNV_TMB rows that contain "microenvironment_classification" == "pro-tumoral"

df1157_CNV_TMB_pro <- df1157_CNV_TMB %>%
  filter(microenvironment_classification == "pro-tumoral")  

# Store row count
n_CNV_TMB_pro <- nrow(df1157_CNV_TMB_pro)  

# Compute percentage relative to df1157_CNV
percent_n_CNV_TMB_pro<- round((n_CNV_TMB_pro / nrow(df1157_CNV_TMB)) * 100, 2)

# Display summary
print(paste("Rows in df1157_CNV_TMB_pro:", n_CNV_TMB_pro, "(", percent_n_CNV_TMB_pro, "% of df1157_CNV_TMB_pro)"))

# Step 9: Filter df1157_CNV_TMB rows that contain "microenvironment_classification" == "anti-tumoral"

df1157_CNV_TMB_anti <- df1157_CNV_TMB %>%
  filter(microenvironment_classification == "anti-tumoral")  

# Store row count
n_CNV_TMB_anti <- nrow(df1157_CNV_TMB_anti)  

# Compute percentage relative to df1157_CNV
percent_n_CNV_TMB_anti<- round((n_CNV_TMB_anti / nrow(df1157_CNV_TMB)) * 100, 2)

# Display summary
print(paste("Rows in df1157_CNV_TMB_anti:", n_CNV_TMB_anti, "(", percent_n_CNV_TMB_anti, "% of df1157_CNV_TMB_anti)"))

# Step 10: Filter df1157_CNV_TMB rows that contain "microenvironment_classification" == "dual"

df1157_CNV_TMB_dual <- df1157_CNV_TMB %>%
  filter(microenvironment_classification == "dual")  

# Store row count
n_CNV_TMB_dual <- nrow(df1157_CNV_TMB_dual)  

# Compute percentage relative to df1157_CNV
percent_n_CNV_TMB_dual<- round((n_CNV_TMB_dual / nrow(df1157_CNV_TMB)) * 100, 2)

# Display summary
print(paste("Rows in df1157_CNV_TMB_dual:", n_CNV_TMB_dual, "(", percent_n_CNV_TMB_dual, "% of df1157_CNV_TMB)"))

# Step 11: filter df1157_CNV_TMB rows that contain "immune_classification" == "Hot"

df1157_CNV_TMB_Hot <- df1157_CNV_TMB %>%
  filter(immune_classification == "Hot")  

# Store row count
n_CNV_TMB_Hot <- nrow(df1157_CNV_TMB_Hot)  

# Compute percentage relative to df1157_CNV
percent_n_CNV_TMB_Hot<- round((n_CNV_TMB_Hot / nrow(df1157_CNV_TMB)) * 100, 2)

# Display summary
print(paste("Rows in df1157_CNV_TMB_Hot:", n_CNV_TMB_Hot, "(", percent_n_CNV_TMB_Hot, "% of df1157_CNV_TMB)"))

# Step 12: filter df1157_CNV_TMB rows that contain "immune_classification" == "Variable"

df1157_CNV_TMB_Variable <- df1157_CNV_TMB %>%
  filter(immune_classification == "Variable")  

# Store row count
n_CNV_TMB_Variable <- nrow(df1157_CNV_TMB_Variable)  

# Compute percentage relative to df1157_CNV
percent_n_CNV_TMB_Variable<- round((n_CNV_TMB_Variable / nrow(df1157_CNV_TMB)) * 100, 2)

# Display summary
print(paste("Rows in df1157_CNV_TMB_Variable:", n_CNV_TMB_Variable, "(", percent_n_CNV_TMB_Variable, "% of df1157_CNV_TMB)"))

# Step 13: filter df1157_CNV rows that contain "immune_classification" == "Cold"

df1157_CNV_Cold <- df1157_CNV %>%
  filter(immune_classification == "Cold")  

# Store row count
n_CNV_Cold <- nrow(df1157_CNV_Cold)  

# Compute percentage relative to df1157_CNV
percent_n_CNV_Cold<- round((n_CNV_Cold / nrow(df1157_CNV)) * 100, 2)

# Display summary
print(paste("Rows in df1157_CNV_Cold:", n_CNV_Cold, "(", percent_n_CNV_Cold, "% of df1157_CNV_anti)"))

# Step 14: filter df1157_CNV rows that contain "immune_classification" == "Variable"

df1157_CNV_Variable <- df1157_CNV %>%
  filter(immune_classification == "Variable")  

# Store row count
n_CNV_Variable <- nrow(df1157_CNV_Variable)  

# Compute percentage relative to df1157_CNV
percent_n_CNV_Variable<- round((n_CNV_Variable / nrow(df1157_CNV)) * 100, 2)

# Display summary
print(paste("Rows in df1157_CNV_Variable:", n_CNV_Variable, "(", percent_n_CNV_Variable, "% of df1157_CNV_anti)"))

# Step 15: filter df1157_CNV rows that contain "microenvironment_classification" == "anti-tumoral"

df1157_CNV_anti <- df1157_CNV %>%
  filter(microenvironment_classification == "anti-tumoral")  

# Store row count
n_CNV_anti <- nrow(df1157_CNV_anti)  

# Compute percentage relative to df1157_CNV
percent_n_CNV_anti<- round((n_CNV_anti / nrow(df1157_CNV)) * 100, 2)

# Display summary
print(paste("Rows in df1157_CNV_anti:", n_CNV_anti, "(", percent_n_CNV_anti, "% of df1157_CNV_anti)"))

# Step 16: filter df1157_CNV rows that contain "microenvironment_classification" == "dual"

df1157_CNV_dual <- df1157_CNV %>%
  filter(microenvironment_classification == "dual")  

# Store row count
n_CNV_dual <- nrow(df1157_CNV_dual)  

# Compute percentage relative to df1157_CNV
percent_n_CNV_dual<- round((n_CNV_dual / nrow(df1157_CNV)) * 100, 2)

# Display summary
print(paste("Rows in df1157_CNV_dual:", n_CNV_dual, "(", percent_n_CNV_dual, "% of df1157_CNV_anti)"))

# Step 17: filter df1157_CNV rows that contain "immune_classification" == "Hot"

df1157_CNV_Hot <- df1157_CNV %>%
  filter(immune_classification == "Hot")  

# Store row count
n_CNV_Hot <- nrow(df1157_CNV_Hot)  

# Compute percentage relative to df1157_CNV
percent_n_CNV_Hot<- round((n_CNV_Hot / nrow(df1157_CNV)) * 100, 2)

# Display summary
print(paste("Rows in df1157_CNV_Hot:", n_CNV_Hot, "(", percent_n_CNV_Hot, "% of df1157_CNV_anti)"))

# Step 18: filter df1157_CNV rows that contain "immune_classification" == "Cold"

df1157_CNV_Cold <- df1157_CNV %>%
  filter(immune_classification == "Cold")  

# Store row count
n_CNV_Cold <- nrow(df1157_CNV_Cold)  

# Compute percentage relative to df1157_CNV
percent_n_CNV_Cold<- round((n_CNV_Cold / nrow(df1157_CNV)) * 100, 2)

# Display summary
print(paste("Rows in df1157_CNV_Cold:", n_CNV_Cold, "(", percent_n_CNV_Cold, "% of df1157_CNV_anti)"))

# Step 19: filter df1157_CNV rows that contain "immune_classification" == "Variable"

df1157_CNV_Variable <- df1157_CNV %>%
  filter(immune_classification == "Variable")  

# Store row count
n_CNV_Variable <- nrow(df1157_CNV_Variable)  

# Compute percentage relative to df1157_CNV
percent_n_CNV_Variable<- round((n_CNV_Variable / nrow(df1157_CNV)) * 100, 2)

# Display summary
print(paste("Rows in df1157_CNV_Variable:", n_CNV_Variable, "(", percent_n_CNV_Variable, "% of df1157_CNV_anti)"))

# Step 20: Further filter rows that contain "Phenotype" == "TSM" and Risky_Count == 4 or Protective_Count == 4

df1157_CNV_TSM_R_P_4 <- df1157_CNV %>%
  filter(Phenotype == "TSM" & (Risky_Count == 4 | Protective_Count == 4))

# Store row count
n_CNV_TSM_R_P_4 <- nrow(df1157_CNV_TSM_R_P_4)  

# Compute percentage relative to df1157_CNV
percent_n_CNV_TSM_R_P_4<- round((n_CNV_TSM_R_P_4 / nrow(df1157_CNV_TSM)) * 100, 2)

# Display summary
print(paste("Rows in df1157_CNV_TSM_R_P_4:", n_CNV_TSM_R_P_4, "(", percent_n_CNV_TSM_R_P_4, "% of df1157_CNV_TSM_R_P_4)"))

##tep 8: Further filter rows that contain "microenvironment_classification" == "anti-tumoral"

df1157_CNV_TSM_R_P_4_anti <- df1157_CNV_TSM_R_P_4 %>%
  filter(microenvironment_classification == "anti-tumoral")

# Store row count
n_CNV_TSM_R_P_4_anti <- nrow(df1157_CNV_TSM_R_P_4_anti)  

# Compute percentage relative to df1157_CNV
percent_n_CNV_TSM_R_P_4_anti<- round((n_CNV_TSM_R_P_4_anti / nrow(df1157_CNV_TSM_R_P_4)) * 100, 2)

# Display summary
print(paste("Rows in df1157_CNV_TSM_R_P_4_anti:", n_CNV_TSM_R_P_4_anti, "(", percent_n_CNV_TSM_R_P_4_anti, "% of df1157_CNV_TMB_R_P_4)"))

## Step 21: Further filter rows that contain "microenvironment_classification" == "pro-tumoral"

df1157_CNV_TSM_R_P_4_pro <- df1157_CNV_TSM_R_P_4 %>%
  filter(microenvironment_classification == "pro-tumoral")

# Store row count
n_CNV_TSM_R_P_4_pro <- nrow(df1157_CNV_TSM_R_P_4_pro)  

# Compute percentage relative to df1157_CNV
percent_n_CNV_TSM_R_P_4_pro<- round((n_CNV_TSM_R_P_4_pro / nrow(df1157_CNV_TSM_R_P_4)) * 100, 2)

# Display summary
print(paste("Rows in df1157_CNV_TSM_R_P_4_pro:", n_CNV_TSM_R_P_4_pro, "(", percent_n_CNV_TSM_R_P_4_pro, "% of df1157_CNV_TSM_R_P_4)"))

## Step 22: Further filter rows that contain "microenvironment_classification" == "dual"

df1157_CNV_TSM_R_P_4_dual <- df1157_CNV_TSM_R_P_4 %>%
  filter(microenvironment_classification == "dual")

# Store row count
n_CNV_TSM_R_P_4_dual <- nrow(df1157_CNV_TSM_R_P_4_dual)  

# Compute percentage relative to df1157_CNV
percent_n_CNV_TSM_R_P_4_dual<- round((n_CNV_TSM_R_P_4_dual / nrow(df1157_CNV_TSM_R_P_4)) * 100, 2)

# Display summary
print(paste("Rows in df1157_CNV_TSM_R_P_4_dual:", n_CNV_TSM_R_P_4_dual, "(", percent_n_CNV_TSM_R_P_4_dual, "% of df1157_CNV_TSM_R_P_4)"))

## Step 23: Further filter rows that contain "immune_classification" == "Hot"

df1157_CNV_TSM_R_P_4_Hot <- df1157_CNV_TSM_R_P_4 %>%
  filter(immune_classification == "Hot")

# Store row count
n_CNV_TSM_R_P_4_Hot <- nrow(df1157_CNV_TSM_R_P_4_Hot)  

# Compute percentage relative to  <- _CNV
percent_n_CNV_TSM_R_P_4_Hot<- round((n_CNV_TSM_R_P_4_Hot / nrow(df1157_CNV_TSM_R_P_4)) * 100, 2)

# Display summary
print(paste("Rows in <- _CNV_TSM_R_P_4_Hot:", n_CNV_TSM_R_P_4_Hot, "(", percent_n_CNV_TSM_R_P_4_Hot, "% of df1157_CNV_TSM_R_P_4)"))

## Step 24: Further filter rows that contain "immune_classification" == "Cold"

df1157_CNV_TSM_R_P_4_Cold <- df1157_CNV_TSM_R_P_4 %>%
  filter(immune_classification == "Cold")

# Store row count
n_CNV_TSM_R_P_4_Cold <- nrow(df1157_CNV_TSM_R_P_4_Cold)  

# Compute percentage relative to df1157_CNV
percent_n_CNV_TSM_R_P_4_Cold<- round((n_CNV_TSM_R_P_4_Cold / nrow(df1157_CNV_TSM_R_P_4)) * 100, 2)

# Display summary
print(paste("Rows in df1157_CNV_TSM_R_P_4_Cold:", n_CNV_TSM_R_P_4_Cold, "(", percent_n_CNV_TSM_R_P_4_Cold, "% of df1157_CNV_TMB_R_P_4)"))

## Step 25: Further filter rows that contain "immune_classification" == "Variable"

df1157_CNV_TSM_R_P_4_Variable <- df1157_CNV_TSM_R_P_4 %>%
  filter(immune_classification == "Variable")

# Store row count
n_CNV_TSM_R_P_4_Variable <- nrow(df1157_CNV_TSM_R_P_4_Variable)  

# Compute percentage relative to df1157_CNV
percent_n_CNV_TSM_R_P_4_Variable<- round((n_CNV_TSM_R_P_4_Variable / nrow(df1157_CNV_TSM_R_P_4)) * 100, 2)

# Display summary
print(paste("Rows in df1157_CNV_TSM_R_P_4_Variable:", n_CNV_TSM_R_P_4_Variable, "(", percent_n_CNV_TSM_R_P_4_Variable, "% of df1157_CNV_TSM_R_P_4)"))

#####
#####
#####
##### Creating most clinical meaningful signatures for validation by PRECOG analysis
#####
#####
##### 

# shell.exec("https://precog.stanford.edu/")
df1157_mRNA_R_P_4_STT <- df1157_mRNA %>%
  filter(Risky_Count == 4 | Protective_Count == 4) %>%
  filter(SMC_ranking != 0) %>%
  filter(TMC_ranking != 0)  %>%
  filter(TIC_ranking != 0)

# Store row count
n_mRNA_R_P_4_STT <- nrow(df1157_mRNA_R_P_4_STT)  

# Compute percentage relative to df1157_mRNA
percent_n_mRNA_R_P_4_STT<- round((n_mRNA_R_P_4_STT / nrow(df1157_mRNA)) * 100, 2)

# Display summary
print(paste("Rows in df1157_mRNA_R_P_4_STT:", n_mRNA_R_P_4_STT, "(", percent_n_mRNA_R_P_4_STT, "% of df1157_mRNA_R_P_4_STT)"))

df1158 <- df1157_mRNA_R_P_4_STT

# Rename variables in df1158
df1158 <- df1158 %>%
  rename(
    Rank = Ranking,
    Elements = Count_source,
    `RCD form` = RCD_types,
    `Omic feature` = Genotype
  )

rio::export(df1158, "df1158.tsv")
rio::export(df1158, "df1158.xlsx")

write_tsv(df1158, "df1158.tsv")

# Step 6: rio::export the cleaned dataframe to an Excel file
write.xlsx(df1158, "df1158.xlsx", overwrite = TRUE)

#######################
#######################
##### Analysis of top signatures for precog ####
#######################
#######################
# setwd("D:/Pré-artigo 5-optosis model/Dataframes all metrics")

data <- import("Drafts_2/gene_info_trancritos_excel2.xlsx")

data <- data %>%
  mutate_all(trimws)

miRNA_database <-  import("Drafts_2/miRNA_database.xlsx")

miRNA_database <- miRNA_database %>%
  mutate_all(trimws)

protein_info <- import("Drafts_2/protein_info.xlsx") 

protein_info <- protein_info %>%
  mutate_all(trimws)

df1158 <-df1158 %>%
  mutate_all(trimws)

df1157 <- df1157 %>%
  mutate_all(trimws)

# Step 0: Filter df1157 based on matching values in the "Nomenclature" column of df1157 and "Nomenclature" column indf1158
df1157_filtered <- df1157 %>%
  filter(Nomenclature %in%df1158$Nomenclature)

# Step 1: Copy the filtered columns "Signature" and "Nomenclature" from df1157_filtered to df38
df38 <- df1157_filtered %>%
  select(Signature, Nomenclature)

# Step 2: Create a new column "Gene_conversion" after "Signature"
df38 <- df38 %>%
  mutate(Gene_conversion = NA)  # Initialize the column with NA values

# Step 3: Populate "Gene_conversion" using left_join for better matching reliability

# For values starting with "ENST", join with the "data" dataframe on "Transcript_ID"
df38 <- df38 %>%
  left_join(data %>% select(Transcript_ID, Display_Name), by = c("Signature" = "Transcript_ID")) %>%
  mutate(Gene_conversion = ifelse(!is.na(Display_Name), Display_Name, Gene_conversion)) %>%
  select(-Display_Name)  # Remove the temporary column

# For values starting with "hsa-miR", join with the "miRNA_database" on "Mature1" and "Mature2"
df38 <- df38 %>%
  left_join(miRNA_database %>% select(Mature1, Gene_Symbol), by = c("Signature" = "Mature1")) %>%
  mutate(Gene_conversion = ifelse(!is.na(Gene_Symbol), Gene_Symbol, Gene_conversion)) %>%
  select(-Gene_Symbol)  # Remove the temporary column

df38 <- df38 %>%
  left_join(miRNA_database %>% select(Mature2, Gene_Symbol), by = c("Signature" = "Mature2")) %>%
  mutate(Gene_conversion = ifelse(!is.na(Gene_Symbol), Gene_Symbol, Gene_conversion)) %>%
  select(-Gene_Symbol)  # Remove the temporary column

# For protein values join with the "protein_info" on "Protein1" 
df38 <- df38 %>%
  left_join(protein_info %>% select(Protein1, Gene_Symbol), by = c("Signature" = "Protein1")) %>%
  mutate(Gene_conversion = ifelse(!is.na(Gene_Symbol), Gene_Symbol, Gene_conversion)) %>%
  select(-Gene_Symbol)  # Remove the temporary column

# Step 1: Extract gene symbols, transcript names, and miRNA from the "Signature" column
df39 <- df38 %>%
  # Use str_replace_all to remove parentheses and separate by "+" into individual elements
  mutate(Signature_cleaned = str_replace_all(Signature, "[()]", "")) %>%  # Remove parentheses
  separate_rows(Signature_cleaned, sep = " \\+ ") %>%  # Split by " + " into separate rows
  distinct() %>%  # Remove duplicate rows
  filter(str_detect(Signature_cleaned, "^ENST|^hsa-miR|^[A-Za-z]")) %>%  # Keep only gene symbols, ENST transcripts, or miRNA
  arrange(Signature_cleaned)  # Order alphabetically

# Step 2: Create the df39 with the unique and ordered entries
df39 <- df39 %>%
  select(Signature_cleaned) %>%
  rename(Extracted_Signature = Signature_cleaned)

# Step 1: Create a new column "Gene_Symbol" in df39 and populate it based on the mapping logic
df39 <- df39 %>%
  rowwise() %>%
  mutate(Gene_Symbol = case_when(
    
    # For values with prefix "ENST", look up in the "data" dataframe using "Transcript_ID"
    str_detect(Extracted_Signature, "^ENST") ~ 
      ifelse(Extracted_Signature %in% data$Transcript_ID, 
             data$Display_Name[which(data$Transcript_ID == Extracted_Signature)], 
             NA_character_),
    
    # For values with prefix "hsa-miR", look up in "miRNA_database" for Mature1
    str_detect(Extracted_Signature, "^hsa-miR") & Extracted_Signature %in% miRNA_database$Mature1 ~ 
      ifelse(Extracted_Signature %in% miRNA_database$Mature1, 
             miRNA_database$Gene_Symbol[which(miRNA_database$Mature1 == Extracted_Signature)], 
             NA_character_),
    
    # For values with prefix "hsa-miR", look up in "miRNA_database" for Mature2
    str_detect(Extracted_Signature, "^hsa-miR") & Extracted_Signature %in% miRNA_database$Mature2 ~ 
      ifelse(Extracted_Signature %in% miRNA_database$Mature2, 
             miRNA_database$Gene_Symbol[which(miRNA_database$Mature2 == Extracted_Signature)], 
             NA_character_),
    
    # If no match is found, leave as NA
    TRUE ~ NA_character_
  ))

# Fill NA values in "Gene_Symbol" with corresponding values from "Extracted_Signature"
df39 <- df39 %>%
  mutate(Gene_Symbol = ifelse(is.na(Gene_Symbol), Extracted_Signature, Gene_Symbol))

# Remove duplicates based on "Gene_Symbol"
df39 <- df39 %>%
  distinct(Gene_Symbol, .keep_all = TRUE)

# Order the dataframe by "Gene_Symbol" alphabetically
df39 <- df39 %>%
  arrange(Gene_Symbol)

# Step 1: Identify the common values between "Protein1" and "Gene_symbol"
common_proteins <- intersect(protein_info$Protein1, df39$Gene_Symbol)

# Step 2: Create a lookup table for replacements
protein_lookup <- protein_info %>%
  filter(Protein1 %in% common_proteins) %>%
  select(Protein1) %>%
  distinct()

# Step 3: Replace the "Gene_symbol" values in df39 with the corresponding "Protein1" values from protein_info
df39 <- df39 %>%
  mutate(Gene_Symbol = ifelse(Gene_Symbol %in% protein_info$Protein1,  # Check if Gene_symbol is in Protein1
                              protein_info$Gene_Symbol[match(Gene_Symbol, protein_info$Protein1)],  # Replace with corresponding Gene_symbol in protein_info
                              Gene_Symbol))  # Keep original Gene_symbol if no match

# Exclude the "Extracted_Signature" column
df39 <- df39 %>%
  select(-Extracted_Signature)

# Assign dataframe to Precog_Gene_Symbols_top_signatures_df1158
Precog_Gene_Symbols_top_signatures_df1158 <- df39

# Order the dataframe alphabetically by the "Gene_Symbol" column
Precog_Gene_Symbols_top_signatures_df1158 <- Precog_Gene_Symbols_top_signatures_df1158[order(Precog_Gene_Symbols_top_signatures_df1158$Gene_Symbol), ]

# Remove leading and trailing spaces globally across all columns
Precog_Gene_Symbols_top_signatures_df1158 <- Precog_Gene_Symbols_top_signatures_df1158 %>%
  mutate_all(trimws)

#  Remove duplicate rows
duplicates_exist <- any(duplicated(Precog_Gene_Symbols_top_signatures_df1158))

if (duplicates_exist) {
  cat("Duplicate rows found. Removing duplicates...\n")
  
  # Remove duplicate rows while retaining one instance
  Precog_Gene_Symbols_top_signatures_df1158 <- Precog_Gene_Symbols_top_signatures_df1158 %>%
    distinct()
  
  cat("Duplicates removed. Updated dataframe:\n")
  print(Precog_Gene_Symbols_top_signatures_df1158)
} else {
  cat("No duplicate rows found in the dataframe.\n")
}

# Ensure the dataframes have been loaded or created: Target_genes and df39

# Filter Target_genes based on matching values in Gene_Symbol from df39
filtered_Target_genes <- Target_genes[Target_genes$Gene_Symbol %in% df39$Gene_Symbol, ]

# Identify rows in df39 where Gene_Symbol is NOT in Target_genes$Gene_Symbol
non_matching_rows <- df39[!df39$Gene_Symbol %in% Target_genes$Gene_Symbol, ]

# View the non-matching rows
head(non_matching_rows)

# Load required libraries
library(readr)    # For writing TSV files
library(openxlsx) # For writing Excel files

# Step 5: rio::export the cleaned dataframe to a TSV file
write_tsv(Precog_Gene_Symbols_top_signatures_df1158, "Precog_Gene_Symbols_top_signatures_df1158.tsv")

# Step 6: rio::export the cleaned dataframe to an Excel file
write.xlsx(Precog_Gene_Symbols_top_signatures_df1158, "Precog_Gene_Symbols_top_signatures_df1158.xlsx", overwrite = TRUE)


# Step 5: rio::export the cleaned dataframe to a TSV file
rio::export(Precog_Gene_Symbols_top_signatures_df1158, "Precog_Gene_Symbols_top_signatures_df1158.tsv")
rio::export(Precog_Gene_Symbols_top_signatures_df1158, "Precog_Gene_Symbols_top_signatures_df1158.xlsx")

cat("Cleaned dataframe rio::exported to 'Precog_Gene_Symbols_top_signatures_df1158.tsv'.\n")

####
####
###### Checking signatures nomenclature and purpose for examples in the manuscript
####
####

df1157_POU5F1 <- df1157 %>%
  filter(grepl("\\bPOU5F1\\b", Signature)) %>%
  filter(Phenotype == "TSM") %>%
  filter(`Omic feature` == "mRNA") %>%
  filter(SCS == "P")

df1157_SOX2 <- df1157 %>%
  filter(grepl("\\bSOX2\\b", Signature)) %>%
  filter(Phenotype == "TSM") %>%
  filter(`Omic feature` == "mRNA") %>%
  filter(SCS == "P")

df1157_NANOG <- df1157 %>%
  filter(grepl("\\bNANOG\\b", Signature)) %>%
  filter(Phenotype == "TSM") %>%
  filter(`Omic feature` == "mRNA") %>%
  filter(SCS == "P")

df1157_TP53 <- df1157 %>%
  filter(grepl("\\bTP53\\b", Signature)) %>%
  filter(Phenotype == "TSM") %>%
  filter(`Omic feature` == "CNV") %>%
  filter(SCS == "N")

df1157_CDH1 <- df1157 %>%
  filter(grepl("\\bCDH1\\b", Signature)) %>%
  filter(Expression == "Overexpression")

df1157_SLC7A11_LGG <- df1157 %>%
  filter(grepl("\\bSLC7A11\\b", Signature))%>%
  filter(CTAB == "LGG") 

df1157_SLC7A11_UCEC <- df1157 %>%
  filter(grepl("\\bSLC7A11\\b", Signature))%>%
  filter(CTAB == "UCEC") 

##### 
##### 
##### Distribution of Omic feature  values with CTAB look up for the CTAB with maximum signatures per `genotyoe`Omic feature``
##### Top Cancer Type for Each `Omic feature`
##### 
##### 
# Group by CTAB and`Omic feature`, then count occurrences
df_distribution <- df1157 %>%
  group_by(CTAB, `Omic feature`) %>%
  summarise(Absolute_Count = n(), .groups = 'drop')

# Determine the maximum count per Genotype and the corresponding CTAB
df_max_counts <- df_distribution %>%
  group_by(`Omic feature`) %>%
  filter(Absolute_Count == max(Absolute_Count)) %>%
  ungroup() %>%
  select(CTAB, `Omic feature`, Absolute_Count) %>%
  rename(Absolute_maximum_Count = Absolute_Count) %>%
  arrange(desc(Absolute_maximum_Count))  # Sort in descending order

# Display the final dataframe
print(df_max_counts)

rio::export(df_max_counts, "Top Cancer Type for Each Genotype.xlsx")

### Renaming variable in df1157 for supplemental global dataset

# Rename columns in df1157
df1157 <- df1157 %>%
  rename(
    Rank = Ranking,
    Elements = Count_source,
    `Omic feature` = Genotype,
    `RCD count` = RCD,
    `RCD form` = RCD_types
  )

rio::export(df1157, "df1157.tsv")

saveRDS(df1157, "Dataset_S1_full.rds")

rio::export(df1157, "Dataset_S1_full.tsv")

####
####
####
# Reanalyzing sub-setting df1157 by `Omic feature` values to created omic feature-specifis dataset_S1_
#### Reanalyzing df1157 by `Omic feature` values to create omic feature-specific datasets ####
library(dplyr)
library(openxlsx)
library(readr)

# Step 1: Get unique "Omic feature" values
unique_omic_feature_df1157 <- unique(df1157$`Omic feature`)

# Step 2: Initialize an empty list to store objects in the environment
df1157_subsets <- list()

# Step 3: Loop through each unique `Omic feature` to create and save separate datasets
for (omic_feature in unique_omic_feature_df1157) {
  # Create a filtered dataframe for the current `Omic feature`
  Omic_feature_df <- df1157 %>% filter(`Omic feature` == omic_feature)
  
  # Store dataframe in the list using `Omic feature` as the key
  df1157_subsets[[omic_feature]] <- Omic_feature_df
  
  # Define file names for both .xlsx and .tsv formats
  file_name_xlsx <- paste0("Dataset_S1_", omic_feature, ".xlsx")
  file_name_tsv  <- paste0("Dataset_S1_", omic_feature, ".tsv")
  
  # Save dataframe to an .xlsx file
  write.xlsx(Omic_feature_df, file = file_name_xlsx, overwrite = TRUE)
  
  # Save dataframe to a .tsv file
  write_tsv(Omic_feature_df, file = file_name_tsv)
  
  # Print confirmation message
  print(paste("Files created:", file_name_xlsx, "and", file_name_tsv))
}

# Step 4: Assign the list of subsets to the global environment
assign("df1157_subsets", df1157_subsets, envir = .GlobalEnv)

# Print final message
print("All `Omic feature`-specific files and objects created successfully!")
print("All subsets are stored in the `df1157_subsets` list in the environment.")

######
######
######
###### Saving the Dataset_S1_ per omic feature as excel safely becaise of cell containing over 30000 charatecters
library(openxlsx)
library(rio)  # For `import()` function

# Define dataset names (must match file names)
dataset_names <- c("CNV", "miRNA", "Methylation", "Protein", "Transcript", "Mutation", "mRNA")

# Loop through dataset names, import them, and assign them as separate objects
for (name in dataset_names) {
  # Construct file name
  file_path <- paste0("Dataset_S1_", name, ".tsv")
  
  # Import data
  df <- import(file_path)
  
  # Dynamically assign dataset to its original name in the global environment
  assign(paste0("Dataset_S1_", name), df, envir = .GlobalEnv)
  
  # Print confirmation
  print(paste("Loaded:", file_path, "as Dataset_S1_", name))
}

# Final confirmation message
print("All datasets have been successfully loaded as separate objects in the environment.")

# Assuming Dataset_S1_CNV is already loaded with dimensions 2442 x 58
df <- Dataset_S1_CNV  # Replace this with your actual dataframe

# Define the chunk size limit for "Signature"
chunk_limit <- 30000

# Check for values in "Signature" that exceed the limit and count them
long_values_count <- sum(nchar(df$Signature) > chunk_limit)

# Create a new workbook
wb <- createWorkbook()

# Create the main sheet and write headers
addWorksheet(wb, "Main_Data")
writeData(wb, "Main_Data", df[0, ], startRow = 1, colNames = TRUE)  # Write only headers initially

if (long_values_count > 0) {
  # Notify about the count of values that exceed the character limit
  message(paste("There are", long_values_count, "values in the 'Signature' column exceeding", chunk_limit, "characters."))
  
  # Loop through each row to check the "Signature" column
  for (i in 1:nrow(df)) {
    if (nchar(df$Signature[i]) > chunk_limit) {
      # Split the "Signature" value into chunks of 30,000 characters to be safer
      signature_chunks <- substring(df$Signature[i], 
                                    seq(1, nchar(df$Signature[i]), 30000), 
                                    seq(30000, nchar(df$Signature[i]), 30000))
      
      # Create a separate sheet for this row's chunks
      chunk_sheet_name <- paste0("Signature_", i)
      addWorksheet(wb, chunk_sheet_name)
      
      # Write headers in the new chunk sheet
      writeData(wb, chunk_sheet_name, df[0, ], startRow = 1, colNames = TRUE)
      
      # Write the chunked "Signature" values in the new sheet, starting from the second row
      for (j in seq_along(signature_chunks)) {
        writeData(wb, chunk_sheet_name, signature_chunks[j], startRow = j + 1, startCol = which(names(df) == "Signature"))
      }
      
      # Write a text reference in the main sheet pointing to the chunk sheet
      reference_text <- paste("See", chunk_sheet_name)
      writeData(wb, "Main_Data", df[i, -which(names(df) == "Signature")], startRow = i + 1, colNames = FALSE)  # Write other columns
      writeData(wb, "Main_Data", reference_text, startRow = i + 1, startCol = which(names(df) == "Signature"))
      
    } else {
      # Write rows directly to main sheet if "Signature" is within the limit
      writeData(wb, "Main_Data", df[i, ], startRow = i + 1, colNames = FALSE)
    }
  }
  
} else {
  # Notify if no values exceed the limit
  message("No values in the 'Signature' column exceed the character limit.")
  
  # Write the entire dataframe directly to the main sheet
  writeData(wb, "Main_Data", df, startRow = 1, colNames = TRUE)
}

# Save the workbook regardless of conditions
saveWorkbook(wb, "Dataset_S1_CNV_safe.xlsx", overwrite = TRUE)

###
###
# Assuming Dataset_S1_Methylation is already loaded with its respective dimensions
df <- Dataset_S1_Methylation  # Use the target dataframe

# Define the chunk size limit for "Signature"
chunk_limit <- 30000

# Check for values in "Signature" that exceed the limit and count them
long_values_count <- sum(nchar(df$Signature) > chunk_limit)

# Create a new workbook
wb <- createWorkbook()

# Create the main sheet and write headers
addWorksheet(wb, "Main_Data")
writeData(wb, "Main_Data", df[0, ], startRow = 1, colNames = TRUE)  # Write only headers initially

if (long_values_count > 0) {
  # Notify about the count of values that exceed the character limit
  message(paste("There are", long_values_count, "values in the 'Signature' column exceeding", chunk_limit, "characters."))
  
  # Loop through each row to check the "Signature" column
  for (i in 1:nrow(df)) {
    if (nchar(df$Signature[i]) > chunk_limit) {
      # Split the "Signature" value into chunks of 30,000 characters
      signature_chunks <- substring(df$Signature[i], 
                                    seq(1, nchar(df$Signature[i]), chunk_limit), 
                                    seq(chunk_limit, nchar(df$Signature[i]), chunk_limit))
      
      # Create a separate sheet for this row's chunks
      chunk_sheet_name <- paste0("Signature_", i)
      addWorksheet(wb, chunk_sheet_name)
      
      # Write headers in the new chunk sheet
      writeData(wb, chunk_sheet_name, df[0, ], startRow = 1, colNames = TRUE)
      
      # Write the chunked "Signature" values in the new sheet, starting from the second row
      for (j in seq_along(signature_chunks)) {
        writeData(wb, chunk_sheet_name, signature_chunks[j], startRow = j + 1, startCol = which(names(df) == "Signature"))
      }
      
      # Write a text reference in the main sheet pointing to the chunk sheet
      reference_text <- paste("See", chunk_sheet_name)
      writeData(wb, "Main_Data", df[i, -which(names(df) == "Signature")], startRow = i + 1, colNames = FALSE)  # Write other columns
      writeData(wb, "Main_Data", reference_text, startRow = i + 1, startCol = which(names(df) == "Signature"))
      
    } else {
      # Write rows directly to main sheet if "Signature" is within the limit
      writeData(wb, "Main_Data", df[i, ], startRow = i + 1, colNames = FALSE)
    }
  }
  
} else {
  # Notify if no values exceed the limit
  message("No values in the 'Signature' column exceed the character limit.")
  
  # Write the entire dataframe directly to the main sheet
  writeData(wb, "Main_Data", df, startRow = 1, colNames = TRUE)
}

# Save the workbook regardless of conditions
saveWorkbook(wb, "Dataset_S1_Methylation_safe.xlsx", overwrite = TRUE)

##
##
# Assuming Dataset_S1_miRNA is already loaded with its respective dimensions
df <- Dataset_S1_miRNA  # Use the target dataframe

# Define the chunk size limit for "Signature"
chunk_limit <- 30000

# Check for values in "Signature" that exceed the limit and count them
long_values_count <- sum(nchar(df$Signature) > chunk_limit)

# Create a new workbook
wb <- createWorkbook()

# Create the main sheet and write headers
addWorksheet(wb, "Main_Data")
writeData(wb, "Main_Data", df[0, ], startRow = 1, colNames = TRUE)  # Write only headers initially

if (long_values_count > 0) {
  # Notify about the count of values that exceed the character limit
  message(paste("There are", long_values_count, "values in the 'Signature' column exceeding", chunk_limit, "characters."))
  
  # Loop through each row to check the "Signature" column
  for (i in 1:nrow(df)) {
    if (nchar(df$Signature[i]) > chunk_limit) {
      # Split the "Signature" value into chunks of 30,000 characters
      signature_chunks <- substring(df$Signature[i], 
                                    seq(1, nchar(df$Signature[i]), chunk_limit), 
                                    seq(chunk_limit, nchar(df$Signature[i]), chunk_limit))
      
      # Create a separate sheet for this row's chunks
      chunk_sheet_name <- paste0("Signature_", i)
      addWorksheet(wb, chunk_sheet_name)
      
      # Write headers in the new chunk sheet
      writeData(wb, chunk_sheet_name, df[0, ], startRow = 1, colNames = TRUE)
      
      # Write the chunked "Signature" values in the new sheet, starting from the second row
      for (j in seq_along(signature_chunks)) {
        writeData(wb, chunk_sheet_name, signature_chunks[j], startRow = j + 1, startCol = which(names(df) == "Signature"))
      }
      
      # Write a text reference in the main sheet pointing to the chunk sheet
      reference_text <- paste("See", chunk_sheet_name)
      writeData(wb, "Main_Data", df[i, -which(names(df) == "Signature")], startRow = i + 1, colNames = FALSE)  # Write other columns
      writeData(wb, "Main_Data", reference_text, startRow = i + 1, startCol = which(names(df) == "Signature"))
      
    } else {
      # Write rows directly to main sheet if "Signature" is within the limit
      writeData(wb, "Main_Data", df[i, ], startRow = i + 1, colNames = FALSE)
    }
  }
  
} else {
  # Notify if no values exceed the limit
  message("No values in the 'Signature' column exceed the character limit.")
  
  # Write the entire dataframe directly to the main sheet
  writeData(wb, "Main_Data", df, startRow = 1, colNames = TRUE)
}

# Save the workbook regardless of conditions
saveWorkbook(wb, "Dataset_S1_miRNA_safe.xlsx", overwrite = TRUE)

###
###
# Assuming Dataset_S1_mRNA is already loaded with its respective dimensions
df <- ######
###### saving the Dataset_S1_ as excel safely becaise of cell containing over 30000 charatecters
library(openxlsx)
library(rio)  # For `import()` function

# Define dataset names (must match file names)
dataset_names <- c("CNV", "miRNA", "Methylation", "Protein", "Transcript", "Mutation", "mRNA")

# Loop through dataset names, import them, and assign them as separate objects
for (name in dataset_names) {
  # Construct file name
  file_path <- paste0("Dataset_S1_", name, ".tsv")
  
  # Import data
  df <- import(file_path)
  
  # Dynamically assign dataset to its original name in the global environment
  assign(paste0("Dataset_S1_", name), df, envir = .GlobalEnv)
  
  # Print confirmation
  print(paste("Loaded:", file_path, "as Dataset_S1_", name))
}

# Final confirmation message
print("All datasets have been successfully loaded as separate objects in the environment.")

# Assuming Dataset_S1_CNV is already loaded with dimensions 2442 x 58
df <- Dataset_S1_CNV  # Replace this with your actual dataframe

# Define the chunk size limit for "Signature"
chunk_limit <- 30000

# Check for values in "Signature" that exceed the limit and count them
long_values_count <- sum(nchar(df$Signature) > chunk_limit)

# Create a new workbook
wb <- createWorkbook()

# Create the main sheet and write headers
addWorksheet(wb, "Main_Data")
writeData(wb, "Main_Data", df[0, ], startRow = 1, colNames = TRUE)  # Write only headers initially

if (long_values_count > 0) {
  # Notify about the count of values that exceed the character limit
  message(paste("There are", long_values_count, "values in the 'Signature' column exceeding", chunk_limit, "characters."))
  
  # Loop through each row to check the "Signature" column
  for (i in 1:nrow(df)) {
    if (nchar(df$Signature[i]) > chunk_limit) {
      # Split the "Signature" value into chunks of 30,000 characters to be safer
      signature_chunks <- substring(df$Signature[i], 
                                    seq(1, nchar(df$Signature[i]), 30000), 
                                    seq(30000, nchar(df$Signature[i]), 30000))
      
      # Create a separate sheet for this row's chunks
      chunk_sheet_name <- paste0("Signature_", i)
      addWorksheet(wb, chunk_sheet_name)
      
      # Write headers in the new chunk sheet
      writeData(wb, chunk_sheet_name, df[0, ], startRow = 1, colNames = TRUE)
      
      # Write the chunked "Signature" values in the new sheet, starting from the second row
      for (j in seq_along(signature_chunks)) {
        writeData(wb, chunk_sheet_name, signature_chunks[j], startRow = j + 1, startCol = which(names(df) == "Signature"))
      }
      
      # Write a text reference in the main sheet pointing to the chunk sheet
      reference_text <- paste("See", chunk_sheet_name)
      writeData(wb, "Main_Data", df[i, -which(names(df) == "Signature")], startRow = i + 1, colNames = FALSE)  # Write other columns
      writeData(wb, "Main_Data", reference_text, startRow = i + 1, startCol = which(names(df) == "Signature"))
      
    } else {
      # Write rows directly to main sheet if "Signature" is within the limit
      writeData(wb, "Main_Data", df[i, ], startRow = i + 1, colNames = FALSE)
    }
  }
  
} else {
  # Notify if no values exceed the limit
  message("No values in the 'Signature' column exceed the character limit.")
  
  # Write the entire dataframe directly to the main sheet
  writeData(wb, "Main_Data", df, startRow = 1, colNames = TRUE)
}

# Save the workbook regardless of conditions
saveWorkbook(wb, "Dataset_S1_CNV_safe.xlsx", overwrite = TRUE)

###
###
# Assuming Dataset_S1_Methylation is already loaded with its respective dimensions
df <- Dataset_S1_Methylation  # Use the target dataframe

# Define the chunk size limit for "Signature"
chunk_limit <- 30000

# Check for values in "Signature" that exceed the limit and count them
long_values_count <- sum(nchar(df$Signature) > chunk_limit)

# Create a new workbook
wb <- createWorkbook()

# Create the main sheet and write headers
addWorksheet(wb, "Main_Data")
writeData(wb, "Main_Data", df[0, ], startRow = 1, colNames = TRUE)  # Write only headers initially

if (long_values_count > 0) {
  # Notify about the count of values that exceed the character limit
  message(paste("There are", long_values_count, "values in the 'Signature' column exceeding", chunk_limit, "characters."))
  
  # Loop through each row to check the "Signature" column
  for (i in 1:nrow(df)) {
    if (nchar(df$Signature[i]) > chunk_limit) {
      # Split the "Signature" value into chunks of 30,000 characters
      signature_chunks <- substring(df$Signature[i], 
                                    seq(1, nchar(df$Signature[i]), chunk_limit), 
                                    seq(chunk_limit, nchar(df$Signature[i]), chunk_limit))
      
      # Create a separate sheet for this row's chunks
      chunk_sheet_name <- paste0("Signature_", i)
      addWorksheet(wb, chunk_sheet_name)
      
      # Write headers in the new chunk sheet
      writeData(wb, chunk_sheet_name, df[0, ], startRow = 1, colNames = TRUE)
      
      # Write the chunked "Signature" values in the new sheet, starting from the second row
      for (j in seq_along(signature_chunks)) {
        writeData(wb, chunk_sheet_name, signature_chunks[j], startRow = j + 1, startCol = which(names(df) == "Signature"))
      }
      
      # Write a text reference in the main sheet pointing to the chunk sheet
      reference_text <- paste("See", chunk_sheet_name)
      writeData(wb, "Main_Data", df[i, -which(names(df) == "Signature")], startRow = i + 1, colNames = FALSE)  # Write other columns
      writeData(wb, "Main_Data", reference_text, startRow = i + 1, startCol = which(names(df) == "Signature"))
      
    } else {
      # Write rows directly to main sheet if "Signature" is within the limit
      writeData(wb, "Main_Data", df[i, ], startRow = i + 1, colNames = FALSE)
    }
  }
  
} else {
  # Notify if no values exceed the limit
  message("No values in the 'Signature' column exceed the character limit.")
  
  # Write the entire dataframe directly to the main sheet
  writeData(wb, "Main_Data", df, startRow = 1, colNames = TRUE)
}

# Save the workbook regardless of conditions
saveWorkbook(wb, "Dataset_S1_Methylation_safe.xlsx", overwrite = TRUE)

##
##
# Assuming Dataset_S1_miRNA is already loaded with its respective dimensions
df <- Dataset_S1_miRNA  # Use the target dataframe

# Define the chunk size limit for "Signature"
chunk_limit <- 30000

# Check for values in "Signature" that exceed the limit and count them
long_values_count <- sum(nchar(df$Signature) > chunk_limit)

# Create a new workbook
wb <- createWorkbook()

# Create the main sheet and write headers
addWorksheet(wb, "Main_Data")
writeData(wb, "Main_Data", df[0, ], startRow = 1, colNames = TRUE)  # Write only headers initially

if (long_values_count > 0) {
  # Notify about the count of values that exceed the character limit
  message(paste("There are", long_values_count, "values in the 'Signature' column exceeding", chunk_limit, "characters."))
  
  # Loop through each row to check the "Signature" column
  for (i in 1:nrow(df)) {
    if (nchar(df$Signature[i]) > chunk_limit) {
      # Split the "Signature" value into chunks of 30,000 characters
      signature_chunks <- substring(df$Signature[i], 
                                    seq(1, nchar(df$Signature[i]), chunk_limit), 
                                    seq(chunk_limit, nchar(df$Signature[i]), chunk_limit))
      
      # Create a separate sheet for this row's chunks
      chunk_sheet_name <- paste0("Signature_", i)
      addWorksheet(wb, chunk_sheet_name)
      
      # Write headers in the new chunk sheet
      writeData(wb, chunk_sheet_name, df[0, ], startRow = 1, colNames = TRUE)
      
      # Write the chunked "Signature" values in the new sheet, starting from the second row
      for (j in seq_along(signature_chunks)) {
        writeData(wb, chunk_sheet_name, signature_chunks[j], startRow = j + 1, startCol = which(names(df) == "Signature"))
      }
      
      # Write a text reference in the main sheet pointing to the chunk sheet
      reference_text <- paste("See", chunk_sheet_name)
      writeData(wb, "Main_Data", df[i, -which(names(df) == "Signature")], startRow = i + 1, colNames = FALSE)  # Write other columns
      writeData(wb, "Main_Data", reference_text, startRow = i + 1, startCol = which(names(df) == "Signature"))
      
    } else {
      # Write rows directly to main sheet if "Signature" is within the limit
      writeData(wb, "Main_Data", df[i, ], startRow = i + 1, colNames = FALSE)
    }
  }
  
} else {
  # Notify if no values exceed the limit
  message("No values in the 'Signature' column exceed the character limit.")
  
  # Write the entire dataframe directly to the main sheet
  writeData(wb, "Main_Data", df, startRow = 1, colNames = TRUE)
}

# Save the workbook regardless of conditions
saveWorkbook(wb, "Dataset_S1_miRNA_safe.xlsx", overwrite = TRUE)

##
##
##
# Assuming Dataset_S1_mRNA is already loaded with its respective dimensions
df <- Dataset_S1_mRNA  # Use the target dataframe

# Define the chunk size limit for "Signature"
chunk_limit <- 30000

# Check for values in "Signature" that exceed the limit and count them
long_values_count <- sum(nchar(df$Signature) > chunk_limit)

# Create a new workbook
wb <- createWorkbook()

# Create the main sheet and write headers
addWorksheet(wb, "Main_Data")
writeData(wb, "Main_Data", df[0, ], startRow = 1, colNames = TRUE)  # Write only headers initially

if (long_values_count > 0) {
  # Notify about the count of values that exceed the character limit
  message(paste("There are", long_values_count, "values in the 'Signature' column exceeding", chunk_limit, "characters."))
  
  # Loop through each row to check the "Signature" column
  for (i in 1:nrow(df)) {
    if (nchar(df$Signature[i]) > chunk_limit) {
      # Split the "Signature" value into chunks of 30,000 characters
      signature_chunks <- substring(df$Signature[i], 
                                    seq(1, nchar(df$Signature[i]), chunk_limit), 
                                    seq(chunk_limit, nchar(df$Signature[i]), chunk_limit))
      
      # Create a separate sheet for this row's chunks
      chunk_sheet_name <- paste0("Signature_", i)
      addWorksheet(wb, chunk_sheet_name)
      
      # Write headers in the new chunk sheet
      writeData(wb, chunk_sheet_name, df[0, ], startRow = 1, colNames = TRUE)
      
      # Write the chunked "Signature" values in the new sheet, starting from the second row
      for (j in seq_along(signature_chunks)) {
        writeData(wb, chunk_sheet_name, signature_chunks[j], startRow = j + 1, startCol = which(names(df) == "Signature"))
      }
      
      # Write a text reference in the main sheet pointing to the chunk sheet
      reference_text <- paste("See", chunk_sheet_name)
      writeData(wb, "Main_Data", df[i, -which(names(df) == "Signature")], startRow = i + 1, colNames = FALSE)  # Write other columns
      writeData(wb, "Main_Data", reference_text, startRow = i + 1, startCol = which(names(df) == "Signature"))
      
    } else {
      # Write rows directly to main sheet if "Signature" is within the limit
      writeData(wb, "Main_Data", df[i, ], startRow = i + 1, colNames = FALSE)
    }
  }
  
} else {
  # Notify if no values exceed the limit
  message("No values in the 'Signature' column exceed the character limit.")
  
  # Write the entire dataframe directly to the main sheet
  writeData(wb, "Main_Data", df, startRow = 1, colNames = TRUE)
}

# Save the workbook regardless of conditions
saveWorkbook(wb, "Dataset_S1_mRNA_safe.xlsx", overwrite = TRUE)

###
###
# Assuming Dataset_S1_Mutation is already loaded with its respective dimensions
df <- Dataset_S1_Mutation  # Use the target dataframe

# Define the chunk size limit for "Signature"
chunk_limit <- 30000

# Check for values in "Signature" that exceed the limit and count them
long_values_count <- sum(nchar(df$Signature) > chunk_limit)

# Create a new workbook
wb <- createWorkbook()

# Create the main sheet and write headers
addWorksheet(wb, "Main_Data")
writeData(wb, "Main_Data", df[0, ], startRow = 1, colNames = TRUE)  # Write only headers initially

if (long_values_count > 0) {
  # Notify about the count of values that exceed the character limit
  message(paste("There are", long_values_count, "values in the 'Signature' column exceeding", chunk_limit, "characters."))
  
  # Loop through each row to check the "Signature" column
  for (i in 1:nrow(df)) {
    if (nchar(df$Signature[i]) > chunk_limit) {
      # Split the "Signature" value into chunks of 30,000 characters
      signature_chunks <- substring(df$Signature[i], 
                                    seq(1, nchar(df$Signature[i]), chunk_limit), 
                                    seq(chunk_limit, nchar(df$Signature[i]), chunk_limit))
      
      # Create a separate sheet for this row's chunks
      chunk_sheet_name <- paste0("Signature_", i)
      addWorksheet(wb, chunk_sheet_name)
      
      # Write headers in the new chunk sheet
      writeData(wb, chunk_sheet_name, df[0, ], startRow = 1, colNames = TRUE)
      
      # Write the chunked "Signature" values in the new sheet, starting from the second row
      for (j in seq_along(signature_chunks)) {
        writeData(wb, chunk_sheet_name, signature_chunks[j], startRow = j + 1, startCol = which(names(df) == "Signature"))
      }
      
      # Write a text reference in the main sheet pointing to the chunk sheet
      reference_text <- paste("See", chunk_sheet_name)
      writeData(wb, "Main_Data", df[i, -which(names(df) == "Signature")], startRow = i + 1, colNames = FALSE)  # Write other columns
      writeData(wb, "Main_Data", reference_text, startRow = i + 1, startCol = which(names(df) == "Signature"))
      
    } else {
      # Write rows directly to main sheet if "Signature" is within the limit
      writeData(wb, "Main_Data", df[i, ], startRow = i + 1, colNames = FALSE)
    }
  }
  
} else {
  # Notify if no values exceed the limit
  message("No values in the 'Signature' column exceed the character limit.")
  
  # Write the entire dataframe directly to the main sheet
  writeData(wb, "Main_Data", df, startRow = 1, colNames = TRUE)
}

# Save the workbook regardless of conditions
saveWorkbook(wb, "Dataset_S1_Mutation_safe.xlsx", overwrite = TRUE)

##
##
# Assuming Dataset_S1_Protein is already loaded with its respective dimensions
df <- Dataset_S1_Protein  # Use the target dataframe

# Define the chunk size limit for "Signature"
chunk_limit <- 30000

# Check for values in "Signature" that exceed the limit and count them
long_values_count <- sum(nchar(df$Signature) > chunk_limit)

# Create a new workbook
wb <- createWorkbook()

# Create the main sheet and write headers
addWorksheet(wb, "Main_Data")
writeData(wb, "Main_Data", df[0, ], startRow = 1, colNames = TRUE)  # Write only headers initially

if (long_values_count > 0) {
  # Notify about the count of values that exceed the character limit
  message(paste("There are", long_values_count, "values in the 'Signature' column exceeding", chunk_limit, "characters."))
  
  # Loop through each row to check the "Signature" column
  for (i in 1:nrow(df)) {
    if (nchar(df$Signature[i]) > chunk_limit) {
      # Split the "Signature" value into chunks of 30,000 characters
      signature_chunks <- substring(df$Signature[i], 
                                    seq(1, nchar(df$Signature[i]), chunk_limit), 
                                    seq(chunk_limit, nchar(df$Signature[i]), chunk_limit))
      
      # Create a separate sheet for this row's chunks
      chunk_sheet_name <- paste0("Signature_", i)
      addWorksheet(wb, chunk_sheet_name)
      
      # Write headers in the new chunk sheet
      writeData(wb, chunk_sheet_name, df[0, ], startRow = 1, colNames = TRUE)
      
      # Write the chunked "Signature" values in the new sheet, starting from the second row
      for (j in seq_along(signature_chunks)) {
        writeData(wb, chunk_sheet_name, signature_chunks[j], startRow = j + 1, startCol = which(names(df) == "Signature"))
      }
      
      # Write a text reference in the main sheet pointing to the chunk sheet
      reference_text <- paste("See", chunk_sheet_name)
      writeData(wb, "Main_Data", df[i, -which(names(df) == "Signature")], startRow = i + 1, colNames = FALSE)  # Write other columns
      writeData(wb, "Main_Data", reference_text, startRow = i + 1, startCol = which(names(df) == "Signature"))
      
    } else {
      # Write rows directly to main sheet if "Signature" is within the limit
      writeData(wb, "Main_Data", df[i, ], startRow = i + 1, colNames = FALSE)
    }
  }
  
} else {
  # Notify if no values exceed the limit
  message("No values in the 'Signature' column exceed the character limit.")
  
  # Write the entire dataframe directly to the main sheet
  writeData(wb, "Main_Data", df, startRow = 1, colNames = TRUE)
}

# Save the workbook regardless of conditions
saveWorkbook(wb, "Dataset_S1_Protein_safe.xlsx", overwrite = TRUE)

##
##
## Saving Dataset_S1_Transcript as excel with values çaege than 3000 characters
##
##### Optimized Signature Data rio::export with Reduced Loop Write Operations #####
# Load data
df165 <- Dataset_S1_Transcript

# Create a new workbook
wb <- createWorkbook()

# Add the main worksheet for primary data
addWorksheet(wb, "MainData")

# Write main data (excluding Signature) to MainData sheet in one go
writeData(wb, "MainData", df165[, -which(names(df165) == "Signature")])

# Define optimized function to split long Signature strings into chunks under Excel's character limit
split_signature <- function(signature_entry, max_length = 30000) {
  ids <- strsplit(signature_entry, " \\+ ")[[1]]
  current_chunk <- ids[1]
  chunks <- vector("list")  # Predefine the list size
  
  for (i in seq_along(ids)[-1]) {  # Start loop from the second element
    if (nchar(paste(current_chunk, ids[i], sep = " + ")) <= max_length) {
      current_chunk <- paste(current_chunk, ids[i], sep = " + ")
    } else {
      chunks[[length(chunks) + 1]] <- current_chunk
      current_chunk <- ids[i]
    }
  }
  chunks[[length(chunks) + 1]] <- current_chunk  # Final chunk
  
  return(unlist(chunks))
}

# Check if overflow is needed, and if so, prepare TranscriptDetails
overflow_needed <- any(nchar(df165$Signature) > 30000)
if (overflow_needed) {
  addWorksheet(wb, "TranscriptDetails")
  writeData(wb, "TranscriptDetails", data.frame(Row_Index = "Row_Index", Nomenclature = "Nomenclature", Signature_Chunk = "Signature_Chunk"), startRow = 1)
}

# Prepare to write data in bulk by storing rows in lists before writing them all at once
main_data_hyperlinks <- vector("list", nrow(df165))
transcript_details_data <- list()
start_row <- 2  # Row in TranscriptDetails to start data (skip headers)

# Loop over each row in df165
for (i in 1:nrow(df165)) {
  if (nchar(df165$Signature[i]) > 30000) {
    chunks <- split_signature(df165$Signature[i])
    
    # Add row information for each chunk to transcript details in a bulk list
    for (chunk in chunks) {
      transcript_details_data[[length(transcript_details_data) + 1]] <- data.frame(
        Row_Index = i + 1,
        Nomenclature = df165$Nomenclature[i],
        Signature_Chunk = chunk
      )
    }
    
    # Prepare hyperlink text for MainData
    summary_text <- paste0("n=", df165$Member_counts[i], " transcript isoforms")
    main_data_hyperlinks[[i]] <- data.frame(
      row = i + 1,
      col = which(names(df165) == "Signature"),
      text = summary_text,
      formula = paste0("HYPERLINK(\"#TranscriptDetails!A", start_row, "\", \"See full list\")")
    )
    
    start_row <- start_row + length(chunks)
  } else {
    # If within limit, write Signature entry directly to MainData
    main_data_hyperlinks[[i]] <- data.frame(
      row = i + 1,
      col = which(names(df165) == "Signature"),
      text = df165$Signature[i],
      formula = NA  # No formula required for short entries
    )
  }
}

# Bulk write transcript details data to reduce function calls
if (overflow_needed && length(transcript_details_data) > 0) {
  transcript_details_df <- do.call(rbind, transcript_details_data)
  writeData(wb, "TranscriptDetails", transcript_details_df, startRow = 2, colNames = FALSE)
}

# Write hyperlinks and text data to MainData in one go
for (item in main_data_hyperlinks) {
  writeData(wb, "MainData", item$text, startRow = item$row, startCol = item$col)
  if (!is.na(item$formula)) {
    writeFormula(wb, "MainData", item$formula, startRow = item$row, startCol = item$col)
  }
}

# Save workbook
saveWorkbook(wb, "Transcript_S1_Optimized_Supplementary_Data_with_Selective_Hyperlinks.xlsx", overwrite = TRUE)

# Import saved workbook to df166 for verification
df166 <- import("Transcript_S1_Optimized_Supplementary_Data_with_Selective_Hyperlinks.xlsx")

##### Optimized Signature Data  rio::export with Reduced Loop Write Operations #####

# Load required packages
library(openxlsx)
# Assuming `import` is from `rio` package

# Load data
df165 <- Dataset_S1_Transcript

# Create a new workbook
wb <- createWorkbook()

# Add the main worksheet for primary data
addWorksheet(wb, "MainData")

# Write main data (excluding Signature) to MainData sheet in one go
writeData(wb, "MainData", df165[, -which(names(df165) == "Signature")])

# Define optimized function to split long Signature strings into chunks under Excel's character limit
split_signature <- function(signature_entry, max_length = 30000) {
  ids <- strsplit(signature_entry, " \\+ ")[[1]]
  current_chunk <- ids[1]
  chunks <- vector("list")  # Predefine the list size
  
  for (i in seq_along(ids)[-1]) {  # Start loop from the second element
    if (nchar(paste(current_chunk, ids[i], sep = " + ")) <= max_length) {
      current_chunk <- paste(current_chunk, ids[i], sep = " + ")
    } else {
      chunks[[length(chunks) + 1]] <- current_chunk
      current_chunk <- ids[i]
    }
  }
  chunks[[length(chunks) + 1]] <- current_chunk  # Final chunk
  
  return(unlist(chunks))
}

# Check if overflow is needed, and if so, prepare TranscriptDetails
overflow_needed <- any(nchar(df165$Signature) > 30000)
if (overflow_needed) {
  addWorksheet(wb, "TranscriptDetails")
  writeData(wb, "TranscriptDetails", data.frame(Row_Index = "Row_Index", Nomenclature = "Nomenclature", Signature_Chunk = "Signature_Chunk"), startRow = 1)
}

# Prepare to write data in bulk by storing rows in lists before writing them all at once
main_data_hyperlinks <- vector("list", nrow(df165))
transcript_details_data <- list()
start_row <- 2  # Row in TranscriptDetails to start data (skip headers)

# Loop over each row in df165
for (i in 1:nrow(df165)) {
  if (nchar(df165$Signature[i]) > 30000) {
    chunks <- split_signature(df165$Signature[i])
    
    # Add row information for each chunk to transcript details in a bulk list
    for (chunk in chunks) {
      transcript_details_data[[length(transcript_details_data) + 1]] <- data.frame(
        Row_Index = i + 1,
        Nomenclature = df165$Nomenclature[i],
        Signature_Chunk = chunk
      )
    }
    
    # Prepare hyperlink text for MainData
    summary_text <- paste0("n=", df165$Member_counts[i], " transcript isoforms")
    main_data_hyperlinks[[i]] <- data.frame(
      row = i + 1,
      col = which(names(df165) == "Signature"),
      text = summary_text,
      formula = paste0("HYPERLINK(\"#TranscriptDetails!A", start_row, "\", \"See full list\")")
    )
    
    start_row <- start_row + length(chunks)
  } else {
    # If within limit, write Signature entry directly to MainData
    main_data_hyperlinks[[i]] <- data.frame(
      row = i + 1,
      col = which(names(df165) == "Signature"),
      text = df165$Signature[i],
      formula = NA  # No formula required for short entries
    )
  }
}

# Bulk write transcript details data to reduce function calls
if (overflow_needed && length(transcript_details_data) > 0) {
  transcript_details_df <- do.call(rbind, transcript_details_data)
  writeData(wb, "TranscriptDetails", transcript_details_df, startRow = 2, colNames = FALSE)
}

# Write hyperlinks and text data to MainData in one go
for (item in main_data_hyperlinks) {
  writeData(wb, "MainData", item$text, startRow = item$row, startCol = item$col)
  if (!is.na(item$formula)) {
    writeFormula(wb, "MainData", item$formula, startRow = item$row, startCol = item$col)
  }
}

# Save workbook
saveWorkbook(wb, "Transcript_S1b_Optimized_Supplementary_Data_with_Selective_Hyperlinks_v02.xlsx", overwrite = TRUE)

# Import saved workbook to df166 for verification
df166_v04 <- import("Transcript_S1b_Optimized_Supplementary_Data_with_Selective_Hyperlinks_v02.xlsx")


#####AQUIAQUIAQUI
#####
#####
######
######
######
##### Optimized Signature Data rio::export with Reduced Loop Write Operations #####
######
######
######
# Load required packages
library(openxlsx)
# Assuming `import` is from `rio` package

# Load data
df165 <- Dataset_S1_Transcript

# Create a new workbook
wb <- createWorkbook()

# Add the main worksheet for primary data
addWorksheet(wb, "MainData")

# Write main data (excluding Signature) to MainData sheet in one go
writeData(wb, "MainData", df165[, -which(names(df165) == "Signature")])

# Define optimized function to split long Signature strings into chunks under Excel's character limit
split_signature <- function(signature_entry, max_length = 30000) {
  ids <- strsplit(signature_entry, " \\+ ")[[1]]
  current_chunk <- ids[1]
  chunks <- vector("list")  # Predefine the list size
  
  for (i in seq_along(ids)[-1]) {  # Start loop from the second element
    if (nchar(paste(current_chunk, ids[i], sep = " + ")) <= max_length) {
      current_chunk <- paste(current_chunk, ids[i], sep = " + ")
    } else {
      chunks[[length(chunks) + 1]] <- current_chunk
      current_chunk <- ids[i]
    }
  }
  chunks[[length(chunks) + 1]] <- current_chunk  # Final chunk
  
  return(unlist(chunks))
}

# Check if overflow is needed, and if so, prepare TranscriptDetails
overflow_needed <- any(nchar(df165$Signature) > 30000)
if (overflow_needed) {
  addWorksheet(wb, "TranscriptDetails")
  writeData(wb, "TranscriptDetails", data.frame(Row_Index = "Row_Index", Nomenclature = "Nomenclature", Signature_Chunk = "Signature_Chunk"), startRow = 1)
}

# Prepare to write data in bulk by storing rows in lists before writing them all at once
main_data_hyperlinks <- vector("list", nrow(df165))
transcript_details_data <- list()
start_row <- 2  # Row in TranscriptDetails to start data (skip headers)

# Loop over each row in df165
for (i in 1:nrow(df165)) {
  if (nchar(df165$Signature[i]) > 30000) {
    chunks <- split_signature(df165$Signature[i])
    
    # Add row information for each chunk to transcript details in a bulk list
    for (chunk in chunks) {
      transcript_details_data[[length(transcript_details_data) + 1]] <- data.frame(
        Row_Index = i + 1,
        Nomenclature = df165$Nomenclature[i],
        Signature_Chunk = chunk
      )
    }
    
    # Prepare hyperlink text for MainData
    summary_text <- paste0("n=", df165$Member_counts[i], " transcript isoforms")
    main_data_hyperlinks[[i]] <- data.frame(
      row = i + 1,
      col = which(names(df165) == "Signature"),
      text = summary_text,
      formula = paste0("HYPERLINK(\"#TranscriptDetails!A", start_row, "\", \"See full list\")")
    )
    
    start_row <- start_row + length(chunks)
  } else {
    # If within limit, write Signature entry directly to MainData
    main_data_hyperlinks[[i]] <- data.frame(
      row = i + 1,
      col = which(names(df165) == "Signature"),
      text = df165$Signature[i],
      formula = NA  # No formula required for short entries
    )
  }
}

# Bulk write transcript details data to reduce function calls
if (overflow_needed && length(transcript_details_data) > 0) {
  transcript_details_df <- do.call(rbind, transcript_details_data)
  writeData(wb, "TranscriptDetails", transcript_details_df, startRow = 2, colNames = FALSE)
}

# Write hyperlinks and text data to MainData in one go
for (item in main_data_hyperlinks) {
  writeData(wb, "MainData", item$text, startRow = item$row, startCol = item$col)
  if (!is.na(item$formula)) {
    writeFormula(wb, "MainData", item$formula, startRow = item$row, startCol = item$col)
  }
}

# Save workbook
saveWorkbook(wb, "Transcript_S1c_Optimized_Supplementary_Data_with_Selective_Hyperlinks_v02.xlsx", overwrite = TRUE)

# Import saved workbook to df166 for verification
df166_v05 <- import("Transcript_S1c_Optimized_Supplementary_Data_with_Selective_Hyperlinks_v02.xlsx")
mRNA  # Use the target dataframe

# Define the chunk size limit for "Signature"
chunk_limit <- 30000

# Check for values in "Signature" that exceed the limit and count them
long_values_count <- sum(nchar(df$Signature) > chunk_limit)

# Create a new workbook
wb <- createWorkbook()

# Create the main sheet and write headers
addWorksheet(wb, "Main_Data")
writeData(wb, "Main_Data", df[0, ], startRow = 1, colNames = TRUE)  # Write only headers initially

if (long_values_count > 0) {
  # Notify about the count of values that exceed the character limit
  message(paste("There are", long_values_count, "values in the 'Signature' column exceeding", chunk_limit, "characters."))
  
  # Loop through each row to check the "Signature" column
  for (i in 1:nrow(df)) {
    if (nchar(df$Signature[i]) > chunk_limit) {
      # Split the "Signature" value into chunks of 30,000 characters
      signature_chunks <- substring(df$Signature[i], 
                                    seq(1, nchar(df$Signature[i]), chunk_limit), 
                                    seq(chunk_limit, nchar(df$Signature[i]), chunk_limit))
      
      # Create a separate sheet for this row's chunks
      chunk_sheet_name <- paste0("Signature_", i)
      addWorksheet(wb, chunk_sheet_name)
      
      # Write headers in the new chunk sheet
      writeData(wb, chunk_sheet_name, df[0, ], startRow = 1, colNames = TRUE)
      
      # Write the chunked "Signature" values in the new sheet, starting from the second row
      for (j in seq_along(signature_chunks)) {
        writeData(wb, chunk_sheet_name, signature_chunks[j], startRow = j + 1, startCol = which(names(df) == "Signature"))
      }
      
      # Write a text reference in the main sheet pointing to the chunk sheet
      reference_text <- paste("See", chunk_sheet_name)
      writeData(wb, "Main_Data", df[i, -which(names(df) == "Signature")], startRow = i + 1, colNames = FALSE)  # Write other columns
      writeData(wb, "Main_Data", reference_text, startRow = i + 1, startCol = which(names(df) == "Signature"))
      
    } else {
      # Write rows directly to main sheet if "Signature" is within the limit
      writeData(wb, "Main_Data", df[i, ], startRow = i + 1, colNames = FALSE)
    }
  }
  
} else {
  # Notify if no values exceed the limit
  message("No values in the 'Signature' column exceed the character limit.")
  
  # Write the entire dataframe directly to the main sheet
  writeData(wb, "Main_Data", df, startRow = 1, colNames = TRUE)
}

# Save the workbook regardless of conditions
saveWorkbook(wb, "Dataset_S1_mRNA_safe.xlsx", overwrite = TRUE)

###
###
# Assuming Dataset_S1_Mutation is already loaded with its respective dimensions
df <- Dataset_S1_Mutation  # Use the target dataframe

# Define the chunk size limit for "Signature"
chunk_limit <- 30000

# Check for values in "Signature" that exceed the limit and count them
long_values_count <- sum(nchar(df$Signature) > chunk_limit)

# Create a new workbook
wb <- createWorkbook()

# Create the main sheet and write headers
addWorksheet(wb, "Main_Data")
writeData(wb, "Main_Data", df[0, ], startRow = 1, colNames = TRUE)  # Write only headers initially

if (long_values_count > 0) {
  # Notify about the count of values that exceed the character limit
  message(paste("There are", long_values_count, "values in the 'Signature' column exceeding", chunk_limit, "characters."))
  
  # Loop through each row to check the "Signature" column
  for (i in 1:nrow(df)) {
    if (nchar(df$Signature[i]) > chunk_limit) {
      # Split the "Signature" value into chunks of 30,000 characters
      signature_chunks <- substring(df$Signature[i], 
                                    seq(1, nchar(df$Signature[i]), chunk_limit), 
                                    seq(chunk_limit, nchar(df$Signature[i]), chunk_limit))
      
      # Create a separate sheet for this row's chunks
      chunk_sheet_name <- paste0("Signature_", i)
      addWorksheet(wb, chunk_sheet_name)
      
      # Write headers in the new chunk sheet
      writeData(wb, chunk_sheet_name, df[0, ], startRow = 1, colNames = TRUE)
      
      # Write the chunked "Signature" values in the new sheet, starting from the second row
      for (j in seq_along(signature_chunks)) {
        writeData(wb, chunk_sheet_name, signature_chunks[j], startRow = j + 1, startCol = which(names(df) == "Signature"))
      }
      
      # Write a text reference in the main sheet pointing to the chunk sheet
      reference_text <- paste("See", chunk_sheet_name)
      writeData(wb, "Main_Data", df[i, -which(names(df) == "Signature")], startRow = i + 1, colNames = FALSE)  # Write other columns
      writeData(wb, "Main_Data", reference_text, startRow = i + 1, startCol = which(names(df) == "Signature"))
      
    } else {
      # Write rows directly to main sheet if "Signature" is within the limit
      writeData(wb, "Main_Data", df[i, ], startRow = i + 1, colNames = FALSE)
    }
  }
  
} else {
  # Notify if no values exceed the limit
  message("No values in the 'Signature' column exceed the character limit.")
  
  # Write the entire dataframe directly to the main sheet
  writeData(wb, "Main_Data", df, startRow = 1, colNames = TRUE)
}

# Save the workbook regardless of conditions
saveWorkbook(wb, "Dataset_S1_Mutation_safe.xlsx", overwrite = TRUE)

##
##
# Assuming Dataset_S1_Protein is already loaded with its respective dimensions
df <- Dataset_S1_Protein  # Use the target dataframe

# Define the chunk size limit for "Signature"
chunk_limit <- 30000

# Check for values in "Signature" that exceed the limit and count them
long_values_count <- sum(nchar(df$Signature) > chunk_limit)

# Create a new workbook
wb <- createWorkbook()

# Create the main sheet and write headers
addWorksheet(wb, "Main_Data")
writeData(wb, "Main_Data", df[0, ], startRow = 1, colNames = TRUE)  # Write only headers initially

if (long_values_count > 0) {
  # Notify about the count of values that exceed the character limit
  message(paste("There are", long_values_count, "values in the 'Signature' column exceeding", chunk_limit, "characters."))
  
  # Loop through each row to check the "Signature" column
  for (i in 1:nrow(df)) {
    if (nchar(df$Signature[i]) > chunk_limit) {
      # Split the "Signature" value into chunks of 30,000 characters
      signature_chunks <- substring(df$Signature[i], 
                                    seq(1, nchar(df$Signature[i]), chunk_limit), 
                                    seq(chunk_limit, nchar(df$Signature[i]), chunk_limit))
      
      # Create a separate sheet for this row's chunks
      chunk_sheet_name <- paste0("Signature_", i)
      addWorksheet(wb, chunk_sheet_name)
      
      # Write headers in the new chunk sheet
      writeData(wb, chunk_sheet_name, df[0, ], startRow = 1, colNames = TRUE)
      
      # Write the chunked "Signature" values in the new sheet, starting from the second row
      for (j in seq_along(signature_chunks)) {
        writeData(wb, chunk_sheet_name, signature_chunks[j], startRow = j + 1, startCol = which(names(df) == "Signature"))
      }
      
      # Write a text reference in the main sheet pointing to the chunk sheet
      reference_text <- paste("See", chunk_sheet_name)
      writeData(wb, "Main_Data", df[i, -which(names(df) == "Signature")], startRow = i + 1, colNames = FALSE)  # Write other columns
      writeData(wb, "Main_Data", reference_text, startRow = i + 1, startCol = which(names(df) == "Signature"))
      
    } else {
      # Write rows directly to main sheet if "Signature" is within the limit
      writeData(wb, "Main_Data", df[i, ], startRow = i + 1, colNames = FALSE)
    }
  }
  
} else {
  # Notify if no values exceed the limit
  message("No values in the 'Signature' column exceed the character limit.")
  
  # Write the entire dataframe directly to the main sheet
  writeData(wb, "Main_Data", df, startRow = 1, colNames = TRUE)
}

# Save the workbook regardless of conditions
saveWorkbook(wb, "Dataset_S1_Protein_safe.xlsx", overwrite = TRUE)

##
##
##### Optimized Signature Data rio::export with Reduced Loop Write Operations #####
# Load data
df165 <- Dataset_S1_Transcript

# Create a new workbook
wb <- createWorkbook()

# Add the main worksheet for primary data
addWorksheet(wb, "MainData")

# Write main data (excluding Signature) to MainData sheet in one go
writeData(wb, "MainData", df165[, -which(names(df165) == "Signature")])

# Define optimized function to split long Signature strings into chunks under Excel's character limit
split_signature <- function(signature_entry, max_length = 30000) {
  ids <- strsplit(signature_entry, " \\+ ")[[1]]
  current_chunk <- ids[1]
  chunks <- vector("list")  # Predefine the list size
  
  for (i in seq_along(ids)[-1]) {  # Start loop from the second element
    if (nchar(paste(current_chunk, ids[i], sep = " + ")) <= max_length) {
      current_chunk <- paste(current_chunk, ids[i], sep = " + ")
    } else {
      chunks[[length(chunks) + 1]] <- current_chunk
      current_chunk <- ids[i]
    }
  }
  chunks[[length(chunks) + 1]] <- current_chunk  # Final chunk
  
  return(unlist(chunks))
}

# Check if overflow is needed, and if so, prepare TranscriptDetails
overflow_needed <- any(nchar(df165$Signature) > 30000)
if (overflow_needed) {
  addWorksheet(wb, "TranscriptDetails")
  writeData(wb, "TranscriptDetails", data.frame(Row_Index = "Row_Index", Nomenclature = "Nomenclature", Signature_Chunk = "Signature_Chunk"), startRow = 1)
}

# Prepare to write data in bulk by storing rows in lists before writing them all at once
main_data_hyperlinks <- vector("list", nrow(df165))
transcript_details_data <- list()
start_row <- 2  # Row in TranscriptDetails to start data (skip headers)

# Loop over each row in df165
for (i in 1:nrow(df165)) {
  if (nchar(df165$Signature[i]) > 30000) {
    chunks <- split_signature(df165$Signature[i])
    
    # Add row information for each chunk to transcript details in a bulk list
    for (chunk in chunks) {
      transcript_details_data[[length(transcript_details_data) + 1]] <- data.frame(
        Row_Index = i + 1,
        Nomenclature = df165$Nomenclature[i],
        Signature_Chunk = chunk
      )
    }
    
    # Prepare hyperlink text for MainData
    summary_text <- paste0("n=", df165$Elements[i], " transcript isoforms")
    main_data_hyperlinks[[i]] <- data.frame(
      row = i + 1,
      col = which(names(df165) == "Signature"),
      text = summary_text,
      formula = paste0("HYPERLINK(\"#TranscriptDetails!A", start_row, "\", \"See full list\")")
    )
    
    start_row <- start_row + length(chunks)
  } else {
    # If within limit, write Signature entry directly to MainData
    main_data_hyperlinks[[i]] <- data.frame(
      row = i + 1,
      col = which(names(df165) == "Signature"),
      text = df165$Signature[i],
      formula = NA  # No formula required for short entries
    )
  }
}

# Bulk write transcript details data to reduce function calls
if (overflow_needed && length(transcript_details_data) > 0) {
  transcript_details_df <- do.call(rbind, transcript_details_data)
  writeData(wb, "TranscriptDetails", transcript_details_df, startRow = 2, colNames = FALSE)
}

# Write hyperlinks and text data to MainData in one go
for (item in main_data_hyperlinks) {
  writeData(wb, "MainData", item$text, startRow = item$row, startCol = item$col)
  if (!is.na(item$formula)) {
    writeFormula(wb, "MainData", item$formula, startRow = item$row, startCol = item$col)
  }
}

# Save workbook
saveWorkbook(wb, "Transcript_S1d_Optimized_Supplementary_Data_with_Selective_Hyperlinks.xlsx", overwrite = TRUE)

# Import saved workbook to df166 for verification
df166 <- import("Transcript_S1d_Optimized_Supplementary_Data_with_Selective_Hyperlinks.xlsx")

##### Optimized Signature Data  rio::export with Reduced Loop Write Operations #####
# Load data
df165 <- Dataset_S1_Transcript

# Create a new workbook
wb <- createWorkbook()

# Add the main worksheet for primary data
addWorksheet(wb, "MainData")

# Write main data (excluding Signature) to MainData sheet in one go
writeData(wb, "MainData", df165[, -which(names(df165) == "Signature")])

# Define optimized function to split long Signature strings into chunks under Excel's character limit
split_signature <- function(signature_entry, max_length = 30000) {
  ids <- strsplit(signature_entry, " \\+ ")[[1]]
  current_chunk <- ids[1]
  chunks <- vector("list")  # Predefine the list size
  
  for (i in seq_along(ids)[-1]) {  # Start loop from the second element
    if (nchar(paste(current_chunk, ids[i], sep = " + ")) <= max_length) {
      current_chunk <- paste(current_chunk, ids[i], sep = " + ")
    } else {
      chunks[[length(chunks) + 1]] <- current_chunk
      current_chunk <- ids[i]
    }
  }
  chunks[[length(chunks) + 1]] <- current_chunk  # Final chunk
  
  return(unlist(chunks))
}

# Check if overflow is needed, and if so, prepare TranscriptDetails
overflow_needed <- any(nchar(df165$Signature) > 30000)
if (overflow_needed) {
  addWorksheet(wb, "TranscriptDetails")
  writeData(wb, "TranscriptDetails", data.frame(Row_Index = "Row_Index", Nomenclature = "Nomenclature", Signature_Chunk = "Signature_Chunk"), startRow = 1)
}

# Prepare to write data in bulk by storing rows in lists before writing them all at once
main_data_hyperlinks <- vector("list", nrow(df165))
transcript_details_data <- list()
start_row <- 2  # Row in TranscriptDetails to start data (skip headers)

# Loop over each row in df165
for (i in 1:nrow(df165)) {
  if (nchar(df165$Signature[i]) > 30000) {
    chunks <- split_signature(df165$Signature[i])
    
    # Add row information for each chunk to transcript details in a bulk list
    for (chunk in chunks) {
      transcript_details_data[[length(transcript_details_data) + 1]] <- data.frame(
        Row_Index = i + 1,
        Nomenclature = df165$Nomenclature[i],
        Signature_Chunk = chunk
      )
    }
    
    # Prepare hyperlink text for MainData
    summary_text <- paste0("n=", df165$Elements[i], " transcript isoforms")
    main_data_hyperlinks[[i]] <- data.frame(
      row = i + 1,
      col = which(names(df165) == "Signature"),
      text = summary_text,
      formula = paste0("HYPERLINK(\"#TranscriptDetails!A", start_row, "\", \"See full list\")")
    )
    
    start_row <- start_row + length(chunks)
  } else {
    # If within limit, write Signature entry directly to MainData
    main_data_hyperlinks[[i]] <- data.frame(
      row = i + 1,
      col = which(names(df165) == "Signature"),
      text = df165$Signature[i],
      formula = NA  # No formula required for short entries
    )
  }
}

# Bulk write transcript details data to reduce function calls
if (overflow_needed && length(transcript_details_data) > 0) {
  transcript_details_df <- do.call(rbind, transcript_details_data)
  writeData(wb, "TranscriptDetails", transcript_details_df, startRow = 2, colNames = FALSE)
}

# Write hyperlinks and text data to MainData in one go
for (item in main_data_hyperlinks) {
  writeData(wb, "MainData", item$text, startRow = item$row, startCol = item$col)
  if (!is.na(item$formula)) {
    writeFormula(wb, "MainData", item$formula, startRow = item$row, startCol = item$col)
  }
}

# Save workbook
saveWorkbook(wb, "Transcript_S1e_Optimized_Supplementary_Data_with_Selective_Hyperlinks_v02.xlsx", overwrite = TRUE)

# Import saved workbook to df166 for verification
df166_v06 <- import("Transcript_S1e_Optimized_Supplementary_Data_with_Selective_Hyperlinks_v02.xlsx")

######
######
##### Optimized Signature Data rio::export with Reduced Loop Write Operations #####
######
######
######
# Load required packages
library(openxlsx)
# Assuming `import` is from `rio` package

# Load data
df165 <- Dataset_S1_Transcript


# Create a new workbook
wb <- createWorkbook()

# Add the main worksheet for primary data
addWorksheet(wb, "MainData")

# Write main data (excluding Signature) to MainData sheet in one go
writeData(wb, "MainData", df165[, -which(names(df165) == "Signature")])

# Define optimized function to split long Signature strings into chunks under Excel's character limit
split_signature <- function(signature_entry, max_length = 30000) {
  ids <- strsplit(signature_entry, " \\+ ")[[1]]
  current_chunk <- ids[1]
  chunks <- vector("list")  # Predefine the list size
  
  for (i in seq_along(ids)[-1]) {  # Start loop from the second element
    if (nchar(paste(current_chunk, ids[i], sep = " + ")) <= max_length) {
      current_chunk <- paste(current_chunk, ids[i], sep = " + ")
    } else {
      chunks[[length(chunks) + 1]] <- current_chunk
      current_chunk <- ids[i]
    }
  }
  chunks[[length(chunks) + 1]] <- current_chunk  # Final chunk
  
  return(unlist(chunks))
}

# Check if overflow is needed, and if so, prepare TranscriptDetails
overflow_needed <- any(nchar(df165$Signature) > 30000)
if (overflow_needed) {
  addWorksheet(wb, "TranscriptDetails")
  writeData(wb, "TranscriptDetails", data.frame(Row_Index = "Row_Index", Nomenclature = "Nomenclature", Signature_Chunk = "Signature_Chunk"), startRow = 1)
}

# Prepare to write data in bulk by storing rows in lists before writing them all at once
main_data_hyperlinks <- vector("list", nrow(df165))
transcript_details_data <- list()
start_row <- 2  # Row in TranscriptDetails to start data (skip headers)

# Loop over each row in df165
for (i in 1:nrow(df165)) {
  if (nchar(df165$Signature[i]) > 30000) {
    chunks <- split_signature(df165$Signature[i])
    
    # Add row information for each chunk to transcript details in a bulk list
    for (chunk in chunks) {
      transcript_details_data[[length(transcript_details_data) + 1]] <- data.frame(
        Row_Index = i + 1,
        Nomenclature = df165$Nomenclature[i],
        Signature_Chunk = chunk
      )
    }
    
    # Prepare hyperlink text for MainData
    summary_text <- paste0("n=", df165$Elements[i], " transcript isoforms")
    main_data_hyperlinks[[i]] <- data.frame(
      row = i + 1,
      col = which(names(df165) == "Signature"),
      text = summary_text,
      formula = paste0("HYPERLINK(\"#TranscriptDetails!A", start_row, "\", \"See full list\")")
    )
    
    start_row <- start_row + length(chunks)
  } else {
    # If within limit, write Signature entry directly to MainData
    main_data_hyperlinks[[i]] <- data.frame(
      row = i + 1,
      col = which(names(df165) == "Signature"),
      text = df165$Signature[i],
      formula = NA  # No formula required for short entries
    )
  }
}

# Bulk write transcript details data to reduce function calls
if (overflow_needed && length(transcript_details_data) > 0) {
  transcript_details_df <- do.call(rbind, transcript_details_data)
  writeData(wb, "TranscriptDetails", transcript_details_df, startRow = 2, colNames = FALSE)
}

# Write hyperlinks and text data to MainData in one go
for (item in main_data_hyperlinks) {
  writeData(wb, "MainData", item$text, startRow = item$row, startCol = item$col)
  if (!is.na(item$formula)) {
    writeFormula(wb, "MainData", item$formula, startRow = item$row, startCol = item$col)
  }
}

# Save workbook
saveWorkbook(wb, "Transcript_S1f_Optimized_Supplementary_Data_with_Selective_Hyperlinks_v02.xlsx", overwrite = TRUE)

# Import saved workbook to df166 for verification
df166_v07 <- import("Transcript_S1f_Optimized_Supplementary_Data_with_Selective_Hyperlinks_v02.xlsx")

######
######
######
###### Remaking MOST_CLINICALLY_meaningful_filtered_signatures_df1157
###### to replace Dataset S1V_NOVO sheet in Dataset S1 (02/23/2025)
###### 167 clinically meaningful signatures across five omic features: 
###### Transcript, mRNA, CNV, Methylation, and Mutation.
# None_effect value_signatures: Select all rows where the value in "HRC_series" is "1A2A3A4A"
filtered_rows_HRC_df1157 <- df1157 %>%
  filter(HRC_series == "1A2A3A4A") %>%  # Filter rows where HRC_series is "1A2A3A4A"
  filter(TNC == "0") %>%  # Exclude rows where TNC is "0"
  filter(SMC != "0") 

filtered_rows_SMC_df1157 <- df1157 %>%
  filter(SMC_series == "1A2A3A4A")%>%  # Filter rows where HRC_series is "1A2A3A4A"
  filter(TNC == "0") %>%  # Exclude rows where TNC is "0"
  filter(HRC != "0") 

#### No effect value signatures: OR
filtered_rows_HRC_df1157_SMC <- df1157 %>%
  filter(HRC_series == "1A2A3A4A" | SMC_series == "1A2A3A4A")

#### No effect value signatures: AND
filtered_rows_HRC_df1157_SMC <- df1157 %>%
  filter(HRC_series == "1A2A3A4A" & SMC_series == "1A2A3A4A") # Logical AND: Both conditions must be true

#### Select all rows that are meaningful by "HRC_series" and "SMC_series"
# Select all rows where the value in "HRC_series" an  ""SMC_series) is different that "1A2A3A4A" 
filtered_rows_meaningful_df1157 <- df1157 %>%
  filter(HRC_series != "1A2A3A4A", SMC_series != "1A2A3A4A") %>%  # Exclude rows where TNC is "0"
  filter(TNC != "0") 

# Select all rows that are meaningful by the value in "HRC_series", "SMC_series", "immune_classification", and "microenvironment_classification") is different that "1A2A3A4A"
meaningful_filtered_signatures_df1157 <- df1157 %>%
  filter(HRC_series != "1A2A3A4A", SMC_series != "1A2A3A4A", immune_classification != "NS", microenvironment_classification != "NS")  %>%  # Exclude rows where TNC is "0"
  filter(TNC != "0") 

# Select all rows that are meaningful by the value in "HRC_series", "SMC_series", "immune_classification", and "microenvironment_classification") is different that "1A2A3A4A"
meaningful_filtered_signatures_df1157 <- df1157 %>%
  filter(HRC_series != "1A2A3A4A", SMC_series != "1A2A3A4A", immune_classification != "NS", microenvironment_classification != "NS",
         Type_log_rank_DSS != "NS", Type_log_rank_DFI != "NS", Type_log_rank_PFI != "NS", Type_log_rank_OS != "NS",
         Type_Cox_DSS != "NS", Type_Cox_DFI != "NS", Type_Cox_PFI != "NS", Type_Cox_OS != "NS" ) %>%  # Exclude rows where TNC is "0"
  filter(TNC != "0") 

#####
#####
##### Remaking TOP MOST CLINICALLY MEANINGFUL SIGNATURES based on df1157
##### Table 6b. Top most clinically meaningful signatures 
##### 
##### 

# Select rows that are meaningful based on specified conditions
MOST_meaningful_filtered_signatures_df1157 <- df1157 %>%
  filter(
    # Exclude rows where specific series values are "1A2A3A4A"
    HRC_series != "1A2A3A4A",
    SMC_series != "1A2A3A4A",
    
    # Exclude rows where classifications are "NS"
    immune_classification != "NS",
    microenvironment_classification != "NS",
    
    # Exclude rows where specific log rank types are "NS"
    Type_log_rank_DSS != "NS",
    Type_log_rank_DFI != "NS",
    Type_log_rank_PFI != "NS",
    Type_log_rank_OS != "NS",
    
    # Exclude rows where specific Cox types are "NS"
    Type_Cox_DSS != "NS",
    Type_Cox_DFI != "NS",
    Type_Cox_PFI != "NS",
    Type_Cox_OS != "NS",
    
    # Exclude rows where TNC is "0"
    TNC != "0"
  ) %>%
  arrange(desc(Rank))  # Order by "Ranking" in descending order

df1157b <- MOST_meaningful_filtered_signatures_df1157

# Select rows where HRC_Rank and SMC_Rank are 11, and TMC_Rank and TIC_Rank are 7
top_MOST_ranked_signatures_df1157b <- df1157b %>%
  filter(SMC_ranking == 11, TMC_ranking == 7, TIC_ranking == 7) %>%
  arrange(desc(Rank))  # Ensure descending order

# Select specified column variables
filtered_columns_top_MOST_ranked_signatures_df1157b <- top_MOST_ranked_signatures_df1157b[, c(1, 2, 4, 20, 21, 36, 39, 42, 45, 47, 49, 51, 53, 54, 56)]

# Reorder specified columns and retain the remaining columns in their original order
reordered_filtered_columns_top_MOST_ranked_signatures_df1157b <- filtered_columns_top_MOST_ranked_signatures_df1157b %>%
  select(1, 3, 2, 5, 4, everything())

# Reorder the columns in reordered_filtered_columns_top_MOST_ranked_signatures_df1157b
reordered_filtered_columns_top_MOST_ranked_signatures_df1157b <- reordered_filtered_columns_top_MOST_ranked_signatures_df1157b %>%
  
  relocate(Type_log_rank_OS, .after = Type_log_rank_PFI) %>%
  relocate(Type_Cox_OS, .after = Type_Cox_PFI)

# Rename columns
reordered_filtered_columns_top_MOST_ranked_signatures_df1157b <- reordered_filtered_columns_top_MOST_ranked_signatures_df1157b %>%
  rename(
    TMC = microenvironment_classification,
    TIC = immune_classification
  )

# Print confirmation message
print("Columns renamed: 'microenvironment_classification' → 'TMC', 'immune_classification' → 'TIC'.")

# rio::export the updated table to Excel
write.xlsx(
  reordered_filtered_columns_top_MOST_ranked_signatures_df1157b, 
  "Table 6b_TOP_MOST_meaningful_filtered_signatures_df1157b_updated_reordered.xlsx", 
  rowNames = FALSE
)
# rio::export ordered and renamed dataframes
rio::export(MOST_meaningful_filtered_signatures_df1157, "Dataset_S1V_NOVO_MOST_meaningful_filtered_signatures_df1157b_updated.tsv") 

# rio::export the reordered dataframe
rio::export(reordered_filtered_columns_top_MOST_ranked_signatures_df1157b, "TOP_MOST_meaningful_filtered_signatures_df1157b_updated.tsv") 
rio::export(reordered_filtered_columns_top_MOST_ranked_signatures_df1157b, "TOP_MOST_meaningful_filtered_signatures_df1157b_updated.xlsx") 

##### 
##### 
##### Validating the number of gene symbols present and absent in signature database (dB)
##### Having had handled "PRKN" to "PARK2" now handle "MRE11A" to "MRE11"
#####
#####
# Step 1: Import `Target_genes.xlsx`
Target_genes <- import("Target_genes.csv")  # Adjust sheet if necessary

# Step 2: Preprocess `Gene_Symbols_signatures_df1157`
Gene_Symbols_signatures_df1157 <- import("Gene_Symbols_signatures_df1157.tsv")

Gene_Symbols_signatures_df1157 <- Gene_Symbols_signatures_df1157 %>%
  mutate(Gene_Symbol = trimws(Gene_Symbol)) %>%  # Trim whitespace
  mutate(Gene_Symbol = ifelse(Gene_Symbol == "MRE11A", "MRE11", Gene_Symbol)) %>%  # Replace PRKN
  distinct()  # Remove duplicate rows

# Print confirmation message
print("Preprocessed `Gene_Symbols_signatures_df1157`: Trimmed whitespace, replaced 'PRKN' with 'PARK2', and removed duplicates.")

# Step 3: Preprocess `Target_genes`
Target_genes <- Target_genes %>%
  mutate(Gene = trimws(Gene))  # Trim whitespace before matching

library(dplyr)

# Step 1: Identify duplicated values in "Gene"
duplicated_genes <- Target_genes %>%
  group_by(Gene) %>%
  filter(n() > 1) %>%
  ungroup()

# Step 2: Count total duplicated occurrences
duplicated_count <- duplicated_genes %>%
  count(Gene, name = "Duplicate_Count") %>%
  arrange(desc(Duplicate_Count))

# Display results
if (nrow(duplicated_count) > 0) {
  print("Duplicated 'Gene' values found:")
  print(duplicated_count)
} else {
  print("No duplicated 'Gene' values found in `Target_genes`.")
}

# Step 4: Add "Present in signature dB" variable after "Gene"
Target_genes <- Target_genes %>%
  mutate(`Present in Signature dB` = ifelse(Gene %in% Gene_Symbols_signatures_df1157$Gene_Symbol, "Yes", "No")) %>%
  relocate(`Present in Signature dB`, .after = Gene)

# Step 5: Count distribution of "Yes" and "No" values
signature_distribution_dB <- Target_genes %>%
  count(`Present in Signature dB`) %>%
  mutate(Percentage = round(n / sum(n) * 100, 2))  # Calculate percentage

# Print the distribution summary
print("Distribution of 'Present in signature dB' values:")
print(signature_distribution_dB)

# Step 6: Create a new dataframe with rows where "Present in signature dB" is "No"
Target_genes_absent_in_signatures <- Target_genes %>%
  filter(`Present in Signature dB` == "No")

# Rename columns in Target_genes
Target_genes <- Target_genes %>%
  rename(
    `Entrez ID` = uid,
    `Gene symbol` = Gene,
    `Present in Signature dB` = `Present in Signature dB`,
    `COSMIC Driver` = Driver,
    `COSMIC Tier` = TIER,
    `Category Search` = Category_Search,
    `RCD count` = Members_Pathway_Categorization,
    `RCD form` = Pathway_Categorization,
    `Description` = description,
    `Summary` = summary,
    `Other aliases` = otheraliases,
    `Chromosome` = chromosome,
    `Start` = chrstart,
    `End` = chrstop,
    `Map location` = maplocation
  )

# Print confirmation message
print("Column names in `Target_genes` have been successfully renamed.")

# Step 7: Save the filtered dataframe to a .tsv file
write_tsv(Target_genes_absent_in_signatures, "Target_genes_absent_in_signature_dB.tsv")

# Step 8: Save the filtered dataframe to a .tsv file
write_tsv(Target_genes, "Target_genes_final_dB.tsv")
 
rio::export(Target_genes, "Target_genes_final_dB.xlsx")

# Print confirmation message
print("Processed `Target_genes`: Added 'Present in Signature dB' column, computed distribution, and saved absent genes to 'Target_genes_absent_in_signatures_dB.tsv'.")

##### Issue saving as Excel because character lenght exceeding 30000 characters
##### Identifying variables whsie values execeed 32767 characters.
##### Example: MOST_meaningful_filtered_signatures_df1157

###
###
### Part A. Oversize values in df1157
### 
### 
# Define Excel's character limit
max_excel_characters <- 32767

# Identify oversized rows for each column
oversized_rows <- df1157 %>%
  mutate(across(where(is.character), ~ nchar(.))) %>%  # Compute character count for all text columns
  rowwise() %>%
  mutate(Row_Number = cur_group_id()) %>%  # Preserve row index
  ungroup() %>%
  pivot_longer(-Row_Number, names_to = "Column", values_to = "Char_Count") %>%
  filter(Char_Count > max_excel_characters)  # Only keep rows exceeding the limit

# Print oversized row report
if (nrow(oversized_rows) > 0) {
  print("⚠️ Warning: The following rows exceed the Excel character limit:")
  print(oversized_rows)
} else {
  print("✅ No rows exceed the Excel character limit.")
}

### Creating df1159 without the oversized variable "immune_score_details"
#### Note: given that the oversized values are all under the variable "immune_score_detail",
#### to split the df into omic-specific Dataset_S1_, the oversized column will be omitted thus creating df df1159 
#### df1159 also has a =n oversized row value under "Signature" which will nee to be handled using split chuncks with "wb <- createWorkbook()
# Create df1159 by excluding "immune_score_details"
df1159 <- df1157 %>% select(-immune_score_details)

# Print confirmation message
print("✅ Created df1159 with 'immune_score_details' removed.")

Dataset_S2 <- df1159

rio::export(df1159, "df1159.tsv")
rio::export(Dataset_S2, "Dataset_S2.tsv")

### 
### 
### Part B. Oversize values in "MOST_meaningful_filtered_signatures_df1157" 
###
###
# Define Excel's character limit
max_excel_characters <- 32767

# Identify oversized rows for each column
oversized_rows <- MOST_meaningful_filtered_signatures_df1157 %>%
  mutate(across(where(is.character), ~ nchar(.))) %>%  # Compute character count for all text columns
  rowwise() %>%
  mutate(Row_Number = cur_group_id()) %>%  # Preserve row index
  ungroup() %>%
  pivot_longer(-Row_Number, names_to = "Column", values_to = "Char_Count") %>%
  filter(Char_Count > max_excel_characters)  # Only keep rows exceeding the limit

# Print oversized row report
if (nrow(oversized_rows) > 0) {
  print("⚠️ Warning: The following rows exceed the Excel character limit:")
  print(oversized_rows)
} else {
  print("✅ No rows exceed the Excel character limit.")
}

##### Create MOST_meaningful_filtered_signatures_df1159 by excluding "immune_score_details"
##### Note: given that the oversized values are msotly under the variable ~="immune_score_detail",
#### to split the df into omic-specific Dataset_S1_, the over-sized column will be omitted thus creating df df1159 
#### df1157 has one oversized value under Signature which will need to hadled using wd  
MOST_meaningful_filtered_signatures_df1159 <- MOST_meaningful_filtered_signatures_df1157 %>% select(-immune_score_details)

# Print confirmation message
print("✅ Created MOST_meaningful_filtered_signatures_df1159 with 'immune_score_details' removed.")

rio::export(MOST_meaningful_filtered_signatures_df1159, "Dataset_S1V_NOVO_MOST_meaningful_filtered_signatures.xlsx")
rio::export(MOST_meaningful_filtered_signatures_df1159, "Dataset_S1V_NOVO_MOST_meaningful_filtered_signatures.tsv")

####
####
#### Making a consolidated dataframe with signatures from tables 1, 2, 3, 6
#### filtering d1157 for df1156 with n = 45 top signatures inf Tables 1, 2, 3 and 6 in the Manuscript
#### Making Dataset S1P_NOVO_df1157_filtered_druggable top signatures
#### Application CancerRCDShiny
#### Note>: Considering applying the login to the set of top clinical meaningful signatures
#### 
# Step 0: Filter df1157 based on matching "Nomenclature" values in df1156
df1157_filtered_druggable <- df1157 %>%
  filter(Nomenclature %in% df1156$Nomenclature) %>%
  arrange(desc(Rank))  # Order by "Rank" in descending order

# Print confirmation message
print("✅ Filtered and ordered df1157_filtered_druggable by 'Rank' in descending order.")

rio::export(df1157_filtered_druggable,"Dataset S1P_NOVO_df1157_filtered_druggable.xlsx")

####
####
#### Sub-setting df1157 by`Omic feature` values for final Dataset in tsv format
####
####

# Step 1: Get unique "Omic feature" values
unique_genotypes <- unique(df1157$`Omic feature`)

# Step 2: Loop through each unique "Omic feature" and save separately
for (genotype in unique_genotypes) {
  # Filter dataframe for the current "Omic feature"
  genotype_df <- df1157 %>% filter(`Omic feature` == genotype)
  
  # Create a sanitized file name (replace spaces/special characters with underscores)
  safe_genotype <- gsub("[^A-Za-z0-9_]", "_", genotype)
  file_name_tsv <- paste0("Dataset_S1_", safe_genotype, ".tsv")
  
  # Save the dataframe as .tsv using rio::export()
  rio::export(genotype_df, file_name_tsv, format = "tsv")
  
  # Print confirmation
  print(paste("File created:", file_name_tsv))
}

# Print completion message
print("✅ All files created successfully!")

####
#### 
#### Distribution of RCD forms in Target genes (Dataset S1L)
####
####
# Step 1: Separate the terms into individual rows using '/' as the delimiter
expanded_df <- Target_genes %>%
  separate_rows(`RCD form`, sep = "/")

# Step 2: Trim any leading or trailing whitespace from the terms
expanded_df <- expanded_df %>%
  mutate(`RCD form` = str_trim(`RCD form`))

# Step 3: Count the occurrences of each term
term_counts <- expanded_df %>%
  group_by(`RCD form`) %>%
  summarise(`Count term-based associations in target genes` = n()) %>%
  ungroup()

# Step 4: Calculate the total number of terms
total_terms <- sum(term_counts$`Count term-based associations in target genes`)

# Step 5: Calculate the percentage for each term
term_counts <- term_counts %>%
  mutate(Percentage = (`Count term-based associations in target genes` / total_terms) * 100)

# Step 6: Sort the dataframe in descending order based on the count
RCD_distribution_in_target_genes <- term_counts %>%
  arrange(desc(`Count term-based associations in target genes`))

rio::export(RCD_distribution_in_target_genes, "RCD_distribution_in_target_genes.xlsx")

####
####
####
#### HLA analysis (03/10/2025)
####
####
####
#### Filtering df1157 for any match (exact or partial) of "HLA" in "Signature" ####

df1157_filtered_HLA <- df1157 %>%
  filter(grepl("HLA", Signature, ignore.case = TRUE) & 
           Combined_Outcome != "Meaningless" & 
           Expression != "Unchanged")

# Print confirmation message
message("Filtered dataframe created: df1157_filtered_HLA")

rio::export(df1157_filtered_HLA, "df1157_filtered_HLA.xlsx")

####
#### 
#### 
######################
######################
# Save the entire workspace to a file
save.image(file = "Rconstructing_workspace.RData")
save.image(file = "D:/Pré-artigo 5-optosis model/Dataframes all metrics/Backup workspace/Rconstructing_workspace.RData")
#####################
#####################
####
####
####
####
####


#####
#####
#####
# Load the saved workspace into the current session
# load("Rconstructing_workspace.RData")
#####
#####
