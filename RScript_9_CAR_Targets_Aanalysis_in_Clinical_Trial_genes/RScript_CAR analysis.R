### CAR analysis from "Dataset_S2_full.tsv"
### PART A - Centered in gene symbol signatures. Transcript and protein signature will be analysis in parts B ad C

setwd("D:/Pré-artigo 5-optosis model/Dataframes all metrics/CAR Targets Analysis in Clinical Trial genes")
### Note: API Key = e89a9427ce52284b89e414389903bd387207 

library(rio)
library(dplyr)
library(stringr)

# Load the data
dfCAR_targets <- import("CAR Targets in Clinical Trials as Identified from Clinicaltrials_gov.xlsx")
Gene_Symbol <- dfCAR_targets$Gene_Symbol

df300 <- dfCAR_targets
df301 <- df300$Gene_Symbol
df302 <-  import("Dataset_S2_full.tsv")
df303 <-import("Target_genes_final_dB.tsv")

####
####
####
#### CAR analysis
#### 
#### 
#### 

setwd("D:/Pré-artigo 5-optosis model/Dataframes all metrics/CAR Targets Analysis in Clinical Trial genes")

# Load the data
dfCAR_targets <- import("CAR Targets in Clinical Trials as Identified from Clinicaltrials_gov.xlsx")
Gene_Symbol <- dfCAR_targets$Gene_Symbol

# Step 1: Process the "Signature" variable in df302
df302 <- df302 %>%
  mutate(Processed_Signatures = Signature) %>%
  mutate(Processed_Signatures = str_replace_all(Processed_Signatures, pattern = "\\s*\\+\\s*", replacement = ",")) %>% # Replace "+" with ","
  mutate(Processed_Signatures = str_remove_all(Processed_Signatures, pattern = "[\\(\\)]")) %>% # Remove parentheses
  mutate(Processed_Signatures = str_remove_all(Processed_Signatures, pattern = "`")) # Remove backticks

# Step 2: Create a mapping function
check_gene_match <- function(gene_symbol, signature) {
  # Check if the gene_symbol exists in the processed signature
  any(str_detect(signature, regex(paste0("\\b", gene_symbol, "\\b"), ignore_case = TRUE)))
}

# Step 3: Create "CAR target present" variable
df302 <- df302 %>%
  rowwise() %>%
  mutate(`CAR target present` = any(sapply(df300$Gene_Symbol, check_gene_match, signature = Processed_Signatures))) %>%
  ungroup()

# Step 4: Create "CAR Matched Gene" variable
df302 <- df302 %>%
  rowwise() %>%
  mutate(`CAR Matched Gene` = paste(df300$Gene_Symbol[sapply(df300$Gene_Symbol, check_gene_match, signature = Processed_Signatures)], collapse = ", ")) %>%
  ungroup()

# Step 5: Generate filtered dataframe
Dataset_S1W <- df302 %>%
  filter(`CAR target present` == TRUE)

# Step 6: Relocate variables in Dataset_S1W
Dataset_S1W <- Dataset_S1W %>%
  relocate(Processed_Signatures, `CAR target present`, `CAR Matched Gene`, .after = Signature)

# Step 7: Extract and list unique CAR matched gene members
unique_car_gene_list <- Dataset_S1W %>%
  pull(`CAR Matched Gene`) %>% # Extract the "CAR Matched Gene" column
  str_split(", ") %>% # Split strings by ", " to separate genes
  unlist() %>% # Flatten the list into a single vector
  unique() # Keep only unique genes

# Print the list of unique CAR matched genes
print("Unique CAR matched gene members:")
print(unique_car_gene_list)

# Step 8: Count unique values in "CAR Matched Gene"
  unique_car_genes <- Dataset_S1W %>%
  pull(`CAR Matched Gene`) %>% # Extract the "CAR Matched Gene" column
  str_split(", ") %>% # Split strings by ", " to separate genes
  unlist() %>% # Flatten the list into a single vector
  unique() %>% # Keep only unique genes
  length() # Count the number of unique genes

# Print the count of unique CAR matched genes
print(paste("Number of unique CAR matched genes:", unique_car_genes))
# Count the number of unique CAR matched genes
unique_car_gene_count <- length(unique_car_gene_list)
print(paste("Number of unique CAR matched genes:", unique_car_gene_count))
# > unique_car_genes
# [1] 27

# Print the count of unique CAR matched genes
print(paste("Number of unique CAR matched genes:", unique_car_genes))

# Step 9: Identify missing "CAR Genes" in the unique_car_gene_list
missing_car_genes <- setdiff(df300$Gene_Symbol, unique_car_gene_list)

# Print the missing CAR genes
print("CAR Genes missing in the unique_car_gene_list:")
print(missing_car_genes)
# > print(missing_car_genes)
# [1] "B4GALNT1" "FOLH1" 
# 
# Check the maximum string length for each column
string_lengths <- Dataset_S1W %>%
  summarise(across(everything(), ~ max(nchar(as.character(.))))) %>%
  pivot_longer(everything(), names_to = "Column", values_to = "Max_Length")# 

rio::export(Dataset_S1W, "Dataset_S1W_final.xlsx" )

# Step 10. Adding the "CAR Matched gene" variable to df303 (Target_genes_final_dB)
# Step 10: Create the "CAR Matched gene" variable in df303
# This out-out is the new final Dataset_S1B for the supplementary material
df303 <- df303 %>%
  mutate(`CAR Matched gene` = ifelse(`Gene symbol` %in% unique_car_gene_list, "Yes", "No"),
         .after = `Present in Signature dB`) # Place the new column immediately after "Present in Signature dB"

rio::export(df303,"Dataset S1B_final.tsv" )
rio::export(df303,"Dataset S1B_final.xlsx" )

####### Part B - CAR Analysis from "Dataset_S4_Transcript.tsv"
####### 
####### 
####### CAR analysis for transcripts to gene symbols in "Gene Symbol members"
setwd("D:/Pré-artigo 5-optosis model/Dataframes all metrics/CAR Targets Analysis in Clinical Trial genes")
### Note: API Key = e89a9427ce52284b89e414389903bd387207

library(rio)
library(dplyr)
library(stringr)
library(tidyr)

# Load the data
dfCAR_targets <- import("CAR Targets in Clinical Trials as Identified from Clinicaltrials_gov.xlsx")
Gene_Symbol <- dfCAR_targets$Gene_Symbol

df300 <- dfCAR_targets
df301 <- df300$Gene_Symbol
df302 <-  import("Dataset_S2_full.tsv")
df303 <-import("Target_genes_final_dB.tsv")
df304 <- import("Dataset_S4_Transcript.tsv")

# Step 1: Process the "Gene Symbol members" variable in df304
df304 <- df304 %>%
  mutate(Processed_Gene_Symbol_members = `Gene Symbol members`) %>%
  mutate(Processed_Gene_Symbol_members = str_replace_all(Processed_Gene_Symbol_members, pattern = "\\s*\\+\\s*", replacement = ",")) %>% # Replace "+" with ","
  mutate(Processed_Gene_Symbol_members = str_remove_all(Processed_Gene_Symbol_members, pattern = "\\([^()]+\\)")) %>% # Remove parentheses and their contents
  mutate(Processed_Gene_Symbol_members = str_remove_all(Processed_Gene_Symbol_members, pattern = "[\\(\\)]")) %>% # Remove any remaining parentheses
  mutate(Processed_Gene_Symbol_members = str_remove_all(Processed_Gene_Symbol_members, pattern = "`")) # Remove backticks

# Step 2: Create a mapping function
check_gene_match <- function(gene_symbol, signature) {
  # Check if the gene_symbol exists in the processed signature
  any(str_detect(signature, regex(paste0("\\b", gene_symbol, "\\b"), ignore_case = TRUE)))
}

# Step 3: Create "CAR target present" variable
df304 <- df304 %>%
  rowwise() %>%
  mutate(`CAR target present` = any(sapply(df300$Gene_Symbol, check_gene_match, signature = Processed_Gene_Symbol_members))) %>%
  ungroup()

# Step 4: Create "CAR Matched Gene" variable
df304 <- df304 %>%
  rowwise() %>%
  mutate(`CAR Matched Gene` = paste(df300$Gene_Symbol[sapply(df300$Gene_Symbol, check_gene_match, signature = Processed_Gene_Symbol_members)], collapse = ", ")) %>%
  ungroup()

# Step 5: Generate filtered dataframe
Dataset_S1W_transcript <- df304 %>%
  filter(`CAR target present` == TRUE)

# Step 6: Relocate variables in Dataset_S1W_transcript
Dataset_S1W_transcript <- Dataset_S1W_transcript %>%
  relocate(Processed_Gene_Symbol_members, `CAR target present`, `CAR Matched Gene`, .after = `Gene Symbol members`)

# Step 7: Extract and list unique CAR matched gene members
unique_car_gene_list <- Dataset_S1W_transcript %>%
  pull(`CAR Matched Gene`) %>% # Extract the "CAR Matched Gene" column
  str_split(", ") %>% # Split strings by ", " to separate genes
  unlist() %>% # Flatten the list into a single vector
  unique() # Keep only unique genes

# Print the list of unique CAR matched genes
print("Unique CAR matched gene members:")
print(unique_car_gene_list)

# Step 8: Count unique values in "CAR Matched Gene"
unique_car_genes <- Dataset_S1W_transcript %>%
  pull(`CAR Matched Gene`) %>% # Extract the "CAR Matched Gene" column
  str_split(", ") %>% # Split strings by ", " to separate genes
  unlist() %>% # Flatten the list into a single vector
  unique() %>% # Keep only unique genes
  length() # Count the number of unique genes

# Print the count of unique CAR matched genes
print(paste("Number of unique CAR matched genes:", unique_car_genes))

# Step 9: Identify missing "CAR Genes" in the unique_car_gene_list
missing_car_genes <- setdiff(df300$Gene_Symbol, unique_car_gene_list)

# Print the missing CAR genes
print("CAR Genes missing in the unique_car_gene_list:")
print(missing_car_genes)
# > print(missing_car_genes)
# [1] "B4GALNT1" "FOLH1"  

# Step 10. check for cell character limit of 32,767 characters
# Check the maximum string length for each column
string_lengths <- Dataset_S1W_transcript %>%
  summarise(across(everything(), ~ max(nchar(as.character(.))))) %>%
  pivot_longer(everything(), names_to = "Column", values_to = "Max_Length")

# Print the results
print(string_lengths)
# Signature and mmune_score_details exceeds the cell character limit!!!! (Threfore, df can not be saved as xlsx)

# Export the final dataframes
rio::export(Dataset_S1W_transcript, "Dataset_S1W_transcript_final.tsv")

######
######
######
######
### PART C - CAR analysis for proteins to gene symbols in "Gene Symbol members"
setwd("D:/Pré-artigo 5-optosis model/Dataframes all metrics/CAR Targets Analysis in Clinical Trial genes")

library(rio)
library(dplyr)
library(stringr)
library(tidyr)

# Load the data
dfCAR_targets <- import("CAR Targets in Clinical Trials as Identified from Clinicaltrials_gov.xlsx")
Gene_Symbol <- dfCAR_targets$Gene_Symbol

df300 <- dfCAR_targets
df301 <- df300$Gene_Symbol
df302 <-import("Dataset_S2_full.tsv")
df303 <-import("Target_genes_final_dB.tsv")
df304 <- import("Dataset_S4_Transcript.tsv")
df305 <- import("Dataset_S7_Protein.tsv")

# Step 1: Process the "Gene Symbol members" variable in df305
df305 <- df305 %>%
  mutate(Processed_Gene_Symbol_members = `Gene Symbol members`) %>%
  mutate(Processed_Gene_Symbol_members = str_replace_all(Processed_Gene_Symbol_members, pattern = "\\s*\\+\\s*", replacement = ",")) %>% # Replace "+" with ","
  mutate(Processed_Gene_Symbol_members = str_remove_all(Processed_Gene_Symbol_members, pattern = "\\([^()]+\\)")) %>% # Remove parentheses and their contents
  mutate(Processed_Gene_Symbol_members = str_remove_all(Processed_Gene_Symbol_members, pattern = "[\\(\\)]")) %>% # Remove any remaining parentheses
  mutate(Processed_Gene_Symbol_members = str_remove_all(Processed_Gene_Symbol_members, pattern = "`")) # Remove backticks

# Step 2: Create a mapping function
check_gene_match <- function(gene_symbol, signature) {
  # Check if the gene_symbol exists in the processed signature
  any(str_detect(signature, regex(paste0("\\b", gene_symbol, "\\b"), ignore_case = TRUE)))
}

# Step 3: Create "CAR target present" variable
df305 <- df305 %>%
  rowwise() %>%
  mutate(`CAR target present` = any(sapply(df300$Gene_Symbol, check_gene_match, signature = Processed_Gene_Symbol_members))) %>%
  ungroup()

# Step 4: Create "CAR Matched Gene" variable
df305 <- df305 %>%
  rowwise() %>%
  mutate(`CAR Matched Gene` = paste(df300$Gene_Symbol[sapply(df300$Gene_Symbol, check_gene_match, signature = Processed_Gene_Symbol_members)], collapse = ", ")) %>%
  ungroup()

# Step 5: Generate filtered dataframe
Dataset_S1W_protein <- df305 %>%
  filter(`CAR target present` == TRUE)

# Step 6: Relocate variables in Dataset_S1W_protein
Dataset_S1W_protein <- Dataset_S1W_protein %>%
  relocate(Processed_Gene_Symbol_members, `CAR target present`, `CAR Matched Gene`, .after = `Gene Symbol members`)

# Step 7: Extract and list unique CAR matched gene members
unique_car_gene_list <- Dataset_S1W_protein %>%
  pull(`CAR Matched Gene`) %>% # Extract the "CAR Matched Gene" column
  str_split(", ") %>% # Split strings by ", " to separate genes
  unlist() %>% # Flatten the list into a single vector
  unique() # Keep only unique genes

# Print the list of unique CAR matched genes
print("Unique CAR matched gene members:")
print(unique_car_gene_list)

# Step 8: Count unique values in "CAR Matched Gene"
unique_car_genes <- Dataset_S1W_protein %>%
  pull(`CAR Matched Gene`) %>% # Extract the "CAR Matched Gene" column
  str_split(", ") %>% # Split strings by ", " to separate genes
  unlist() %>% # Flatten the list into a single vector
  unique() %>% # Keep only unique genes
  length() # Count the number of unique genes

# Print the count of unique CAR matched genes
print(paste("Number of unique CAR matched genes:", unique_car_genes))

# Step 9: Identify missing "CAR Genes" in the unique_car_gene_list
missing_car_genes <- setdiff(df300$Gene_Symbol, unique_car_gene_list)

# Print the missing CAR genes
print("CAR Genes missing in the unique_car_gene_list:")
print(missing_car_genes)
# > print(missing_car_genes)
# [1] "B4GALNT1" "FOLH1"  

# Step 10. check for cell character limit of 32,767 characters
# Check the maximum string length for each column
string_lengths <- Dataset_S1W_protein %>%
  summarise(across(everything(), ~ max(nchar(as.character(.))))) %>%
  pivot_longer(everything(), names_to = "Column", values_to = "Max_Length")

# Print the results
print(string_lengths)

# Export the final dataframes
rio::export(Dataset_S1W_protein, "Dataset_S1W_protein_final.")

#######
#######
####### 
####### Binding Dataset_S1W_transcript and Dataset_S1W_protein_final
#######
#######
#######
# Load necessary libraries
library(dplyr)

# Step 1: Check if variable names and order are identical
identical_names <- identical(names(Dataset_S1W_protein), names(Dataset_S1W_transcript))

# Step 2: Identify differences if they exist
if (!identical_names) {
  # Find differences in variable names
  diff_names <- setdiff(names(Dataset_S1W_protein), names(Dataset_S1W_transcript))
  
  # Find differences in variable order
  order_mismatch <- which(names(Dataset_S1W_protein) != names(Dataset_S1W_transcript))
  
  # Print results
  if (length(diff_names) > 0) {
    cat("Differences in variable names:\n")
    print(diff_names)
  } else {
    cat("Variable names are identical.\n")
  }
  
  if (length(order_mismatch) > 0) {
    cat("Differences in variable order at positions:\n")
    print(order_mismatch)
  } else {
    cat("Variable order is identical.\n")
  }
} else {
  cat("Variable names and order are identical.\n")
}

# Step 3: Check if variable names and order are identical
identical_names <- identical(names(Dataset_S1W_protein), names(Dataset_S1W_transcript))

# Step 4: If identical, combine the dataframes
if (identical_names) {
  cat("Variable names and order are identical. Combining dataframes...\n")
  Dataset_S1W_transcript_protein <- rbind(Dataset_S1W_protein, Dataset_S1W_transcript)
  cat("Dataframes successfully combined into 'Dataset_S1W_transcript_protein'.\n")
} else {
  # Identify differences if they exist
  diff_names <- setdiff(names(Dataset_S1W_protein), names(Dataset_S1W_transcript))
  order_mismatch <- which(names(Dataset_S1W_protein) != names(Dataset_S1W_transcript))
  
  # Print results
  if (length(diff_names) > 0) {
    cat("Differences in variable names:\n")
    print(diff_names)
  } else {
    cat("Variable names are identical.\n")
  }
  
  if (length(order_mismatch) > 0) {
    cat("Differences in variable order at positions:\n")
    print(order_mismatch)
  } else {
    cat("Variable order is identical.\n")
  }
  
  cat("Cannot combine dataframes due to differences in variable names or order.\n")
}

rio::export(Dataset_S1W_transcript_protein, "Dataset_S1W_transcript_protein.tsv")


###
###
###
# Step 1: Rename the variable
Dataset_S1W_transcript_protein <- Dataset_S1W_transcript_protein %>%
  rename(Processed_Signatures = Processed_Gene_Symbol_members)

# Step 2: Export the updated dataframe (optional)
rio::export(Dataset_S1W_transcript_protein, "Dataset_S1W_transcript_protein_renamed.tsv")

## Making names(Dataset_S1W) identical to names in Dataset_S1W_transcript_protein for rbinding
# Step 1: Duplicate and rename the variable
Dataset_S1W <- Dataset_S1W %>%
  mutate(`Gene Symbol members` = Processed_Signatures, .after = Signature)

# Step 1: Check if variable names and order are identical
identical_names <- identical(names(Dataset_S1W_transcript_protein), names(Dataset_S1W))

# Step 2: Identify differences if they exist
if (!identical_names) {
  # Find differences in variable names
  diff_names <- setdiff(names(Dataset_S1W_transcript_protein), names(Dataset_S1W))
  
  # Find differences in variable order
  order_mismatch <- which(names(Dataset_S1W_transcript_protein) != names(Dataset_S1W))
  
  # Print results
  if (length(diff_names) > 0) {
    cat("Differences in variable names:\n")
    print(diff_names)
  } else {
    cat("Variable names are identical.\n")
  }
  
  if (length(order_mismatch) > 0) {
    cat("Differences in variable order at positions:\n")
    print(order_mismatch)
  } else {
    cat("Variable order is identical.\n")
  }
} else {
  cat("Variable names and order are identical.\n")
}





# Step 2: Identify differences if they exist
if (!identical_names) {
  # Find differences in variable names
  diff_names <- setdiff(names(Dataset_S1W), names(Dataset_S1W_transcript_protein))
  
  # Find differences in variable order
  order_mismatch <- which(names(Dataset_S1W) != names(Dataset_S1W_transcript_protein))
  
  # Print results
  if (length(diff_names) > 0) {
    cat("Differences in variable names:\n")
    print(diff_names)
  } else {
    cat("Variable names are identical.\n")
  }
  
  if (length(order_mismatch) > 0) {
    cat("Differences in variable order at positions:\n")
    print(order_mismatch)
  } else {
    cat("Variable order is identical.\n")
  }
} else {
  cat("Variable names and order are identical.\n")
}

# Step 3: Check if variable names and order are identical
identical_names <- identical(names(Dataset_S1W), names(Dataset_S1W_transcript_protein))

# Step 4: If identical, rbinf combine the dataframes
if (identical_names) {
  cat("Variable names and order are identical. Combining dataframes...\n")
  Dataset_S1W_mRNA_transcript_protein <- rbind(Dataset_S1W, Dataset_S1W_transcript_protein)
  cat("Dataframes successfully combined into 'Dataset_S1W_mRNA_transcript_protein'.\n")
} else {
  # Identify differences if they exist
  diff_names <- setdiff(names(Dataset_S1W), names(Dataset_S1W_transcript_protein))
  order_mismatch <- which(names(Dataset_S1W) != names(Dataset_S1W_transcript_protein))
  
  # Print results
  if (length(diff_names) > 0) {
    cat("Differences in variable names:\n")
    print(diff_names)
  } else {
    cat("Variable names are identical.\n")
  }
  
  if (length(order_mismatch) > 0) {
    cat("Differences in variable order at positions:\n")
    print(order_mismatch)
  } else {
    cat("Variable order is identical.\n")
  }
  
  cat("Cannot combine dataframes due to differences in variable names or order.\n")
}


##### Idnetifying duplicate vakues under the varibake "Nomemclature"
##### 
# Step 1: Identify duplicated values in the "Nomenclature" variable
duplicated_values <- Dataset_S1W_mRNA_transcript_protein %>%
  group_by(Nomenclature) %>%
  filter(n() > 1) %>%
  summarise(Duplicated_Count = n(), .groups = 'drop')

# Step 2: Report the duplicated values
if (nrow(duplicated_values) > 0) {
  cat("Duplicated values in 'Nomenclature':\n")
  print(duplicated_values)
} else {
  cat("No duplicated values found in 'Nomenclature'.\n")
}

# Step 1: Identify fully duplicated rows
fully_duplicated_rows <- Dataset_S1W_mRNA_transcript_protein %>%
  group_by(across(everything())) %>% # Group by all columns
  filter(n() > 1) %>% # Keep only fully duplicated rows
  ungroup()

# Step 2: Exclude one duplicated row in each instance
Dataset_S1W_mRNA_transcript_protein_cleaned <- Dataset_S1W_mRNA_transcript_protein %>%
  distinct(across(everything()), .keep_all = TRUE) # Keep only unique rows

# Step 3: Report the results
if (nrow(fully_duplicated_rows) > 0) {
  cat("Fully duplicated rows found. Excluding one duplicated row in each instance...\n")
  cat("Number of fully duplicated rows removed:", nrow(fully_duplicated_rows) - nrow(Dataset_S1W_mRNA_transcript_protein_cleaned), "\n")
} else {
  cat("No fully duplicated rows found.\n")
}

# Step 4: View the cleaned dataframe
View(Dataset_S1W_mRNA_transcript_protein_cleaned)

rio::export(Dataset_S1W_mRNA_transcript_protein_cleaned, "Dataset_S1W_consolidated_CAR_final.tsv")

rio::export(Dataset_S1W_mRNA_transcript_protein_cleaned, "Dataset_S10_CAR_NEW.tsv")

#######
#######
#######
#######
#######
# Save the entire workspace to a file
save.image(file = "my_workspace.RData")
#######
#######
#######
#######
#######

###
###
###
###
###
# Load the workspace from a file
load("my_workspace.RData")
###
###
###
###
