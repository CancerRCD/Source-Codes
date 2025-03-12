# Counting effective transcript per gene locus, converting transcript signatures in genes symbols with frequency
# and choosing cases for figures and manuscript

setwd("D:/Pré-artigo 5-optosis model/Dataframes all metrics_good/Transcript count Emanuell")

# Load necessary libraries
library(dplyr)
library(stringr)
library(tidyr)
library(networkD3)
library(htmlwidgets)
library(webshot)
library(magick)
library(purrr)

# Define the list of required packages
required_packages <- c("dplyr", "tidyr", "rio", "stringr", "readxl", 
                       "networkD3", "webshot", "webshot", "htmlwidgets", "magick", "pdftools")

# Step 1. Creating  function to conver transcript ID  to gene symbols
# Function to check, install, and load missing packages
install_and_load_packages <- function(packages) {
  # Identify packages that are not installed
  missing_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  
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

# Carregar os dados
Dataset_S1_Transcript <- import("Dataset_S1_Transcript.tsv")
gene_info_transcritos <- import("gene_info_trancritos_excel2.xlsx")
df160 <- Dataset_S1_Transcript

# Remover espaços extras nos IDs das assinatiras e dos transcritos
Dataset_S1_Transcript$Signature <- trimws(Dataset_S1_Transcript$Signature)
gene_info_transcritos$Transcript_ID <- trimws(gene_info_transcritos$Transcript_ID)

# Função para mapear transcritos a genes e contar transcritos únicos
map_transcripts_to_genes <- function(transcript_string, gene_info) {
  # Exibir a assinatura que está sendo processada
  cat("Processando assinatura:", transcript_string, "\n")
  
  # Separar os transcritos
  transcripts <- unlist(strsplit(gsub("[()]", "", transcript_string), " \\+ "))
  
  # Mapear cada transcrito ao gene correspondente
  gene_list <- sapply(transcripts, function(transcript) {
    gene <- gene_info %>% 
      filter(Transcript_ID == transcript) %>% 
      pull(Display_Name)
    
    # Verificar se o transcrito foi encontrado
    if(length(gene) == 0) {
      cat("Transcrito não encontrado:", transcript, "\n")
      return("N/A")
    } else {
      cat("Transcrito encontrado:", transcript, "->", gene, "\n")
      return(gene)
    }
  })
  
  # Contar os genes e formatar a saída
  gene_table <- table(gene_list)
  gene_counts <- sapply(names(gene_table), function(gene) {
    total_transcripts <- nrow(gene_info %>% filter(Display_Name == gene))
    paste0("(", gene, "(", gene_table[gene], "/", total_transcripts, ")", ")")
  })
  
  # Unir os genes de volta em uma string
  paste(gene_counts, collapse = " + ")
}

# Aplicar a função a cada linha da tabela de assinaturas
Dataset_S1_Transcript$member_gene <- sapply(Dataset_S1_Transcript$Signature, 
                                           map_transcripts_to_genes, 
                                           gene_info = gene_info_transcritos)

# Selecting rows where "member_gene" contains "(NA)" as a separate entity or is explicitly "NA"
df_filtered_NA <- Dataset_S1_Transcript %>%
  filter(is.na(member_gene) | grepl("\\(NA\\)", member_gene))

# Print the resulting dataframe
print(df_filtered_NA)

# Step 2: Add backticks only to terms with hyphens that don't already have them
Dataset_S1_Transcript$member_gene <- gsub("(?<!`)\\b([[:alnum:]]+-[[:alnum:]]+)\\b(?!`)", "`\\1`", Dataset_S1_Transcript$member_gene, perl = TRUE)

# Step 3: Remove numeric details inside parentheses, but keep the backtick-flanked terms intact
extracted_terms <- gsub("\\([0-9]+/[0-9]+\\)", "", Dataset_S1_Transcript$member_gene)

# Step 4: Remove extra spaces and split the terms into individual components
extracted_terms <- trimws(unlist(strsplit(extracted_terms, " \\+ ")))

# Rename 'member_gene' and relocate it immediately after 'Signature'
Dataset_S1_Transcript <- Dataset_S1_Transcript %>%
  rename(`Gene Symbol members` = member_gene) %>%
  relocate(`Gene Symbol members`, .after = Signature)



# Rename columns in Dataset_S1_Transcript
Dataset_S4_Transcript <- Dataset_S1_Transcript %>%
  rename(
    "Rank" = "Ranking",
    "Elements" = "Count_source",
    "RCD form" = "RCD_types",
    "Omic feature" = "Genotype"
  )%>%
  relocate(`Omic feature`, .after = Elements)

# Continue with the rest of your code
# (e.g., loading datasets, filtering, creating Sankey diagrams, etc.)

# Salvar o resultado em um novo arquivo 
rio::export(Dataset_S1_Transcript, "Dataset_S1_Transcript_with_counted_genes.tsv")
Dataset_S1_Transcript <- import("Dataset_S1_Transcript_with_counted_genes.tsv")

# Salvar o resultado em um novo arquivo _ Note: MUST use this df for the manuscript files as Dataset for ttancript
rio::export(Dataset_S4_Transcript, "Dataset_S4_Transcript.tsv")

# Step 5: Define a function to extract genes and their n/d values, handling backticks
extract_genes_and_ratios <- function(gene_string) {
  # Use regular expression to extract gene symbols (with or without backticks) and n/d pairs
  matches <- str_extract_all(gene_string, "(`?\\w+-?\\w+`?)\\((\\d+)/(\\d+)\\)")[[1]]
  
  # Create a dataframe with the extracted information
  df <- data.frame(
    gene = sapply(matches, function(x) str_match(x, "`?(\\w+-?\\w+)?`?")[,2]),  # Extract gene symbols, ignoring backticks
    n = as.numeric(sapply(matches, function(x) str_match(x, "\\((\\d+)/")[,2])),  # Extract n values
    d = as.numeric(sapply(matches, function(x) str_match(x, "/(\\d+)\\)")[,2]))   # Extract d values
  )
  
  # Calculate the n/d ratio
  df$ratio <- df$n / df$d
  return(df)
}

# Step 6: Apply the extraction function to the `Gene Symbol members` column and bind results
df_extracted <- do.call(rbind, lapply(Dataset_S1_Transcript$`Gene Symbol members`, extract_genes_and_ratios))

# Apply trimws() to all columns of the dataframe to remove leading and trailing spaces
df_extracted <- df_extracted %>% 
  mutate(across(everything(), ~ trimws(.)))

# Remove duplicate rows based on the "gene" column and then sort alphabetically by "gene"
df_extracted_no_duplicates <- df_extracted %>%
  distinct(gene, .keep_all = TRUE) %>%
  arrange(gene)

# Step 7: Rename the 'gene' column to 'Genes'
colnames(df_extracted_no_duplicates)[colnames(df_extracted_no_duplicates) == "gene"] <- "Genes"
df_extracted$n <- as.numeric(df_extracted$n)
df_extracted$d <- as.numeric(df_extracted$d)

# Step 8: Order df_extracted by descending values in "n" and then by "d"
df_extracted <- df_extracted %>%
  arrange(desc(n), desc(d))

# Step 9: Goal 1 - Identify genes where n/d = 1. Here we are looking for gebes with the highest number of transcript equal the total transcripts in a given signature
perfect_matches <- df_extracted %>% filter(ratio == 1)

# Step 10: Ensure the "n" variable is numeric
perfect_matches$n <- as.numeric(perfect_matches$n)

# Step 11: Order df_extracted by descending values in the "n" variable
perfect_matches <- perfect_matches %>%
  arrange(desc(n))

# Step 12: Goal 2 - Identify genes with the highest n/d values, including not perfect match
# Sorting by ratio descending, and selecting the highest n/d ratio
highest_nd_genes <- df_extracted %>% arrange(desc(ratio))

# Step 13: Create a dataframe for imperfect matches (n/d != 1)
imperfect_matches <- df_extracted %>% filter(ratio != 1)

# Step 13a: Create a dataframe for perfect matches (n/d == 1)
perfect_matches <- df_extracted %>% filter(ratio == 1)

######### Common transcript genes
df161 <- df_extracted_no_duplicates

# Step 14: Ensure the "n" variable is numeric
df161$n <- as.numeric(df161$n)

# Step 14a: Ensure the "d" variable is numeric
df161$d <- as.numeric(df161$d)

# Trim leading and trailing spaces
df161$Genes <- trimws(df161$Genes)

target_genes <- import("Target_genes.csv")

target_genes$Gene <- trimws(target_genes$Gene)

# Create a list "transcript_genes_6" with the values under the "Genes" column
transcript_genes_6 <- as.list(df161$Genes)

# Check which values from "transcript_genes_6" are present in the "Genes" column of df "Effective_target_genes"
common_genes_6 <- transcript_genes_6[transcript_genes_6 %in% target_genes$Gene]

# Identify genes that are NOT common
not_common_genes_6 <- transcript_genes_6[!transcript_genes_6 %in% target_genes$Gene]

# Calculate the percentage of common values
total_genes_6 <- length(transcript_genes_6)  # Total number of values in transcript_genes_6
common_genes_count_6 <- length(common_genes_6)  # Number of common values

# Estimate the percentage of common values
percentage_common_6 <- (common_genes_count_6 / total_genes_6) * 100

# Print the percentage
print(paste("Percentage of common genes:", round(percentage_common_6, 2), "%"))

# Print the genes that are NOT common
print("Genes that are NOT common:")
print(not_common_genes_6)
#[1] "PRKN"

#### Step 15. Common genes
df161 <- df_extracted_no_duplicates
# Trim leading and trailing spaces in the 'Gene_symbol_EMA' column
df161$Genes <- trimws(df161$Genes)

target_genes <- import("Target_genes.csv")

# Create a list "transcript_genes_6" with the values under the "Genes" column
transcript_genes_6 <- as.list(df161$Genes)

# Check which values from "transcript_genes_6" are present in the "Genes" column of df "Effective_target_genes"
common_genes_6 <- transcript_genes_6[transcript_genes_6 %in% target_genes$Gene]

# Calculate the percentage of common values
total_genes_6 <- length(transcript_genes_6)  # Total number of values in transcript_genes_6
common_genes_count_6 <- length(common_genes_6)  # Number of common values

# Estimate the percentage of common values
percentage_common_6 <- (common_genes_count_6 / total_genes_6) * 100

# Print the percentage
percentage_common_6

#### Step 16. Uncommon genes
# Trim leading and trailing spaces in the 'Genes' column
df161$Genes <- trimws(df161$Genes)

# Create a list "transcript_genes_7" with the values under the "Genes" column
transcript_genes_7 <- as.list(df161$Genes)

# Identify the uncommon gene symbols (those in transcript_genes_7 but not in target_genes)
uncommon_genes_7 <- transcript_genes_7[!(transcript_genes_7 %in% target_genes$Gene)]

Genes_transcripts <- unique(gene_info_transcritos$Display_Name)

#####
#####
##### Step 17. MAPK10 transcript analysis
#####
#####
df160 <-  Dataset_S1_Transcript

## selecting a gene term
# Selecting rows where "MAPK10" occurs anywhere in the row, in any column
df160_MAPK10 <- df160[apply(df160, 1, function(row) any(grepl("MAPK10", row, ignore.case = TRUE))), ]

# Creating a dataframe from the unique values of the 'CTAB' column
df160_MAPK10_list <- data.frame(CTAB = unique(df160_MAPK10$CTAB))

# Checking for duplicate values and removing them if found
df160_MAPK10_list <- df160_MAPK10_list[!duplicated(df160_MAPK10_list$CTAB), ]

# Display the without duplicates
df160_MAPK10_list

# Filter rows that contain the value "MAPK10" in the Gene_Symbol column
df_MAPK10 <- gene_info_transcritos[gene_info_transcritos$Display_Name == "MAPK10", ]

unique(df_MAPK10$Transcript_ID)

# Rename columns in Dataset_S1_Transcript
df160_MAPK10 <- df160_MAPK10 %>%
  rename(
    "Rank" = "Ranking",
    "Elements" = "Count_source",
    "RCD form" = "RCD_types",
    "Omic feature" = "Genotype"
  )%>%
  relocate(`Omic feature`, .after = Elements)

rio::export(df160_MAPK10, "df160_MAPK10.tsv")
rio::export(df160_MAPK10, "df160_MAPK10.xlsx")

rio::export(df160_MAPK10, "Dataset S1N.tsv")
rio::export(df160_MAPK10, "Dataset S1N.xlsx")
######
######
######
# Step 18. Estimating gene usage in transcript signatures:
# To identify top-ten most frequently contributing genes 
# Assuming the dataframe is 'df_extracted' and the variable of interest is 'gene'

# Count occurrences and calculate frequencies (as percentages)
gene_frequencies <- df_extracted %>%
  dplyr::group_by(gene) %>%
  dplyr::summarize(count = n()) %>%
  dplyr::mutate(frequency = (count / sum(count)) * 100)

export(gene_frequencies, "gene_frequencies_in_transcript_signatures.xlsx")

# Count occurrences and calculate frequencies (as percentages)
top_ten_gene_frequencies <- df_extracted %>%
  dplyr::group_by(gene) %>%
  dplyr::summarize(count = n()) %>%
  dplyr::mutate(frequency = (count / sum(count)) * 100) %>%
  dplyr::arrange(desc(count)) %>%   # Arrange in descending order by count
  dplyr::slice(1:10)  # Select the top 10 most frequent gene symbols

# Create a comma-separated list of the top 10 genes
top_ten_gene_list <- paste(top_ten_gene_frequencies$gene, collapse = ", ")

print(top_ten_gene_list)

# > top_ten_gene_list
# [1] "EFEMP2, ABI3BP, TPM1, ELN, FN1, COL1A1, DCN, PDLIM7, TCF4, COL1A2"

#####
#####
#####

# Convert gene_list to a vector of individual gene symbols
top_genes <- unlist(strsplit(top_ten_gene_list, ", "))

# Create a new dataframe by querying 'target_genes' for rows where 'Genes' matches any in the top-ten gene list
top_ten_filtered_df <- target_genes %>%
  dplyr::filter(Gene %in% top_genes)

export(top_ten_filtered_df, "top_ten_frequent_filtered_df_for_transcript_signatures.xlsx")

# Split the values in "Pathway_Categorization" by "/"
top_ten_filtered_df <- top_ten_filtered_df %>%
  dplyr::mutate(Pathway_Categorization_Split = strsplit(as.character(Pathway_Categorization), "/"))

# Relocate the "Pathway_Categorization_Split" variable after the "Gene" variable
top_ten_filtered_df <- top_ten_filtered_df %>%
  dplyr::relocate(Pathway_Categorization_Split, .after = Gene)

# If you want each split value to be in a separate row, use tidyr::unnest
top_ten_split_df <- top_ten_filtered_df %>%
  tidyr::unnest(Pathway_Categorization_Split)

# To which RCD types the top ten Genes with the highest contribution to the df
# Extract unique values, sort them alphabetically, and join them into a comma-separated string
unique_values <- top_ten_split_df %>%
  dplyr::pull(Pathway_Categorization_Split) %>%  # Extract the values from the column
  unique() %>%                                   # Get unique values
  sort() %>%                                     # Sort alphabetically
  paste(collapse = ", ")                         # Combine into a comma-separated string

print(unique_values)

#####
#####
#####
##### Step 19. Mapping Transcripts to genes on chromosomes X or Y or XY
#####
#####

target_genes <- import("Target_genes.csv")
Target_genes_X <- target_genes %>%
  filter(chromosome == "X" |chromosome == "Y" | chromosome == "X, Y")

export(Target_genes_X, "X_Y_linked_Target_genes.xlsx")

df160 <- Dataset_S1_Transcript
df173 <- Target_genes_X

# Step 1: Create a vector of unique Gene values from df173
unique_target_genes <- unique(df173$Gene)

# Step 2: Define a function to extract gene symbols from `Gene Symbol members` values in df160
extract_genes <- function(`Gene Symbol members`) {
  # Remove parentheses and fraction numbers, split by "+" to get individual gene symbols
  genes <- str_extract_all(`Gene Symbol members`, "`?([A-Za-z0-9-]+)`?")[[1]]
  # Remove backticks if present and return the clean gene symbols
  genes <- str_replace_all(genes, "`", "")
  return(genes)
}

# Step 3: Filter df160 based on matching gene symbols with df173
df160_filtered <- df160 %>%
  rowwise() %>%
  mutate(Matched_X_Y_linked_gene = list(intersect(extract_genes(`Gene Symbol members`), unique_target_genes))) %>%
  filter(length(Matched_X_Y_linked_gene) > 0) %>%
  mutate(Matched_X_Y_linked_gene = paste(Matched_X_Y_linked_gene, collapse = ", ")) %>%
  ungroup() %>%
  select(1:which(names(df160) == "Gene Symbol members"), Matched_X_Y_linked_gene, everything())

export(df160_filtered, "Transcript signatures on chr X_Y_XY.tsv")

##### Analysis of miRNA signatures
df174 <-import("Dataset_S1_miRNA.tsv")
df175 <-import("miRNA_database.xlsx") 

# Step 19a. Example dataframe df175
# Assuming df175 has columns Gene_Symbol, OriginalName, Mature1, and Mature2

# Create new rows for each non-NA value in Mature2 by duplicating corresponding rows
df175_expanded <- df175 %>%
  # Select only rows where Mature2 is not NA
  filter(!is.na(Mature2)) %>%
  # Duplicate these rows with Mature2 as the main variant
  mutate(Mature = Mature2) %>%
  # Select relevant columns
  select(Gene_Symbol, OriginalName, Mature)

# Create rows for the original Mature1 values
df175_original_mature1 <- df175 %>%
  mutate(Mature = Mature1) %>%
  select(Gene_Symbol, OriginalName, Mature)

# Combine the rows from Mature1 and Mature2 into a single dataframe
df175_final <- bind_rows(df175_original_mature1, df175_expanded)

df176 <- df175_final

# Trim whitespace in all columns of df176
df176[] <- lapply(df176, trimws)

# Trim whitespace in specific columns of df173
df173$Gene <- trimws(df173$Gene)
df173$chromosome <- trimws(df173$chromosome)

# This code trims whitespace from all columns in df176 and specific columns in df173.
# It then creates a 'chromosome' column in df176 based on matched values in df173.
# Load necessary libraries

# Merge chromosome values from df173 into df176 based on matching gene symbols
df176$chromosome <- df173$chromosome[match(df176$Gene_Symbol, df173$Gene)]

# Creating a new dataframe with rows where 'chromosome' is "X", "Y", or "XY"
df176_filtered_XY <- subset(df176, chromosome %in% c("X", "Y", "XY"))

# Creating a new dataframe with rows where 'chromosome' is equal to "X"
df176_filtered_X <- subset(df176, chromosome == "X")

# miRNA transcript map to chrX only
df177 <- df176_filtered_X
df178 <- df176_filtered_XY


######
######
###### Counting effective mature miRNA per gene locus, converting miRNA signatures in genes symbols with frequency
###### Remaking of Dataset_SXX_miRNA for the manuscriot
###### 
###### 
######
# Load necessary libraries
library(dplyr)
library(tidyr)
library(stringr)  # For str_split


# Step 20 analysing miRNA dataset
df174 <-import("Dataset_S1_miRNA.tsv")

# Create 'Unnested_miRNA' while ensuring row count remains unchanged
df174 <- df174 %>%
  mutate(
    # If 'Signature' contains nested values (inside parentheses), extract and format them
    Unnested_miRNA = ifelse(
      str_detect(Signature, "^\\(.*\\)$"),  # Detect parentheses indicating multiple values
      str_extract_all(Signature, "`([^`]+)`") %>%  # Extract values inside backticks
        lapply(function(x) paste(gsub("`", "", x), collapse = ", ")),  # Remove backticks & join values
      Signature  # Copy solo values as is
    )
  ) %>%
  mutate(Unnested_miRNA = unlist(Unnested_miRNA)) %>%  # Ensure column remains vectorized (not a list)
  relocate(Unnested_miRNA, .after = Signature)  # Move 'Unnested_miRNA' immediately after 'Signature'


######
######
###### Mapping mature miRNAs to Gene_symbols
######
######
# Load necessary libraries
library(dplyr)
library(tidyr)
library(stringr)  # For str_split

# Step 1: Ensure df175_final has no duplicate Mature values (take first occurrence of Gene_Symbol)
df175_final_unique <- df175_final %>%
  mutate(Mature = tolower(trimws(Mature))) %>%  # Trim and convert to lowercase
  group_by(Mature) %>%
  summarise(Gene_Symbol = first(Gene_Symbol), .groups = "drop")  # Take first occurrence of Gene_Symbol

# Step 2: Create a function to map miRNA values to Gene_Symbol and count occurrences
map_miRNA_to_gene <- function(miRNA_list, mapping_df) {
  # Split the miRNA list by commas and trim whitespace
  miRNA_values <- str_split(miRNA_list, ",\\s*")[[1]] %>% trimws()
  
  # Map each miRNA value to its corresponding Gene_Symbol
  gene_symbols <- sapply(miRNA_values, function(miRNA) {
    miRNA_lower <- tolower(miRNA)  # Convert to lowercase for case-insensitive matching
    matched_gene <- mapping_df$Gene_Symbol[mapping_df$Mature == miRNA_lower]
    if (length(matched_gene) > 0) matched_gene else NA
  })
  
  # Count occurrences of each gene symbol
  gene_counts <- table(gene_symbols) %>% as.data.frame()
  colnames(gene_counts) <- c("Gene_Symbol", "Count")
  
  # Create the formatted Gene_component string
  formatted_genes <- sapply(unique(gene_symbols), function(gene) {
    if (is.na(gene)) {
      return(NA)
    } else {
      count <- gene_counts$Count[gene_counts$Gene_Symbol == gene]
      paste0(gene, "(", count, "/", count, ")")
    }
  })
  
  # Combine the formatted gene symbols into a single string (comma-separated)
  paste(na.omit(formatted_genes), collapse = ", ")
}

# Step 3: Apply the mapping function to the Unnested_miRNA column
df174 <- df174 %>%
  rowwise() %>%
  mutate(Gene_component = map_miRNA_to_gene(Unnested_miRNA, df175_final_unique)) %>%
  ungroup()

# Step 4: Relocate Gene_component immediately after Unnested_miRNA
df174 <- df174 %>%
  relocate(Gene_component, .after = Unnested_miRNA)

# View the final df174 with the new Gene_component column
head(df174)

Dataset_S5_miRNA <-  df174

# Rename columns in Dataset_S5_miRNA
Dataset_S5_miRNA <- Dataset_S5_miRNA %>%
  rename(
    "Rank" = "Ranking",
    "Elements" = "Count_source",
    "RCD form" = "RCD_types",
    "Omic feature" = "Genotype"
  ) %>%
  relocate(`Omic feature`, .after = Elements)

rio::export(Dataset_S5_miRNA, "Dataset_S5_miRNA.tsv")

######
######
####### Create a lookup table with "Mature" and "Gene_Symbol" from df177
gene_map <- df177 %>%
  select(Mature, Gene_Symbol) %>%
  distinct()

# Map each unlisted miRNA symbol to its corresponding "Gene_Symbol" in df177
# # Explanation:
# - This code cleans up extracted miRNA symbols, maps them to "Gene_Symbol" values in df177,
#   removes duplicates, and formats them as "(MIRLET7A1 + MIRLET7B)".
# - The `relocate()` function positions "X_linked_Gene_symbols_in_signature" immediately after "Signature".

df174 <- df174 %>%
  rowwise() %>%
  mutate(
    X_linked_Gene_symbols_in_signature = list(
      Unnested_miRNA %>% 
        lapply(function(symbol) gene_map$Gene_Symbol[match(symbol, gene_map$Mature)]) %>% 
        unlist() %>% 
        unique() %>%                           # Remove any duplicate values
        paste(collapse = " + ")                # Format as "(MIRLET7A1 + MIRLET7B)"
    )
  ) %>%
  ungroup() %>%
  mutate(
    X_linked_Gene_symbols_in_signature = ifelse(X_linked_Gene_symbols_in_signature == "", NA, 
                                                paste0("(", X_linked_Gene_symbols_in_signature, ")")) # Add parentheses
  ) %>%
  # Step 4: Move the new column after 'Signature' and drop the intermediate column
  select(everything(), X_linked_Gene_symbols_in_signature = X_linked_Gene_symbols_in_signature, -Unnested_miRNA) %>%
  relocate(X_linked_Gene_symbols_in_signature, .after = Signature)

# Calculate the total number of rows
total_rows <- nrow(df174)

# Count rows where "X_linked_Gene_symbols_in_signature" contains "NA" as part of the string
na_count <- sum(grepl("\\bNA\\b", df174$X_linked_Gene_symbols_in_signature))

# Calculate the percentage of rows with "NA" in the specified column
na_percentage <- (na_count / total_rows) * 100

# Print the results
cat("Number of rows containing 'NA' within values:", na_count, "\n")
cat("Percentage of rows containing 'NA' in 'X_linked_Gene_symbols_in_signature':", round(na_percentage, 2), "%\n")

####
####
####
#### Step 21. Sankey diagram for MAPK10 relevant transcripts per Phenotype and CTAB
####
#### 
#### 

df180 <- df160_MAPK10
df_mapk10 <- df_MAPK10

# Extract the list of MAPK10 transcripts
mapk10_transcripts <- df_mapk10$Transcript_ID

# Unnest and filter the MAPK10 main dataset (df180) for MAPK10-specific transcripts
df_filtered <- df180 %>%
  mutate(Signature = str_remove_all(Signature, "[()]")) %>%  # Remove parentheses
  separate_rows(Signature, sep = "\\s*\\+\\s*") %>%  # Separate transcripts by '+'
  filter(Signature %in% mapk10_transcripts)  # Keep only MAPK10-specific transcripts

# Define nodes for the Sankey plot based on filtered data
unique_genes <- "MAPK10"
unique_transcripts <- unique(df_filtered$Signature)
unique_phenotypes <- unique(df_filtered$Phenotype)
unique_CTABs <- unique(df_filtered$CTAB)

# Create a combined node dataframe
nodes <- data.frame(name = c(unique_genes, unique_transcripts, unique_phenotypes, unique_CTABs))

# Create links for Gene -> Transcripts
gene_transcript_links <- df_filtered %>%
  mutate(source = match("MAPK10", nodes$name) - 1,
         target = match(Signature, nodes$name) - 1) %>%
  select(source, target) %>%
  mutate(value = 1)

# Create links for Transcripts -> Phenotypes
transcript_phenotype_links <- df_filtered %>%
  mutate(source = match(Signature, nodes$name) - 1,
         target = match(Phenotype, nodes$name) - 1) %>%
  select(source, target) %>%
  mutate(value = 1)

# Create links for Phenotypes -> CTABs
phenotype_cancer_links <- df_filtered %>%
  mutate(source = match(Phenotype, nodes$name) - 1,
         target = match(CTAB, nodes$name) - 1) %>%
  select(source, target) %>%
  mutate(value = 1)

# Combine all links
links <- bind_rows(gene_transcript_links, transcript_phenotype_links, phenotype_cancer_links)

# Generate the Sankey plot and save as an HTML widget
sankey_plot <- sankeyNetwork(Links = links, Nodes = nodes, Source = "source", Target = "target",
                             Value = "value", NodeID = "name", fontSize = 12, nodeWidth = 30)

# Generate the Sankey plot with custom CSS for font styling
sankey_plot <- sankeyNetwork(Links = links, Nodes = nodes, Source = "source", Target = "target",
                             Value = "value", NodeID = "name", fontSize = 12, nodeWidth = 30)

# Generate the Sankey plot with custom CSS for font styling
sankey_plot <- sankeyNetwork(Links = links, Nodes = nodes, Source = "source", Target = "target",
                             Value = "value", NodeID = "name", fontSize = 12, nodeWidth = 30)

# Add custom CSS to change the font style for MAPK10 node
sankey_plot <- htmlwidgets::onRender(sankey_plot, "
  function(el, x) {
    d3.select(el).selectAll('.node text')
      .style('font-family', 'Arial, Helvetica, sans-serif')
      .style('font-size', function(d) {
        return (d.name === 'MAPK10') ? '18px' : '12px';
      })
      .style('font-style', function(d) {
        return (d.name === 'MAPK10') ? 'italic' : 'normal';
      });
  }
")

sankey_plot

# Save the HTML widget with the updated font
html_file <- "MAPK10_sankey_plot.html"
saveWidget(sankey_plot, file = html_file)

library(webshot2)

# Convert to PDF and TIFF using webshot with higher zoom for better resolution
webshot(html_file, file = "MAPK10_sankey_plot.pdf", vwidth = 1100, vheight = 595, zoom = 6, delay = 0.5)  # A4 landscape at high zoom

# Save as PNG with webshot
webshot(html_file, file = "MAPK10_sankey_plot.png", vwidth = 1100, vheight = 595, zoom = 4, delay = 2)

# Convert PNG to TIFF using magick
library(magick)
png_image <- image_read("MAPK10_sankey_plot.png")
image_write(png_image, path = "MAPK10_sankey_plot.tiff", format = "tiff", density = 600)

####
####
####
#### Step 22. Sankey diagram for MAPK10 relevant transcripts per Phenotype and CTAB with connection numbers in parentheses
####
####
####

library(dplyr)
library(stringr)
library(tidyr)
library(networkD3)
library(htmlwidgets)
library(webshot)
library(magick)

# Load datasets
df180 <- df160_MAPK10
df_mapk10 <- df_MAPK10

# Extract the list of MAPK10 transcripts
mapk10_transcripts <- df_mapk10$Transcript_ID

# Unnest and filter the MAPK10 main dataset (df180) for MAPK10-specific transcripts
df_filtered <- df180 %>%
  mutate(Signature = str_remove_all(Signature, "[()]") %>% str_trim()) %>%
  separate_rows(Signature, sep = "\\s*\\+\\s*") %>%
  filter(Signature %in% mapk10_transcripts)

# Define nodes for the Sankey plot based on filtered data
unique_genes <- "MAPK10"
unique_transcripts <- unique(df_filtered$Signature)

# Convert transcripts to a dataframe, extract numbers, and sort in descending order
unique_transcripts <- data.frame(Signature = unique_transcripts, stringsAsFactors = FALSE) %>%
  mutate(Transcript_Num = as.numeric(str_extract(Signature, "\\d+"))) %>%
  arrange(desc(Transcript_Num)) %>%
  pull(Signature)  # Convert back to a sorted vector

unique_phenotypes <- unique(df_filtered$Phenotype)
unique_CTABs <- sort(unique(df_filtered$CTAB))  # Sort CTABs alphabetically

# Create a combined node dataframe ensuring CTAB nodes appear last in sorted order
nodes <- data.frame(name = c(unique_genes, unique_transcripts, unique_phenotypes, unique_CTABs),
                    stringsAsFactors = FALSE) %>%
  mutate(node_id = row_number() - 1)  # Assign numeric node IDs after sorting

# Update link definitions to match the correctly sorted node order

# Create links for Gene -> Transcripts
gene_transcript_links <- df_filtered %>%
  mutate(source = match("MAPK10", nodes$name) - 1,
         target = match(Signature, nodes$name) - 1) %>%
  select(source, target) %>%
  mutate(value = 1)

# Create links for Transcripts -> Phenotypes
transcript_phenotype_links <- df_filtered %>%
  mutate(source = match(Signature, nodes$name) - 1,
         target = match(Phenotype, nodes$name) - 1) %>%
  select(source, target) %>%
  mutate(value = 1)

# Create links for Phenotypes -> CTABs (Ensuring alphabetical order is preserved)
phenotype_cancer_links <- df_filtered %>%
  mutate(source = match(Phenotype, nodes$name) - 1,
         target = match(CTAB, nodes$name) - 1) %>%
  select(source, target) %>%
  mutate(value = 1)

# Combine all links
links <- bind_rows(gene_transcript_links, transcript_phenotype_links, phenotype_cancer_links)

# Adjust node sizes based on number of connections
cancer_counts <- links %>% group_by(target) %>% summarise(value = sum(value))

# Ensure all nodes have a default value before merging
nodes <- nodes %>%
  mutate(value = 1) %>%  # Initialize all nodes with value = 1
  left_join(cancer_counts, by = c("node_id" = "target")) %>%
  mutate(value = coalesce(value.y, value.x)) %>%
  select(-value.x, -value.y)

# Append identifiers to nodes while maintaining parentheses globally
nodes <- nodes %>%
  mutate(name = paste0(name, " (", value, ")"))

# Generate the Sankey plot and save as an HTML widget
sankey_plot <- sankeyNetwork(Links = links, Nodes = nodes, Source = "source", Target = "target",
                             Value = "value", NodeID = "name", fontSize = 12, nodeWidth = 30)

# Add custom CSS to change the font style for MAPK10 node
sankey_plot <- htmlwidgets::onRender(sankey_plot, "
  function(el, x) {
    d3.select(el).selectAll('.node text')
      .style('font-family', 'Arial, Helvetica, sans-serif')
      .style('font-size', function(d) {
        return (d.name.includes('MAPK10')) ? '18px' : '12px';
      })
      .style('font-style', function(d) {
        return (d.name.includes('MAPK10')) ? 'italic' : 'normal';
      });
  }
")

sankey_plot

# Save the HTML widget
html_file <- "MAPK10_sankey_plot.html"
saveWidget(sankey_plot, file = html_file)

# Convert to PDF and TIFF using webshot
webshot(html_file, file = "MAPK10_sankey_plot.pdf", vwidth = 1100, vheight = 595, zoom = 6, delay = 0.5)
webshot(html_file, file = "MAPK10_sankey_plot.png", vwidth = 1100, vheight = 595, zoom = 4, delay = 2)

# Convert PNG to TIFF using magick
png_image <- image_read("MAPK10_sankey_plot.png")
image_write(png_image, path = "MAPK10_sankey_plot.tiff", format = "tiff", density = 600)


####
####
####
#### Step 23. Filtering rows by exact match in the Display_Name column
####
####
####
# Define the genes of interest
genes_of_interest <- c("COL1A1", "UMOD", "CD36", "ABI3BP", "TCF4")
data <- gene_info_transcritos

# Create separate dataframes for each gene based on exact match in 'Display_Name'
for (gene in genes_of_interest) {
  assign(paste0("df_", gene), data[data$Display_Name == gene, ])
}

# Display each dataframe for Part 1
print(df_COL1A1)
print(df_UMOD)
print(df_CD36)
print(df_ABI3BP)
print(df_TCF4)

# Part 2: Filtering rows for specific gene symbols where gene appears anywhere in a row in 'df160'

# Define the same genes of interest
genes_of_interest <- c("COL1A1", "UMOD", "CD36", "ABI3BP", "TCF4")

# Define the dataframe
df160 <- Dataset_S1_Transcript

# Create separate dataframes for each gene based on presence in any column
for (gene in genes_of_interest) {
  assign(paste0("df160_", gene), df160[apply(df160, 1, function(row) any(grepl(gene, row, ignore.case = TRUE))), ])
}

# Display each dataframe for Part 2
print(df160_COL1A1)
print(df160_UMOD)
print(df160_CD36)
print(df160_ABI3BP)
print(df160_TCF4)

#######
#######
#######
####### Step 24. Generate a corrected Sankey data for COL1A1 and UMOD with proportional cancer type node sizes
#######
#######
#######

# Adding the specific signature identifiers
# Note: The thickness of the connective stroke lines is cumulative,
# based on all the signatures that contain the specific gene-transcript per CTAB.
# In the examples provided, only one signature had either the 13/13 (COL1A1) or 12/12 (UMOD)
# transcripts in each signature.

# Load necessary libraries
library(dplyr)
library(stringr)
library(tidyr)
library(networkD3)
library(htmlwidgets)
library(webshot)
library(magick)

# Load datasets
data <- gene_info_transcritos
df160 <- Dataset_S1_Transcript

# Initialize global node list and link list
nodes_list <- list()
links_list <- list()
global_node_counter <- 0

# Function to generate nodes and links for a given gene
create_sankey_data <- function(test_gene) {
  message("Processing gene: ", test_gene)
  
  # Retrieve known transcripts
  df_known_transcripts <- data %>% filter(Display_Name == test_gene)
  known_transcripts <- df_known_transcripts$Transcript_ID
  
  # Retrieve gene signature dataset
  df_gene_signatures <- df160 %>%
    filter(if_any(everything(), ~ str_detect(., regex(test_gene, ignore_case = TRUE))))
  
  # Expand and filter dataset for gene-specific transcripts
  df_filtered <- df_gene_signatures %>%
    mutate(Signature = str_remove_all(Signature, "[()]") %>% str_trim()) %>%
    separate_rows(Signature, sep = "\\s*\\+\\s*") %>%
    mutate(Signature = str_remove_all(Signature, "`")) %>%
    filter(Signature %in% known_transcripts)
  
  # Validate correct cancer type mapping
  if (test_gene == "UMOD") {
    df_filtered <- df_filtered %>%
      mutate(CTAB = ifelse(CTAB == "CHOL", "KICH", CTAB))
  }
  
  # Define nodes
  gene_node <- data.frame(name = test_gene, group = "Gene", node = global_node_counter, value = 1)
  global_node_counter <<- global_node_counter + 1
  
  transcript_nodes <- data.frame(
    name = unique(df_filtered$Signature), group = "Transcript",
    node = global_node_counter:(global_node_counter + length(unique(df_filtered$Signature)) - 1),
    value = table(df_filtered$Signature)[unique(df_filtered$Signature)]
  )
  global_node_counter <<- global_node_counter + nrow(transcript_nodes)
  
  phenotype_nodes <- data.frame(
    name = unique(df_filtered$Phenotype), group = "Phenotype",
    node = global_node_counter:(global_node_counter + length(unique(df_filtered$Phenotype)) - 1),
    value = 1
  )
  global_node_counter <<- global_node_counter + nrow(phenotype_nodes)
  
  CTAB_nodes <- data.frame(
    name = unique(df_filtered$CTAB), group = "Cancer Type",
    node = global_node_counter:(global_node_counter + length(unique(df_filtered$CTAB)) - 1),
    value = table(df_filtered$CTAB)[unique(df_filtered$CTAB)]
  ) %>% arrange(name)
  
  global_node_counter <<- global_node_counter + nrow(CTAB_nodes)
  
  nodes <- bind_rows(gene_node, transcript_nodes, phenotype_nodes, CTAB_nodes)
  nodes_list <<- bind_rows(nodes_list, nodes)
  
  # Define links
  gene_transcript_links <- df_filtered %>%
    mutate(source = gene_node$node,
           target = transcript_nodes$node[match(Signature, transcript_nodes$name)],
           value = 1) %>%
    select(source, target, value)
  
  transcript_phenotype_links <- df_filtered %>%
    mutate(source = transcript_nodes$node[match(Signature, transcript_nodes$name)],
           target = phenotype_nodes$node[match(Phenotype, phenotype_nodes$name)],
           value = 1) %>%
    select(source, target, value)
  
  phenotype_cancer_links <- df_filtered %>%
    mutate(source = phenotype_nodes$node[match(Phenotype, phenotype_nodes$name)],
           target = CTAB_nodes$node[match(CTAB, CTAB_nodes$name)],
           value = 1) %>%
    select(source, target, value)
  
  links <- bind_rows(gene_transcript_links, transcript_phenotype_links, phenotype_cancer_links)
  links_list <<- bind_rows(links_list, links)
  
  message("Finished processing: ", test_gene)
}

# Generate data for both genes
create_sankey_data("COL1A1")
create_sankey_data("UMOD")

# Convert lists to dataframes
nodes <- nodes_list
links <- links_list

# Adjust node sizes based on number of connections
cancer_counts <- links %>% group_by(target) %>% summarise(value = sum(value))

# Ensure all nodes have a default value before merging
nodes <- nodes %>%
  mutate(value = 1) %>%  # Initialize all nodes with value = 1
  left_join(cancer_counts, by = c("node" = "target")) %>%
  mutate(value = coalesce(value.y, value.x)) %>%
  select(-value.x, -value.y)

# Append identifiers to specific CTAB names while maintaining parentheses globally, but excluding them from genes
nodes <- nodes %>%
  mutate(name = case_when(
    name == "HNSC" ~ paste0("(HNSC-308.5.3.N.3.0.0.3.2.3) HNSC (", value, ")"),
    name == "KICH" ~ paste0("(KICH-117.5.3.N.2.0.0.2.4.3) KICH (", value, ")"),
    group == "Gene" ~ name,  # Exclude parentheses for gene nodes
    TRUE ~ paste0(name, " (", value, ")")
  ))

# Define color mapping
color_mapping <- 'd3.scaleOrdinal()
                      .domain(["Gene", "Signature", "Phenotype", "CTAB", "TSM", "MSI", "TMB"])
                      .range(["#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#90ee90", "#ff69b4", "#ff0000"])'

# Create Sankey plot
sankey_plot <- sankeyNetwork(
  Links = links, Nodes = nodes, Source = "source", Target = "target",
  Value = "value", NodeID = "name", fontSize = 12, nodeWidth = 30,
  colourScale = color_mapping
)

# Add custom CSS to change the font for specific gene nodes to a larger, italicized style
sankey_plot <- htmlwidgets::onRender(sankey_plot, "
      function(el, x) {
        d3.select(el).selectAll('.node text')
          .style('font-family', 'Arial, Helvetica, sans-serif')
          .style('font-size', function(d) {
            // Check if the node name is COL1A1 or UMOD
            return (d.name === 'COL1A1' || d.name === 'UMOD') ? '18px' : '12px';
          })
          .style('font-style', function(d) {
            // Apply italic style to COL1A1 and UMOD nodes
            return (d.name === 'COL1A1' || d.name === 'UMOD') ? 'italic' : 'normal';
          });
      }
    ")

sankey_plot

# Save output
html_file <- "COL1A1_UMOD_corrected_sankey_plot.html"
saveWidget(sankey_plot, file = html_file)
webshot(html_file, file = "COL1A1_UMOD_corrected_sankey_plot.png", vwidth = 1100, vheight = 595, zoom = 4, delay = 2)
png_image <- image_read("COL1A1_UMOD_corrected_sankey_plot.png")
image_write(png_image, path = "COL1A1_UMOD_corrected_sankey_plot.tiff", format = "tiff", density = 600)
# Convert to PDF using webshot with higher zoom for better resolution
webshot(html_file, file = "COL1A1_UMOD_corrected_sankey_plot.pdf", vwidth = 1100, vheight = 595, zoom = 6, delay = 0.5)  # A4 landscape at high zoom


####
#### Step 25. Sankey diagram for MAPK10 relevant transcripts per Phenotype and CTAB with connection numbers in parentheses
#### Transcript nodes sorted by descending order; CTAB nodes sorted by alphabetical order
#### Selected signatsures ctab_identifiers look up and plotted
#### Figure 7. Sankey diagram illustrating transcript-specific associations of the MAPK10 gene across various phenotypes and cancer types
####
####
library(dplyr)
library(stringr)
library(tidyr)
library(networkD3)
library(htmlwidgets)
library(webshot)
library(magick)
library(purrr)

# Load datasets
df180 <- df160_MAPK10
df_mapk10 <- df_MAPK10

# Extract the list of MAPK10 transcripts
mapk10_transcripts <- df_mapk10$Transcript_ID

# Unnest and filter the MAPK10 main dataset (df180) for MAPK10-specific transcripts
df_filtered <- df180 %>%
  mutate(Signature = str_remove_all(Signature, "[()]") %>% str_trim()) %>%
  separate_rows(Signature, sep = "\\s*\\+\\s*") %>%
  filter(Signature %in% mapk10_transcripts)

# Define nodes for the Sankey plot based on filtered data
unique_genes <- "MAPK10"
unique_transcripts <- unique(df_filtered$Signature)

# Extract the numeric part of the transcript ID and sort in ascending order
unique_transcripts <- data.frame(Signature = unique_transcripts, stringsAsFactors = FALSE) %>%
  mutate(Transcript_Num = as.numeric(str_extract(Signature, "\\d+"))) %>%
  arrange(Transcript_Num) %>%
  pull(Signature)  # Convert back to a sorted vector

unique_phenotypes <- unique(df_filtered$Phenotype)
unique_CTABs <- sort(unique(df_filtered$CTAB))  # Sort CTABs alphabetically

# Create a combined node dataframe ensuring CTAB nodes appear last in sorted order
nodes <- data.frame(name = c(unique_genes, unique_transcripts, unique_phenotypes, unique_CTABs),
                    stringsAsFactors = FALSE) %>%
  mutate(node_id = row_number() - 1)  # Assign numeric node IDs after sorting

# Update link definitions to match the correctly sorted node order

# Create links for Gene -> Transcripts
gene_transcript_links <- df_filtered %>%
  mutate(source = match("MAPK10", nodes$name) - 1,
         target = match(Signature, nodes$name) - 1) %>%
  select(source, target) %>%
  mutate(value = 1)

# Create links for Transcripts -> Phenotypes
transcript_phenotype_links <- df_filtered %>%
  mutate(source = match(Signature, nodes$name) - 1,
         target = match(Phenotype, nodes$name) - 1) %>%
  select(source, target) %>%
  mutate(value = 1)

# Create links for Phenotypes -> CTABs (Ensuring alphabetical order is preserved)
phenotype_cancer_links <- df_filtered %>%
  mutate(source = match(Phenotype, nodes$name) - 1,
         target = match(CTAB, nodes$name) - 1) %>%
  select(source, target) %>%
  mutate(value = 1)

# Combine all links
links <- bind_rows(gene_transcript_links, transcript_phenotype_links, phenotype_cancer_links)

# Adjust node sizes based on number of connections
cancer_counts <- links %>% group_by(target) %>% summarise(value = sum(value))

# Ensure all nodes have a default value before merging
nodes <- nodes %>%
  mutate(value = 1) %>%  # Initialize all nodes with value = 1
  left_join(cancer_counts, by = c("node_id" = "target")) %>%
  mutate(value = coalesce(value.y, value.x)) %>%
  select(-value.x, -value.y)

# Define a lookup table for CTAB names and their corresponding identifiers
ctab_identifiers <- c(
  "LUAD-350.5.3.N.2.0.0.1.4.2" = "(LUAD-350.5.3.N.2.0.0.1.4.2)",
  "LUSC-1549.5.2.P.1.4.0.4.4.2" = "(LUSC-1549.5.2.P.1.4.0.4.4.2)",
  "LUAD-1824.5.1.N.1.0.0.3.4.2" = "(LUAD-1824.5.1.N.1.0.0.3.4.2)",
  "LGG-1814.5.3.P.3.93.72.2.3.2" = "(LGG-1814.5.3.P.3.93.72.2.3.2)",
  "STAD-1718.5.3.N.1.44.0.3.4.2" = "(STAD-1718.5.3.N.1.44.0.3.4.2)"
)

# Append identifiers to specific CTAB names while maintaining parentheses globally, but excluding them from genes
nodes <- nodes %>%
  mutate(
    ctab_identifiers_matched = map_chr(name, function(ctab) {
      matches <- ctab_identifiers[grepl(paste0("^", ctab), names(ctab_identifiers))]  # Match identifiers that start with CTAB
      if (length(matches) > 0) {
        if (ctab == "LUAD-350.5.3.N.2.0.0.1.4.2" || ctab == "LUAD-1824.5.1.N.1.0.0.3.4.2") {
          paste0(matches, collapse = "\n")  # Append all matching identifiers with newline for specific CTABs
        } else {
          paste0(matches, collapse = " ")  # Append all matching identifiers
        }
      } else {
        NA_character_
      }
    }),
    name = case_when(
      !is.na(ctab_identifiers_matched) ~ paste0(ctab_identifiers_matched, " ", name, " (", value, ")"), # Append identifiers
      name == "MAPK10" ~ name,  # Keep MAPK10 node unchanged
      TRUE ~ paste0(name, " (", value, ")") # Default formatting
    )
  ) %>%
  select(-ctab_identifiers_matched)  # Remove temporary column

# Generate the Sankey plot and save as an HTML widget
sankey_plot <- sankeyNetwork(Links = links, Nodes = nodes, Source = "source", Target = "target",
                             Value = "value", NodeID = "name", fontSize = 12, nodeWidth = 30)

# Add custom CSS to change the font style for MAPK10 node
sankey_plot <- htmlwidgets::onRender(sankey_plot, "
  function(el, x) {
    d3.select(el).selectAll('.node text')
      .style('font-family', 'Arial, Helvetica, sans-serif')
      .style('font-size', function(d) {
        return (d.name === 'MAPK10') ? '18px' : '12px';
      })
      .style('font-style', function(d) {
        return (d.name === 'MAPK10') ? 'italic' : 'normal';
      });
  }
")

sankey_plot

# Save the HTML widget
html_file <- "MAPK10_sankey_plot.html"
saveWidget(sankey_plot, file = html_file)

# Convert to PDF and TIFF using webshot
webshot(html_file, file = "MAPK10_sankey_plot.pdf", vwidth = 1500, vheight = 650, zoom = 6, delay = 0.5)
webshot(html_file, file = "MAPK10_sankey_plot.png", vwidth = 1500, vheight = 650, zoom = 4, delay = 2)

# Convert PNG to TIFF using magick
png_image <- image_read("MAPK10_sankey_plot.png")
image_write(png_image, path = "MAPK10_sankey_plot.tiff", format = "tiff", density = 600)

### 
### 
### 
### Step 26. Dynamic network plot
### MAPK10_interactive_network_plot
### 
### 
### 

# Install and load necessary libraries
if (!requireNamespace("visNetwork", quietly = TRUE)) {
  install.packages("visNetwork")
}
library(visNetwork)
library(dplyr)
library(htmlwidgets)
library(tidyr)

# Prepare edges for each relationship type
# Signature -> Signature
display_transcript_edges <- df_filtered %>%
  select(Signature, Signature) %>%
  rename(from = Signature, to = Signature)

# Signature -> Phenotype
transcript_phenotype_edges <- df_filtered %>%
  select(Signature, Phenotype) %>%
  rename(from = Signature, to = Phenotype)

# Phenotype -> CTAB
phenotype_cancer_edges <- df_filtered %>%
  select(Phenotype, CTAB) %>%
  rename(from = Phenotype, to = CTAB)

# Create MAPK10 -> Signature connections
mapk10_edges <- df_filtered %>%
  select(Signature) %>%
  distinct() %>%  # Ensure unique signatures
  mutate(from = "MAPK10", to = Signature)  # Explicitly assign MAPK10 as source

# Count the number of unique signatures MAPK10 is connected to
mapk10_degree <- n_distinct(mapk10_edges$to)  # Count unique "to" values (Signatures)

# Combine all edges and remove NA values
all_edges <- bind_rows(display_transcript_edges, transcript_phenotype_edges, phenotype_cancer_edges, mapk10_edges) %>%
  drop_na(from, to) %>%  # Ensure no missing values
  distinct()

# Ensure MAPK10 is explicitly included in the nodes dataframe
nodes <- data.frame(id = unique(c(all_edges$from, all_edges$to)))  # Unique node IDs

# Calculate node degree (connection count)
degree_data <- all_edges %>%
  tidyr::gather(key = "type", value = "node", from, to) %>%
  count(node)

# Assign node attributes (size, color, group)
nodes <- nodes %>%
  left_join(degree_data, by = c("id" = "node")) %>%
  mutate(
    label = id,
    size = ifelse(id == "MAPK10", mapk10_degree * 2 + 10, n * 2 + 10),  # Dynamically scale MAPK10 size
    group = case_when(
      id == "MAPK10" ~ "MAPK10",
      id %in% df_filtered$Signature ~ "Signature",
      id %in% df_filtered$Phenotype ~ "Phenotype",
      id %in% df_filtered$CTAB ~ "CTAB",
      TRUE ~ "Other"
    ),
    color = case_when(
      id == "MAPK10" ~ "red",  # MAPK10 in red for visibility
      group == "Signature" ~ "skyblue",
      group == "Phenotype" ~ "lightgreen",
      group == "CTAB" ~ "yellow",
      TRUE ~ "grey"
    )
  )

# Create edges dataframe with standardized format
edges <- all_edges %>%
  rename(from = from, to = to) %>%
  mutate(width = 1)  # Uniform width for simplicity

# Generate interactive network plot
network_plot <- visNetwork(nodes, edges, height = "800px", width = "100%") %>%
  visNodes(size = nodes$size) %>%
  visEdges(smooth = TRUE) %>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
  visInteraction(dragNodes = TRUE, dragView = TRUE, zoomView = TRUE) %>%
  visLayout(randomSeed = 123)

# Display the network
network_plot

# Save as HTML
saveWidget(network_plot, file = "MAPK10_interactive_network_plot.html", selfcontained = TRUE)

#####
##### Step 27. Lookup for some specific genes
#####
# Select all rows that contain "SOX2" in any column
df182 <- df181[apply(df181, 1, function(row)  any(grepl("\\(SOX2\\(\\d+/\\d+\\)\\)", row))), ]

# Select all rows that contain "NANOG" in any column
df183 <- df181[apply(df181, 1, function(row)  any(grepl("\\(NANOG\\(\\d+/\\d+\\)\\)", row))), ]

# Select all rows that contain "POU5F1" in any column
df184 <- df181[apply(df181, 1, function(row)  any(grepl("\\(POU5F1\\(\\d+/\\d+\\)\\)", row))), ]

# Select all rows that contain "POU5F1" in any column
df185 <- df181[apply(df181, 1, function(row)  any(grepl("\\(TP53\\(\\d+/\\d+\\)\\)", row))), ]

# Select all rows that contain "POU5F1" in any column
df186 <- df181[apply(df181, 1, function(row)  any(grepl("\\(CDH1\\(\\d+/\\d+\\)\\)", row))), ]


#####
#####
#####
#####
# Load necessary libraries
library(dplyr)
library(tidyr)
library(stringr)

Dataset_S7_Protein <- import("Dataset_S7_Protein.tsv")
potein_info <- import("protein_info.xlsx")
df190 <- Dataset_S7_Protein
df191 <- protein_info


# Step 1: Clean and unnest the "Signature" column in df190
df190_unnested <- df190 %>%
  mutate(Signature_clean = str_replace_all(Signature, "[()]", ""), # Remove parentheses
         Signature_clean = str_replace_all(Signature_clean, " ", "")) %>% # Remove spaces
  separate_rows(Signature_clean, sep = "\\+") # Split by "+" and unnest

# Step 2: Map Protein1 values to Gene_Symbol using df191
df190_unnested <- df190_unnested %>%
  left_join(df191, by = c("Signature_clean" = "Protein1")) # Map proteins to gene symbols

# Step 3: Group by the original row index and concatenate Gene_Symbols
df190_unnested <- df190_unnested %>%
  group_by(across(-c(Signature_clean, Gene_Symbol))) %>% # Group by all columns except Signature_clean and Gene_Symbol
  summarise(`Gene Symbol members` = paste(Gene_Symbol, collapse = ", "), .groups = "drop") # Concatenate gene symbols

# Step 4: Count occurrences of each gene symbol in the string
df190_unnested <- df190_unnested %>%
  mutate(`Gene Symbol members` = sapply(str_split(`Gene Symbol members`, ", "), function(x) {
    counts <- table(x) # Count occurrences of each gene symbol
    paste(names(counts), ifelse(counts > 1, paste0("(", counts, ")"), ""), sep = "", collapse = ", ")
  }))

# Step 5: Insert the "Gene Symbol members" column before the "Elements" column
df190_unnested <- df190_unnested %>%
  relocate(`Gene Symbol members`, .before = Elements)

Dataset_S7_Protein_final <- df190_unnested
rio::export(Dataset_S7_Protein_final, "Dataset_S7_Protein_final.tsv")
###
###
###
###

######
######
######
# Save the current workspace to an RData file
save.image("workspace.RData")
######
######
######


###
###
###
###
# Load the saved workspace into the current session
load("workspace.RData")
###
###
###
###
