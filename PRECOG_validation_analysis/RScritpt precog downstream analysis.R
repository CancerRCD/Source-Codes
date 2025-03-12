####
####
#### PRECOG validation
#### 
#### 

setwd()

library(data.table)


genes_df = fread("Precog_Gene_Symbols_top_signatures_df1158.tsv", data.table = F)
precog_metaZ = fread("precog_metazscore.csv",data.table = F)
top_sig_precog = precog_metaZ[precog_metaZ$Gene %in% genes_df$Gene_Symbol,]
#write.table(top_sig_precog, "PRECOG_Top_signature_df1158_intersection.tsv", sep = "\t", row.names = F, quote = F)

# dataset 1 has the genes in the rows and columns the PRECOG meta-z score for different tumors in the columns
names(top_sig_precog)

# dataset 2 has the genes combination (signatures) in the rows and columns are metrics but i just care about the combination of genes in the signature
sig_df = fread("df1158.tsv", data.table = F)
head(sig_df$Signature, 20)

# Load datasets
sig_df <- fread("df1158.tsv")
top_sig_precog <- fread("top_sig_precog.tsv")

# Function to extract individual genes from signatures
extract_genes <- function(signature) {
  signature <- gsub("[()]", "", signature) # Remove parentheses
  genes <- unlist(strsplit(signature, " \\+ "))
  return(genes)
}

# Initialize result list
aggregated_meta_z_list <- list()

# Prepare dataset for meta-Z calculations
sel_sig <- as.data.frame(top_sig_precog)
rownames(sel_sig) <- sel_sig$Gene
sel_sig$Gene <- sel_sig$Name <- NULL

# Compute mean and standard deviation of meta-Z scores per tumor type using apply
meta_z_means <- apply(sel_sig, 2, function(x) mean(x, na.rm = TRUE))
meta_z_sds <- apply(sel_sig, 2, function(x) sd(x, na.rm = TRUE))

# Iterate over each signature
for (signature in sig_df$Signature) {
  genes <- extract_genes(signature)
  
  # Filter dataset 1 to get rows matching the genes in the signature
  matching_genes <- sel_sig[rownames(sel_sig) %in% genes, , drop = FALSE]
  
  if (nrow(matching_genes) > 0) {
    # Compute the mean meta-Z score for each tumor type
    mean_meta_z <- colMeans(matching_genes, na.rm = TRUE) 
    
    # Compute Z-score: (X - mean) / sd, handling division by zero
    z_scores <- ifelse(meta_z_sds == 0, 0, (mean_meta_z - meta_z_means) / meta_z_sds)
    
    # Store the result
    aggregated_meta_z_list[[signature]] <- data.table(Tumor = names(z_scores), Z_Score = z_scores)
  }
}

# Convert result list to data.table
aggregated_meta_z_dt <- rbindlist(lapply(names(aggregated_meta_z_list), function(sig) {
  data.table(Signature = sig, aggregated_meta_z_list[[sig]])
}), fill = TRUE)

names(aggregated_meta_z_dt)
aggregated_meta_z_dt = dcast(aggregated_meta_z_dt, Signature~Tumor, value.var = "Z_Score")

dd = sig_df[!(sig_df$Signature %in% aggregated_meta_z_dt$Signature),]
dd = as.data.frame(table(sig_df$Signature))

# Save or display
write.table(aggregated_meta_z_dt, "PRECOG_Signatures_aggregated_meta_z_scores.tsv", sep = "\t",
            row.names = F, quote = F)


#########
#
# Median
#
##########

# Load datasets
genes_df = fread("Precog_Gene_Symbols_top_signatures_df1158.tsv", data.table = F)
precog_metaZ = fread("precog_metazscore.csv", data.table = F)
top_sig_precog = precog_metaZ[precog_metaZ$Gene %in% genes_df$Gene_Symbol,]

sig_df = fread("df1158.tsv", data.table = F)

# Function to extract individual genes from signatures
extract_genes <- function(signature) {
  signature <- gsub("[()]", "", signature) # Remove parentheses
  genes <- unlist(strsplit(signature, " \\+ "))
  return(genes)
}

# Initialize result list
aggregated_meta_z_list <- list()

# Prepare dataset for median meta-Z calculations
sel_sig <- as.data.frame(top_sig_precog)
rownames(sel_sig) <- sel_sig$Gene
sel_sig$Gene <- sel_sig$Name <- NULL

# Iterate over each signature
for (signature in sig_df$Signature) {
  genes <- extract_genes(signature)
  
  # Filter dataset to get rows matching the genes in the signature
  matching_genes <- sel_sig[rownames(sel_sig) %in% genes, , drop = FALSE]
  
  if (nrow(matching_genes) > 0) {
    # Compute the median meta-Z score for each tumor type
    median_meta_z <- apply(matching_genes, 2, median, na.rm = TRUE)
    
    # Store the result
    aggregated_meta_z_list[[signature]] <- data.table(Tumor = names(median_meta_z), Median_Z = median_meta_z)
  }
}

# Convert result list to data.table
aggregated_meta_z_dt <- rbindlist(lapply(names(aggregated_meta_z_list), function(sig) {
  data.table(Signature = sig, aggregated_meta_z_list[[sig]])
}), fill = TRUE)

# Reshape to wide format
aggregated_meta_z_dt = dcast(aggregated_meta_z_dt, Signature ~ Tumor, value.var = "Median_Z")

# Save output
write.table(aggregated_meta_z_dt, "PRECOG_Signatures_aggregated_median_meta_z_scores.tsv", sep = "\t",
            row.names = F, quote = F)

aggregated_meta_z_dt = aggregated_meta_z_dt[with(aggregated_meta_z_dt, order(aggregated_meta_z_dt$Unweighted_meta.Z_of_all_cancers)),]
order_all = unique(aggregated_meta_z_dt$Signature)


library(ggplot2)
plot1 = ggplot(aggregated_meta_z_dt, aes(x = factor(Signature, levels = order_all), 
                                         y = Unweighted_meta.Z_of_all_cancers, 
                                         fill = Unweighted_meta.Z_of_all_cancers)) + 
  geom_bar(stat = "identity") + 
  theme_minimal() +
  
  # Horizontal reference lines at -3 and 3
  geom_hline(yintercept = c(-3, 3), linetype = "dashed", color = "black") +
  
  # Gradient fill from blue (-3) to white (0) to red (3)
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0, limits = c(-3, 3), oob = scales::squish) +
  
  theme(
    axis.text.x = element_blank(),  # Hide gene labels
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    panel.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank()
  ) +
  
  labs(title = "PRECOG Meta Z Scores across all cancers per signature", 
       y = "PRECOG Meta Z Score", fill = "Meta Z Score")


png("PRECOG_Signature_prognostic_markers.png", 
    width = 35, height = 20, units = "cm", res = 300, bg = "white")
plot1
dev.off()



###
###
### PRECOG downstream analysis on top most clinically meaningful signatures
### "Meaningful Risky" and "Meaningful Protective"
### Precog pre-analysis of the top most clinically meaningful signatures done by Ronaldo
#shell.exec("https://precog.stanford.edu/")

##### PART A - finding which Signatures values are missing from the precog analysis datasets.

# # Define required packages
required_packages <- c("ggplot2", "reshape2", "dplyr", "tidyr", "readxl", 
                       "gridExtra", "grid", "magick", "rio")

# Install only missing packages
missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(missing_packages) > 0) {
  install.packages(missing_packages)
}

# Load all required packages
invisible(lapply(required_packages, library, character.only = TRUE))

# Print confirmation message
cat("\n✅ All required packages are installed (if needed) and loaded successfully!\n")

library(rio)        # For importing and exporting files
library(ggplot2)    # For plotting
library(reshape2)   # For data transformation
library(dplyr)      # For data manipulation
library(tidyr)      # For data tidying
library(readxl)     # For reading Excel files
library(gridExtra)  # For arranging multiple plots
library(grid)       # For graphical layout elements
library(magick)     # For image processing

#### PART B: Load the Datasets ####
# Load TCGA CTAB cancer types dataset
CTAB <- import("Dataset_S1I.xlsx")

df3001 <- import("PRECOG_Signatures_aggregated_median_meta_z_scores.tsv")
df3002 <- import("PRECOG_Signatures_aggregated_meta_z_scores.tsv")
df1158 <- import("df1158.tsv") # from df1158 <-df1157_mRNA_R_P_4_STT/ n = 126
## NOte: most clinical meaningful signatures for validation by PRECOG analysis
df3003 <- unique(df1158$Signature)

# Create a dataframe with column names of df3001
df3001_columns <- data.frame(Column_Names = names(df3001))

# Remove the first four and last columns
df3001_columns <- df3001_columns[-c(1,41), , drop = FALSE]

# Call PRECOG cancer types dataset
precog_df <- df3001_columns

rio::export(df3001_columns, "df3001_precog_cancer_types.xlsx")

#### PART C: Load Required Libraries ####
library(readxl)
library(dplyr)
library(tidyr)
library(rio)

# Extract TCGA cancer type names and abbreviations
tcga_cancer_types <- CTAB[[2]]  # First column: Full Cancer Name
tcga_ctab <- CTAB[[1]]  # Second column: TCGA Abbreviation

unique_CTAB_df1158 <- unique(df1158$CTAB)

# Select rows where CTAB values match the specified cancer types
selected_rows_CTAB_df1158 <- df1158 %>%
  filter(CTAB %in% c("LUAD", "BRCA", "CESC", "LUSC", "KIRP", "HNSC", "LGG", "BLCA", "PRAD", "STAD", "ACC"))

# Extract PRECOG cancer type names
precog_cancer_types <- precog_df[[1]]  # First column: PRECOG Cancer Name

#### PART C: Define the Correspondence Mapping ####
tcga_precog_mapping <- list(
  "Adrenocortical carcinoma" = "Adrenocortical_cancer",
  "Bladder Urothelial Carcinoma" = "Bladder_cancer",
  "Brain Lower Grade Glioma" = c("Brain_cancer_Astrocytoma", "Brain_cancer_Glioma"),
  "Glioblastoma multiforme" = "Brain_cancer_Glioblastoma",
  "Breast invasive carcinoma" = "Breast_cancer",
  "Colon adenocarcinoma" = "Colon_cancer",
  "Stomach adenocarcinoma" = "Gastric_cancer",
  "Testicular Germ Cell Tumors" = "Germ_cell_tumors",
  "Head and Neck squamous cell carcinoma" = c("Head_and_neck_cancer", "Head_and_neck_cancer_Hypopharyngeal_cancer", "Head_and_neck_cancer_Oral_SCC"),
  "Esophageal carcinoma" = "Head_and_neck_cancer_Oesophageal_cancer",
  "Acute Myeloid Leukemia" = "Hematopoietic_cancer_AML",
  "Lymphoid Neoplasm Diffuse Large B-cell Lymphoma" = "Hematopoietic_cancer_DLBCL",
  "Kidney renal clear cell carcinoma" = "Kidney_cancer",
  "Kidney renal papillary cell carcinoma" = "Kidney_cancer",
  "Kidney Chromophobe" = "Kidney_cancer",
  "Liver hepatocellular carcinoma" = c("Liver_cancer", "Liver_cancer_Primary"),
  "Lung adenocarcinoma" = "Lung_cancer_ADENO",
  "Lung squamous cell carcinoma" = "Lung_cancer_SCC",
  "Skin Cutaneous Melanoma" = c("Melanoma", "Melanoma_Metastasis"),
  "Mesothelioma" = "Mesothelioma",
  "Ovarian serous cystadenocarcinoma" = "Ovarian_cancer",
  "Pancreatic adenocarcinoma" = "Pancreatic_cancer",
  "Prostate adenocarcinoma" = "Prostate_cancer",
  "Sarcoma" = c("Sarcoma_Ewing_sarcoma", "Sarcoma_Osteosarcoma")
)
#### PART D: Create the `CTAB` DataFrame ####
# Convert mapping into a structured dataframe
tcga_mapping_list <- lapply(names(tcga_precog_mapping), function(tcga_type) {
  precog_types <- tcga_precog_mapping[[tcga_type]]
  ctab_value <- tcga_ctab[which(tcga_cancer_types == tcga_type)]
  
  # If multiple PRECOG correspondences, create multiple rows
  if (is.vector(precog_types)) {
    data.frame(CTAB = ctab_value, TCGA_Cancer_Type = tcga_type, PRECOG_Cancer_Type = precog_types)
  } else {
    data.frame(CTAB = ctab_value, TCGA_Cancer_Type = tcga_type, PRECOG_Cancer_Type = precog_types)
  }
})

# Combine all mappings into a single dataframe
CTAB_final <- bind_rows(tcga_mapping_list)

# Select rows where CTAB values match the specified cancer types
selected_rows_CTAB_final  <- CTAB_final %>%
  filter(CTAB %in% c("LUAD", "BRCA", "CESC", "LUSC", "KIRP", "HNSC", "LGG", "BLCA", "PRAD", "STAD", "ACC"))

#### PART D: Verification of Mapping ####
# Check for missing TCGA types
missing_tcga <- setdiff(names(tcga_precog_mapping), tcga_cancer_types)

# Check for missing PRECOG types
missing_precog <- setdiff(unlist(tcga_precog_mapping), precog_cancer_types)

# Create verification dataframe
verification_results <- data.frame(
  Category = c("Missing TCGA Cancer Types", "Missing PRECOG Cancer Types"),
  Missing_Values = c(paste(missing_tcga, collapse = ", "), paste(missing_precog, collapse = ", "))
)

# Print verification results
print("Verification Results:")
print(verification_results)

#### PART E: Save the Finalized `CTAB` ####

rio::export(CTAB_final,"Final_TCGA_PRECOG_Mapping.xlsx" )
rio::export(CTAB_final,"Final_TCGA_PRECOG_Mapping.tsv" )

rio::export(CTAB_final,"Dataset_S1I_Final_TCGA_PRECOG_Mapping.xlsx" )
rio::export(CTAB_final,"Dataset_S1I_Final_TCGA_PRECOG_Mapping.tsv" )

# Extract TCGA cancer type names and abbreviations
tcga_cancer_types <- CTAB_final[[2]]  # First column: Full Cancer Name
tcga_CTAB_finaltcga_CTAB_final <- CTAB_final[[1]]  # Second column: TCGA Abbreviation

# Extract PRECOG cancer type names
precog_cancer_types <- precog_df[[1]]  # First column: PRECOG Cancer Name
#
#### PART F: Define the Correspondence Mapping ####
tcga_precog_mapping <- list(
  "Adrenocortical carcinoma" = "Adrenocortical_cancer",
  "Bladder Urothelial Carcinoma" = "Bladder_cancer",
  "Brain Lower Grade Glioma" = c("Brain_cancer_Astrocytoma", "Brain_cancer_Glioma"),
  "Glioblastoma multiforme" = "Brain_cancer_Glioblastoma",
  "Breast invasive carcinoma" = "Breast_cancer",
  "Colon adenocarcinoma" = "Colon_cancer",
  "Stomach adenocarcinoma" = "Gastric_cancer",
  "Testicular Germ Cell Tumors" = "Germ_cell_tumors",
  "Head and Neck squamous cell carcinoma" = c("Head_and_neck_cancer", "Head_and_neck_cancer_Hypopharyngeal_cancer", "Head_and_neck_cancer_Oral_SCC"),
  "Esophageal carcinoma" = "Head_and_neck_cancer_Oesophageal_cancer",
  "Acute Myeloid Leukemia" = "Hematopoietic_cancer_AML",
  "Lymphoid Neoplasm Diffuse Large B-cell Lymphoma" = "Hematopoietic_cancer_DLBCL",
  "Kidney renal clear cell carcinoma" = "Kidney_cancer",
  "Kidney renal papillary cell carcinoma" = "Kidney_cancer",
  "Kidney Chromophobe" = "Kidney_cancer",
  "Liver hepatocellular carcinoma" = c("Liver_cancer", "Liver_cancer_Primary"),
  "Lung adenocarcinoma" = "Lung_cancer_ADENO",
  "Lung squamous cell carcinoma" = "Lung_cancer_SCC",
  "Skin Cutaneous Melanoma" = c("Melanoma", "Melanoma_Metastasis"),
  "Mesothelioma" = "Mesothelioma",
  "Ovarian serous cystadenocarcinoma" = "Ovarian_cancer",
  "Pancreatic adenocarcinoma" = "Pancreatic_cancer",
  "Prostate adenocarcinoma" = "Prostate_cancer",
  "Sarcoma" = c("Sarcoma_Ewing_sarcoma", "Sarcoma_Osteosarcoma")
)

####
####
#### PART G: Create the `CTAB_final` DataFrame ####
# Convert mapping into a structured dataframe
tcga_mapping_list <- lapply(names(tcga_precog_mapping), function(tcga_type) {
  precog_types <- tcga_precog_mapping[[tcga_type]]
  CTAB_final_value <- tcga_CTAB_final[which(tcga_cancer_types == tcga_type)]
  
  # If multiple PRECOG correspondences, create multiple rows
  if (is.vector(precog_types)) {
    data.frame(CTAB_final = CTAB_final_value, TCGA_Cancer_Type = tcga_type, PRECOG_Cancer_Type = precog_types)
  } else {
    data.frame(CTAB_final = CTAB_final_value, TCGA_Cancer_Type = tcga_type, PRECOG_Cancer_Type = precog_types)
  }
})

# Combine all mappings into a single dataframe
CTAB_final <- bind_rows(tcga_mapping_list)

#### PART H: Verification of Mapping ####
# Check for missing TCGA types
missing_tcga <- setdiff(names(tcga_precog_mapping), tcga_cancer_types)

# Check for missing PRECOG types
missing_precog <- setdiff(unlist(tcga_precog_mapping), precog_cancer_types)

# Create verification dataframe
verification_results <- data.frame(
  Category = c("Missing TCGA Cancer Types", "Missing PRECOG Cancer Types"),
  Missing_Values = c(paste(missing_tcga, collapse = ", "), paste(missing_precog, collapse = ", "))
)

# Print verification results
print("Verification Results:")
print(verification_results)

#### PART I: Save the Finalized `CTAB_final` ####

rio::export(CTAB_final,"Final_TCGA_PRECOG_Mapping.xlsx" )
rio::export(CTAB_final,"Final_TCGA_PRECOG_Mapping.tsv" )

precog_query_genes <- import("Precog_Gene_Symbols_top_signatures_df1158.tsv")


# Ensure both dataframes exist
if (exists("df1158") & exists("df3002")) {
  
  # Trim whitespace from the "Signature" column in both dataframes
  df1158$Signature <- trimws(df1158$Signature)
  df3002$Signature <- trimws(df3002$Signature)
  
  # Identify missing values (present in df1158 but absent in df3002)
  missing_signatures <- setdiff(df1158$Signature, df3002$Signature)
  df1158_missing_in_precog <- df1158[df1158$Signature %in% missing_signatures, ]
  
  # Identify common values (present in both df1158 and df3002)
  common_signatures <- intersect(df1158$Signature, df3002$Signature)
  df1158_present_in_precog <- df1158[df1158$Signature %in% common_signatures, ]
  
  # Display the resulting dataframes
  print("Missing values in df1158 (not in df3002):")
  print(df1158_missing_in_precog)
  
  print("Common values in df1158 and df3002:")
  print(df1158_present_in_precog)
  
} else {
  stop("One or both of the dataframes (df1158, df3002) do not exist.")
}

print(missing_signatures)

#> missing_signatures
#[1] "ADAMTS9-AS1"
#
##### PART J - adding CTAB and Nomenclature values to df3002

# Ensure both dataframes exist
if (exists("df1158") & exists("df3002")) {
  
  # Trim whitespace from the "Signature" column in both dataframes
  df1158$Signature <- trimws(df1158$Signature)
  df3002$Signature <- trimws(df3002$Signature)
  
  df1158$Nomenclature <- trimws(df1158$Nomenclature)
  
  # Perform the lookup and handle duplicate values by merging
  df3002_expanded <- merge(df3002, df1158[, c("Signature", "CTAB", "Nomenclature", "Combined_Outcome")], 
                           by = "Signature", all.x = TRUE)
  
  # Relocate 'CTAB', 'Nomenclature', and 'Combined_Outcome' to be immediately after 'Signature'
  desired_order <- c("Signature", "CTAB", "Nomenclature", "Combined_Outcome", setdiff(names(df3002_expanded), c("Signature", "CTAB", "Nomenclature", "Combined_Outcome")))
  df3002_expanded <- df3002_expanded[, desired_order]
  
  # Assign the expanded dataframe to df3002
  df3002 <- df3002_expanded
  
  # Display the updated dataframe
  print("Updated df3002 with expanded CTAB, Nomenclature, and Combined_Outcome:")
  print(df3002)
  
} else {
  stop("One or both of the dataframes (df1158, df3002) do not exist.")
}

# Create a dataframe with column names of df3002
df3002_columns <- data.frame(Column_Names = names(df3002))

# Remove the first four rows
df3002_columns <- df3002_columns[-c(1:4,41), , drop = FALSE]

rio::export(df3002_columns, "df3002_precog_cancer_types.xlsx")

##### PART K - TCGA Mapping

library(dplyr)

# Define the selected columns
selected_columns <- c(
  "Signature", "CTAB", "Nomenclature", "Combined_Outcome",
  "Adrenocortical_cancer", "Bladder_cancer", "Breast_cancer", "Colon_cancer",
  "Hematopoietic_cancer_DLBCL", "Head_and_neck_cancer_Oesophageal_cancer",
  "Brain_cancer_Glioblastoma", "Brain_cancer_Astrocytoma", "Brain_cancer_Glioma",
  "Head_and_neck_cancer", "Head_and_neck_cancer_Hypopharyngeal_cancer",
  "Head_and_neck_cancer_Oral_SCC", "Kidney_cancer", "Kidney_cancer", "Kidney_cancer",
  "Hematopoietic_cancer_AML", "Liver_cancer", "Liver_cancer_Primary",
  "Lung_cancer_ADENO", "Lung_cancer_SCC", "Melanoma", "Melanoma_Metastasis",
  "Mesothelioma", "Ovarian_cancer", "Pancreatic_cancer", "Prostate_cancer",
  "Sarcoma_Ewing_sarcoma", "Sarcoma_Osteosarcoma"
)

# Ensure unique column names to prevent duplicates
selected_columns <- unique(selected_columns)

# Select the specified columns from df3002 to create df3004
df3004 <- df3002 %>% select(all_of(selected_columns))

# Export df3004 for further use
rio::export(df3004, "df3004_selected_columns.xlsx")

####
####
####
#### PART L -Heatmap plot
#### 
#### 
#### 

# Copy df3004 to df3005 for modifications
df3005 <- df3004  

# Select only the columns with meta-z scores (variables 5 to 30)
meta_z_columns <- colnames(df3005)[5:30]

# Ensure meta_z_columns exist in df3005
if (!all(meta_z_columns %in% colnames(df3005))) {
  stop("Some meta-z columns are missing in df3005.")
}

# Filter for meta-z scores greater than 3.09 OR less than -3.09
df_filtered <- df3005 %>%
  select(Nomenclature, all_of(meta_z_columns)) %>%
  filter_at(vars(all_of(meta_z_columns)), any_vars(. > 3.09 | . < -3.09))

# Check if the dataframe is empty after filtering
if (nrow(df_filtered) == 0) {
  stop("No rows satisfy the |meta-z| > 3.09 or < -3.09 filter. Check your filtering criteria.")
}

# Reshape data from wide to long format
df_long <- df_filtered %>%
  pivot_longer(cols = all_of(meta_z_columns), names_to = "Cancer_Type", values_to = "Meta_Z_Score")

# **Force alphabetical order of Nomenclature (y-axis)**
df_long$Nomenclature <- factor(df_long$Nomenclature, levels = sort(unique(df_long$Nomenclature)))

# **Identify cells that should be highlighted (|Meta_Z_Score| > 3.09 or < -3.09)**
df_long$Highlight <- ifelse(df_long$Meta_Z_Score > 3.09 | df_long$Meta_Z_Score < -3.09, "Yes", "No")

# Define Meta-Z limits dynamically to avoid NA issues
meta_z_max <- max(abs(df_long$Meta_Z_Score), na.rm = TRUE)

# Create the heatmap with FULL black borders for significant cells
heatmap_plot <- ggplot(df_long, aes(x = Cancer_Type, y = Nomenclature, fill = Meta_Z_Score)) +
  geom_tile(color = "white") +  # Keeps gridlines visible
  geom_tile(data = df_long %>% filter(Highlight == "Yes"), fill = NA, color = "black", linewidth = 0.8) +  # Full black border for significant cells
  scale_fill_gradient2(low = "blue", mid = "grey98", high = "red", midpoint = 0, 
                       limits = c(-meta_z_max, meta_z_max),
                       name = "Meta-Z Score") +
  theme_minimal() +
  labs(title = "Heatmap of Prognostic Meta-Z Scores", 
       x = "Cancer Type", 
       y = "Signatures") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right")

# Print the heatmap
print(heatmap_plot)

# Define file paths for saving outputs
pdf_file <- "df3002_Heatmap_Prognostic_MetaZ_Highlighted.pdf"
tiff_file <- "df3002_Heatmap_Prognostic_MetaZ_Highlighted.tiff"

# Save as high-resolution PDF
ggsave(filename = pdf_file, plot = heatmap_plot, width = 14, height = 6, dpi = 600, device = cairo_pdf)

# Save as high-resolution TIFF (600 DPI, LZW compression)
ggsave(filename = tiff_file, plot = heatmap_plot, width = 14, height = 6, 
       dpi = 600, device = "tiff", compression = "lzw")

#######
#######
####### PART M - df3001-based heatmap  plot
####### 
####### 
#####

### Precog downstream analysis

##### PART M - finding which Signatures values are missing from the precog analysis datasets.

# Ensure both dataframes exist
if (exists("df1158") & exists("df3001")) {
  
  # Trim whitespace from the "Signature" column in both dataframes
  df1158$Signature <- trimws(df1158$Signature)
  df3001$Signature <- trimws(df3001$Signature)
  
  # Identify missing values (present in df1158 but absent in df3001)
  missing_signatures <- setdiff(df1158$Signature, df3001$Signature)
  df1158_missing_in_precog <- df1158[df1158$Signature %in% missing_signatures, ]
  
  # Identify common values (present in both df1158 and df3001)
  common_signatures <- intersect(df1158$Signature, df3001$Signature)
  df1158_present_in_precog <- df1158[df1158$Signature %in% common_signatures, ]
  
  # Display the resulting dataframes
  print("Missing values in df1158 (not in df3001):")
  print(df1158_missing_in_precog)
  
  print("Common values in df1158 and df3001:")
  print(df1158_present_in_precog)
  
} else {
  stop("One or both of the dataframes (df1158, df3001) do not exist.")
}

print(missing_signatures)

#> missing_signatures
#[1] "ADAMTS9-AS1"
#
##### PART N - adding CTAB and Nomenclature values to df3001

# Ensure both dataframes exist
if (exists("df1158") & exists("df3001")) {
  
  # Trim whitespace from the "Signature" column in both dataframes
  df1158$Signature <- trimws(df1158$Signature)
  df3001$Signature <- trimws(df3001$Signature)
  
  df1158$Nomenclature <- trimws(df1158$Nomenclature)
  
  # Perform the lookup and handle duplicate values by merging
  df3001_expanded <- merge(df3001, df1158[, c("Signature", "CTAB", "Nomenclature", "Combined_Outcome")], 
                           by = "Signature", all.x = TRUE)
  
  # Relocate 'CTAB', 'Nomenclature', and 'Combined_Outcome' to be immediately after 'Signature'
  desired_order <- c("Signature", "CTAB", "Nomenclature", "Combined_Outcome", setdiff(names(df3001_expanded), c("Signature", "CTAB", "Nomenclature", "Combined_Outcome")))
  df3001_expanded <- df3001_expanded[, desired_order]
  
  # Assign the expanded dataframe to df3001
  df3001 <- df3001_expanded
  
  # Display the updated dataframe
  print("Updated df3001 with expanded CTAB, Nomenclature, and Combined_Outcome:")
  print(df3001)
  
} else {
  stop("One or both of the dataframes (df1158, df3001) do not exist.")
}

# Select the specified columns from df3001 to create df3006
df3006 <- df3001 %>% select(all_of(selected_columns))

# Export df3006 for further use
rio::export(df3006, "df3006_selected_columns.xlsx")

####
####
####
#### PART O - Heatmap plot based on df3001, precog median meta Z sacores
#### 
#### 
#### 

# Copy df3001 to df3007 for modifications
df3007 <- df3006  

# Select only the columns with meta-z scores (variables 5 to 30)
meta_z_columns <- colnames(df3007)[5:30]

# Ensure meta_z_columns exist in df3007
if (!all(meta_z_columns %in% colnames(df3007))) {
  stop("Some meta-z columns are missing in df3007.")
}

# Filter for meta-z scores greater than 3.09 OR less than -3.09
df_filtered <- df3007 %>%
  select(Nomenclature, all_of(meta_z_columns)) %>%
  filter_at(vars(all_of(meta_z_columns)), any_vars(. > 3.09 | . < -3.09))

# Check if the dataframe is empty after filtering
if (nrow(df_filtered) == 0) {
  stop("No rows satisfy the |meta-z| > 3.09 or < -3.09 filter. Check your filtering criteria.")
}

# Reshape data from wide to long format
df_long <- df_filtered %>%
  pivot_longer(cols = all_of(meta_z_columns), names_to = "Cancer_Type", values_to = "Meta_Z_Score")

# **Force alphabetical order of Nomenclature (y-axis)**
df_long$Nomenclature <- factor(df_long$Nomenclature, levels = sort(unique(df_long$Nomenclature)))

# **Identify cells that should be highlighted (|Meta_Z_Score| > 3.09 or < -3.09)**
df_long$Highlight <- ifelse(df_long$Meta_Z_Score > 3.09 | df_long$Meta_Z_Score < -3.09, "Yes", "No")

# Define Meta-Z limits dynamically to avoid NA issues
meta_z_max <- max(abs(df_long$Meta_Z_Score), na.rm = TRUE)

# Create the heatmap with FULL black borders for significant cells
heatmap_plot <- ggplot(df_long, aes(x = Cancer_Type, y = Nomenclature, fill = Meta_Z_Score)) +
  geom_tile(color = "white") +  # Keeps gridlines visible
  geom_tile(data = df_long %>% filter(Highlight == "Yes"), fill = NA, color = "black", linewidth = 0.8) +  # Full black border for significant cells
  scale_fill_gradient2(low = "blue", mid = "grey98", high = "red", midpoint = 0, 
                       limits = c(-meta_z_max, meta_z_max),
                       name = "Meta-Z Score") +
  theme_minimal() +
  labs(title = "Heatmap of Prognostic Meta-Z Scores", 
       x = "Cancer Type", 
       y = "Signatures") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right")

# Print the heatmap
print(heatmap_plot)

# Define file paths for saving outputs
pdf_file <- "df3001_Heatmap_Prognostic_MetaZ_Highlighted.pdf"
tiff_file <- "df3001_Heatmap_Prognostic_MetaZ_Highlighted.tiff"

# Save as high-resolution PDF
ggsave(filename = pdf_file, plot = heatmap_plot, width = 16, height = 30, dpi = 600, device = cairo_pdf)

# Save as high-resolution TIFF (600 DPI, LZW compression)
ggsave(filename = tiff_file, plot = heatmap_plot, width = 16, height = 30, 
       dpi = 600, device = "tiff", compression = "lzw")


######
######
######
######
#######
####### PART P - df3001-based heatmap plot with asterisk for matching CTAB and PREGOC cancer types
####### 
####### 
#####

### Precog downstream analysis

##### PART P - finding which Signatures values are missing from the precog analysis datasets.

# Ensure both dataframes exist
if (exists("df1158") & exists("df3001")) {
  
  # Trim whitespace from the "Signature" column in both dataframes
  df1158$Signature <- trimws(df1158$Signature)
  df3001$Signature <- trimws(df3001$Signature)
  
  # Identify missing values (present in df1158 but absent in df3001)
  missing_signatures <- setdiff(df1158$Signature, df3001$Signature)
  df1158_missing_in_precog <- df1158[df1158$Signature %in% missing_signatures, ]
  
  # Identify common values (present in both df1158 and df3001)
  common_signatures <- intersect(df1158$Signature, df3001$Signature)
  df1158_present_in_precog <- df1158[df1158$Signature %in% common_signatures, ]
  
  # Display the resulting dataframes
  print("Missing values in df1158 (not in df3001):")
  print(df1158_missing_in_precog)
  
  print("Common values in df1158 and df3001:")
  print(df1158_present_in_precog)
  
} else {
  stop("One or both of the dataframes (df1158, df3001) do not exist.")
}

print(missing_signatures)

#> missing_signatures
#[1] "ADAMTS9-AS1"
#
##### PART Q - adding CTAB and Nomenclature values to df3001

# Ensure both dataframes exist
if (exists("df1158") & exists("df3001")) {
  
  # Trim whitespace from the "Signature" column in both dataframes
  df1158$Signature <- trimws(df1158$Signature)
  df3001$Signature <- trimws(df3001$Signature)
  
  df1158$Nomenclature <- trimws(df1158$Nomenclature)
  
  # Perform the lookup and handle duplicate values by merging
  df3001_expanded <- merge(df3001, df1158[, c("Signature", "CTAB", "Nomenclature", "Combined_Outcome")], 
                           by = "Signature", all.x = TRUE)
  
  # Relocate 'CTAB', 'Nomenclature', and 'Combined_Outcome' to be immediately after 'Signature'
  desired_order <- c("Signature", "CTAB", "Nomenclature", "Combined_Outcome", setdiff(names(df3001_expanded), c("Signature", "CTAB", "Nomenclature", "Combined_Outcome")))
  df3001_expanded <- df3001_expanded[, desired_order]
  
  # Assign the expanded dataframe to df3001
  df3001 <- df3001_expanded
  
  # Display the updated dataframe
  print("Updated df3001 with expanded CTAB, Nomenclature, and Combined_Outcome:")
  print(df3001)
  
} else {
  stop("One or both of the dataframes (df1158, df3001) do not exist.")
}

# Select the specified columns from df3001 to create df3006
df3006 <- df3001 %>% select(all_of(selected_columns))

# Export df3006 for further use
rio::export(df3006, "df3006_selected_columns.xlsx")

####
####
####
#### PART R - Heatmap plot based on df3001, precog median meta Z scores
#### 
#### 
#### 

# Copy df3001 to df3007 for modifications
df3007 <- df3006  

# Select only the columns with meta-z scores (variables 5 to 30)
meta_z_columns <- colnames(df3007)[5:30]

# Ensure meta_z_columns exist in df3007
if (!all(meta_z_columns %in% colnames(df3007))) {
  stop("Some meta-z columns are missing in df3007.")
}

# Filter for meta-z scores greater than 3.09 OR less than -3.09
df_filtered <- df3007 %>%
  select(Nomenclature, CTAB, all_of(meta_z_columns)) %>%
  filter_at(vars(all_of(meta_z_columns)), any_vars(. > 3.09 | . < -3.09))

# Check if the dataframe is empty after filtering
if (nrow(df_filtered) == 0) {
  stop("No rows satisfy the |meta-z| > 3.09 or < -3.09 filter. Check your filtering criteria.")
}

# Reshape data from wide to long format
df_long <- df_filtered %>%
  pivot_longer(cols = all_of(meta_z_columns), names_to = "Cancer_Type", values_to = "Meta_Z_Score")

# **Force alphabetical order of Nomenclature (y-axis)**
df_long$Nomenclature <- factor(df_long$Nomenclature, levels = sort(unique(df_long$Nomenclature)))

# **Identify cells that should be highlighted (|Meta_Z_Score| > 3.09 or < -3.09)**
df_long$Highlight <- ifelse(df_long$Meta_Z_Score > 3.09 | df_long$Meta_Z_Score < -3.09, "Yes", "No")

# **Check CTAB correspondence using CTAB_final dataframe**
# Merge df_long with CTAB_final to check correspondence
df_long <- df_long %>%
  left_join(CTAB_final, by = c("Cancer_Type" = "PRECOG_Cancer_Type")) %>%
  mutate(Asterisk = ifelse(CTAB.x == CTAB.y & Highlight == "Yes", "*", ""))

# Define Meta-Z limits dynamically to avoid NA issues
meta_z_max <- max(abs(df_long$Meta_Z_Score), na.rm = TRUE)

# Create the heatmap with FULL black borders for significant cells and asterisks
heatmap_plot <- ggplot(df_long, aes(x = Cancer_Type, y = Nomenclature, fill = Meta_Z_Score)) +
  geom_tile(color = "white") +  # Keeps gridlines visible
  geom_tile(data = df_long %>% filter(Highlight == "Yes"), fill = NA, color = "black", linewidth = 0.8) +  # Full black border for significant cells
  geom_text(aes(label = Asterisk), color = "black", size = 6, vjust = 0.8) +  # Add asterisks to highlighted cells
  scale_fill_gradient2(low = "blue", mid = "grey98", high = "red", midpoint = 0, 
                       limits = c(-meta_z_max, meta_z_max),
                       name = "Meta-Z Score") +
  theme_minimal() +
  labs(title = "Heatmap of Prognostic Meta-Z Scores", 
       x = "Cancer Type", 
       y = "Signatures") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right")

# Print the heatmap
print(heatmap_plot)

# Define file paths for saving outputs
pdf_file <- "df3001_Heatmap_Prognostic_MetaZ_Highlighted.pdf"
tiff_file <- "df3001_Heatmap_Prognostic_MetaZ_Highlighted.tiff"

# Save as high-resolution PDF
ggsave(filename = pdf_file, plot = heatmap_plot, width = 16, height = 30, dpi = 600, device = cairo_pdf)

# Save as high-resolution TIFF (600 DPI, LZW compression)
ggsave(filename = tiff_file, plot = heatmap_plot, width = 16, height = 30, 
       dpi = 600, device = "tiff", compression = "lzw")

# PART S - **Verification of CTAB Correspondence**
# Check if the mapping between CTAB and Cancer_Type is accurate
verification <- df_long %>%
  filter(Asterisk == "*") %>%  # Select rows with asterisks
  mutate(Prognosis = ifelse(Meta_Z_Score > 3.09, "Risky", 
                            ifelse(Meta_Z_Score < -3.09, "Protective", NA))) %>%  # Add Prognosis column
  select(Nomenclature, Cancer_Type, CTAB.x, CTAB.y, Meta_Z_Score, Prognosis)  # Display relevant columns

# Print the verification results
print("Verification of CTAB Correspondence:")
print(verification)

# **Count the Number of Cells with Asterisks**
# Count the number of cells that meet the Meta_Z_Score requirements and have asterisks
asterisk_count <- sum(df_long$Asterisk == "*", na.rm = TRUE)

####  Verification table
# **Add Combined_Outcome to Verification Dataframe**
# Merge verification with df1158 to add Combined_Outcome
verification <- verification %>%
  left_join(df1158 %>% select(Nomenclature, Combined_Outcome), 
            by = "Nomenclature")

# Print the updated verification results
print("Verification of CTAB Correspondence with Combined_Outcome:")
print(verification)
# Print the count
print(paste("Number of cells with asterisks:", asterisk_count))

rio::export(verification, "PRECOG_verification.xlsx")

#####
#####
#####
#####
#####
# Save the entire R workspace
save.image(file = "my_workspace_2.RData")
####
####
####
####
  
#####
#####
#####
#####
# Load the previously saved workspace if needed
load("my_workspace_2.RData")

