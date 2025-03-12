# Load the necessary libraries
library(dplyr)
library(rio)

setwd("C:/Users/Emamnuell/Desktop/panoptosis/genes/genes_24_optosis/codigos/4°_signature/transcript")

# Load the data
final_cox_results_transcript <- import("final_cox_results_transcript.xlsx")
gene_info_trancritos_excel <- import("gene_info_trancritos.xlsx")
big_table_13C_manually_curated <- import("big_table_13C_manually_curated.xlsx")

# Check the structure of the tables
str(final_cox_results_transcript)
str(gene_info_trancritos_excel)

# Merge the tables based on the columns Transcript and Transcript_ID
merged_data <- merge(final_cox_results_transcript, 
                     gene_info_trancritos_excel[, c("Transcript_ID", "Display_Name")], 
                     by.x = "Transcript", 
                     by.y = "Transcript_ID", 
                     all.x = TRUE)

# Rename the column Display_Name to Gene_Symbol
names(merged_data)[names(merged_data) == "Display_Name"] <- "Gene_Symbol"

# Reorder the columns to place Gene_Symbol first
merged_data <- merged_data[, c("Gene_Symbol", setdiff(names(merged_data), "Gene_Symbol"))]

# Merge the tables merged_data and big_table_13C_manually_curated
merged_data_final <- merge(merged_data, 
                           big_table_13C_manually_curated[, c("name", "Driver", "Members_Pathway_Categorization", "Pathway_Categorization")], 
                           by.x = "Gene_Symbol", 
                           by.y = "name", 
                           all.x = TRUE)

# Define the column order as specified
column_order <- c("Gene_Symbol", "Transcript" ,"Driver", "Members_Pathway_Categorization", "Pathway_Categorization", "var_genotipica", "var_fenotipica", "cancer_types",   
                  "correlation", "Correlation_p.adj", "log10_correlation", 
                  "Expression", "Expression_p.sgnif", "Expression_p.adj", "log10_p.adj",  
                  "p.value_Cox_OS", "Type_Cox_OS", "p.value_Cox_DSS", "Type_Cox_DSS", 
                  "p.value_Cox_DFI", "Type_Cox_DFI", "p.value_Cox_PFI", "Type_Cox_PFI")

# Reorder the columns in the final data frame
merged_data_final <- merged_data_final[, column_order]

# Check the result of the merging
head(merged_data_final)

# Save the merged table
export(merged_data_final,"final_cox_results_transcript_com_vias.xlsx")
