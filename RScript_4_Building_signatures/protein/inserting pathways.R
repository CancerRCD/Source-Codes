library(rio)
setwd("C:/Users/Emamnuell/Desktop/panoptosis/genes/genes_24_optosis/codigos/4°_signature/protein")

# Assuming the data frames are already loaded into the R environment
final_cox_results_protein <- import("final_cox_results_protein.xlsx")
big_table_13C <- import("RPPA_proteins_in_xena_mapping.tsv")

# View the first rows of the data frames
head(big_table_13C)
head(final_cox_results_protein)

# Merge the data frames
merged_df <- merge(final_cox_results_protein, big_table_13C, by.x = "genes", by.y = "Protein", all.x = TRUE)
colnames(merged_df)[colnames(merged_df) == "genes"] <- "Protein"

# Match the "Protein" column in merged_df with big_table_13C and add the corresponding Gene_Symbol
merged_df$genes <- big_table_13C$Gene_Symbol[match(merged_df$Protein, big_table_13C$Protein)]

# Select the desired columns
final_cox_results_Protein <- merged_df[, c("genes", "Protein", "cancer_types", "Members_Pathway_Categorization", "Pathway_Categorization",
                                           "correlation", "Correlation_p.adj", "log10_correlation", "Expression",
                                           "Expression_p.sgnif", "Expression_p.adj", "log10_p.adj", "var_fenotipica",
                                           "var_genotipica", "p.value_Cox_OS", "Type_Cox_OS", "p.value_Cox_DSS",
                                           "Type_Cox_DSS", "p.value_Cox_DFI", "Type_Cox_DFI", "p.value_Cox_PFI", "Type_Cox_PFI")]

# View the first rows of the resulting data frame
head(final_cox_results_Protein)
export(final_cox_results_Protein,"final_cox_results_protein_com_vias.xlsx")
