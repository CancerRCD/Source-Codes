library(rio)
setwd("C:/Users/Emamnuell/Desktop/panoptosis/genes/genes_24_optosis/codigos/4°_signature/genes")

# Assuming the data frames are already loaded into the R environment
final_cox_results_genes <- import("final_cox_results_genes.xlsx")
big_table_13C <- import("big_table_13C.xlsx")

# View the first rows of the data frames
head(big_table_13C)
head(final_cox_results_genes)

# Merge the data frames
merged_df <- merge(final_cox_results_genes, big_table_13C, by.x = "genes", by.y = "name", all.x = TRUE)

# Select the desired columns
final_cox_results_genes <- merged_df[, c("genes", "cancer_types", "Members_Pathway_Categorization", "Pathway_Categorization",
                                         "correlation", "Correlation_p.adj", "log10_correlation", "Expression",
                                         "Expression_p.sgnif", "Expression_p.adj", "log10_p.adj", "var_fenotipica",
                                         "var_genotipica", "p.value_Cox_OS", "Type_Cox_OS", "p.value_Cox_DSS",
                                         "Type_Cox_DSS", "p.value_Cox_DFI", "Type_Cox_DFI", "p.value_Cox_PFI", "Type_Cox_PFI")]

# View the first rows of the resulting data frame
head(final_cox_results_genes)
export(final_cox_results_genes, "final_cox_results_genes_com_vias.xlsx")
