library(writexl)
library(readxl)

# Set the working directory and import the Excel file
setwd("C:/Users/Emamnuell/Desktop/panoptosis/genes/genes_24_optosis/codigos/4°_signature/miRNA")
signature_genes <- read_excel("signature_miRNA.xlsx")

# Lists of possible Phenotype and Genotype categories
phenotype_categories <- c("TMB", "MSI", "Stemness")
genotype_categories <- c("Mutation", "mRNA", "Methylation", "CNV", "Protein", "Transcript", "miRNA")

# Export each dataframe containing results
for (phenotype in phenotype_categories) {
  for (genotype in genotype_categories) {
    # Dataframe name
    df_name <- paste("miRNA_signature", phenotype, genotype, sep = "_")
    
    # Filter the original dataframe (signature_genes) based on Phenotype and Genotype
    filtered_df <- subset(signature_genes, Phenotype == phenotype & Genotype == genotype)
    
    # Check if the filtered dataframe has results
    if (nrow(filtered_df) > 0) {
      # Build the full path of the Excel file
      file_path <- file.path("C:/Users/Emamnuell/Desktop/panoptosis/genes/genes_24_optosis/codigos/4°_signature/miRNA", paste0(df_name, ".xlsx"))
      
      # Export the filtered dataframe to an Excel file using writexl
      write_xlsx(filtered_df, file_path)
      
      # Optional: Check the structure of the created dataframe
      cat("Exported data frame:", df_name, ".xlsx\n")
      print(head(filtered_df))  # Show the first rows of the dataframe
      cat("\n")
    } else {
      cat("No results found for:", df_name, "\n\n")
    }
  }
}
