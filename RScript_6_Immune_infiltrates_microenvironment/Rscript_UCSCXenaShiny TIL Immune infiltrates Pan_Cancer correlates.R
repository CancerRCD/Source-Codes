## UCSCXenaShiny TIL Immune infiltrates Pan_Cancer correlates
## Authors: Higor Almeida
## Last update: 16/07/2024

# Loading the necessary libraries
library(UCSCXenaShiny)
library(UCSCXenaTools)
library(dplyr)
library(tidyr)
library(purrr)
library(rio)

# Setting the directory
setwd("C:/Users/Emamnuell/Desktop/laboratorio/monografia/dados/9°_TIL")

# Importing the file
Table <- import("signature_transcript.tsv")

# Removing duplicates
Target_genes <- Table %>%
  distinct(Transcript)

# Function to process compound and simple gene names
process_gene_name <- function(gene) {
  if (grepl("\\+", gene)) {
    return(gene)  # Keep the gene as it is if there is a plus sign
  } else {
    return(gsub("`|\\(|\\)", "", gene))  # Remove backticks and parentheses if there is no plus sign
  }
}

# Initializing a vector to track genes that failed due to lack of data
null_genes <- c()

# Initializing a data frame to store correlation results
til_results <- data.frame()

# Loop to analyze TIL gene correlations
for (i in seq_along(Target_genes$Transcript)) {
  gene <- Target_genes$Transcript[i]
  gene_processed <- process_gene_name(gene)
  
  tryCatch({
    p <- vis_gene_TIL_cor(
      Gene = gene_processed,
      cor_method = "spearman",
      data_type = "transcript",
      sig = c("B cell memory_CIBERSORT",
              "B cell naive_CIBERSORT",
              "B cell plasma_CIBERSORT",
              "Cancer associated fibroblast_XCELL",
              "Class-switched memory B cell_XCELL",
              "Common lymphoid progenitor_XCELL",
              "Endothelial cell_XCELL",
              "Eosinophil_CIBERSORT",
              "Granulocyte-monocyte progenitor_XCELL",
              "Hematopoietic stem cell_XCELL",
              "Macrophage M0_CIBERSORT",
              "Macrophage M1_CIBERSORT",
              "Macrophage M2_CIBERSORT",
              "Mast cell activated_CIBERSORT",
              "Monocyte_CIBERSORT",
              "Myeloid dendritic cell activated_CIBERSORT",
              "Myeloid dendritic cell resting_CIBERSORT",
              "Neutrophil_CIBERSORT",
              "NK cell activated_CIBERSORT",
              "NK cell resting_CIBERSORT",
              "T cell CD4+ memory activated_CIBERSORT",
              "T cell CD4+ memory resting_CIBERSORT",
              "T cell CD4+ naive_CIBERSORT",
              "T cell CD4+ Th1_XCELL",
              "T cell CD4+ Th2_XCELL",
              "T cell CD8+_CIBERSORT",
              "T cell follicular helper_CIBERSORT",
              "T cell gamma delta_CIBERSORT",
              "T cell regulatory (Tregs)_CIBERSORT"),
      Plot = "TRUE"
    )
    
    if (!is.null(p) && !is.null(p[["data"]])) {
      # Adding correlation data to the results data frame
      til_results <- bind_rows(til_results, p[["data"]])
    } else {
      warning("No data for gene: ", gene)
      null_genes <- c(null_genes, gene)
    }
  }, error = function(e) {
    cat("Error processing gene: ", gene, " with error message: ", e$message, "\n")
    null_genes <- c(null_genes, gene)
  })
}

# Checking if there are genes that failed and writing them to a file
if (length(null_genes) > 0) {
  null_genes_file_path <- file.path(getwd(), "null_genes_that_failed.txt")
  write.table(null_genes, file = null_genes_file_path, col.names = FALSE, row.names = FALSE, quote = FALSE)
  cat("List of genes with no data or errors written to ", null_genes_file_path, "\n")
}

# # Apply the Bonferroni-Holm correction to the p-values
# til_results <- til_results %>%
#   mutate(padj = p.adjust(p.value, method = "holm"))

# Filtering the data to include only rows where the FDR is less than or equal to 0.00000005
til_filtered <- til_results %>%
  filter(p.value <= 0.001)

# Exporting the filtered results to a file
export(til_results, "transcritos/TIL_results_transcript.tsv")
export(til_filtered, "genes/TIL_filtered_results_transcript.tsv")

# Applying the function to clean gene names in the Target_genes_cancer data frame
Target <- Table %>%
  mutate(Transcript = sapply(Transcript, process_gene_name))

# Selecting Signatures and cancer types
Target_genes_cancer <- Target %>%
  distinct(Cancer_type, Transcript)

# Merging unique columns back into the original data frame
Target_genes_cancer <- Target_genes_cancer %>%
  left_join(Target %>% select(Transcript, Cancer_type, Expression), 
            by = c("Cancer_type", "Transcript"))


# Extracting data from the immune_cells, cor, and p.value columns from the til_filtered data frame 
# and merging it with the Target_genes data frame
data <- merge(
  Target_genes_cancer, 
  til_filtered %>% select(cancer, gene, immune_cells, cor, p.value), 
  by.x = c("Cancer_type", "Transcript"), 
  by.y = c("cancer", "gene"), 
  all.x = TRUE
)

# Exporting the filtered results to a file
export(data, "transcritos/TIL_filtered_results_merge_transcript.tsv")

cat("Analysis complete")
