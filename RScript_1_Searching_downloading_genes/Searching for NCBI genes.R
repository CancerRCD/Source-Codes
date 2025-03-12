# Load the necessary packages
library(rentrez)
library(rio)

# Set the working directory
setwd("")

# Import gene names
genes <- import("genes_simplificado_model_optosis_excel.xlsx")

# Define gene symbols
gene_symbols <- genes$name  # Replace with your genes of interest

# Define your API key
api_key <- "YOUR API KEY"

# Function to fetch gene information
get_gene_info <- function(gene_symbol) {
  # Message indicating which gene is being processed
  cat("Fetching information for gene:", gene_symbol, "\n")
  
  # Search for Gene ID using the gene symbol and pass the API key
  search_results <- entrez_search(db = "gene", term = paste(gene_symbol, "[Gene Name] AND Homo sapiens[Organism]"), api_key = api_key)
  
  # Check if any results were found
  if (search_results$count == 0) {
    return(data.frame(gene_symbol = gene_symbol, uid = NA, name = NA, maplocation = NA, otheraliases = NA, 
                      chrstart = NA, chrstop = NA))
  }
  
  gene_id <- search_results$ids[1]  # Select the first returned ID
  
  # Retrieve gene details and pass the API key
  gene_summary <- entrez_summary(db = "gene", id = gene_id, api_key = api_key)
  
  # Extract genomic location information
  genomic_info <- gene_summary$genomicinfo
  
  # Extract necessary information
  chrstart <- ifelse("chrstart" %in% names(genomic_info), genomic_info[["chrstart"]], NA)
  chrstop <- ifelse("chrstop" %in% names(genomic_info), genomic_info[["chrstop"]], NA)
  
  data.frame(
    uid = gene_id,
    name = gene_summary$name,
    otheraliases = gene_summary$otheraliases,
    chrstart = chrstart,
    chrstop = chrstop,
    maplocation = gene_summary$maplocation
  )
}

# Load or create an .rds file to store progress
if (file.exists("progresso_genes.rds")) {
  progresso <- readRDS("progresso_genes.rds")
  gene_info_list <- readRDS("gene_info_list.rds")
} else {
  progresso <- 1
  gene_info_list <- list()
}

# Apply the function to all gene symbols from the saved progress
for (i in progresso:length(gene_symbols)) {
  try({
    gene_info <- get_gene_info(gene_symbols[i])
    gene_info_list[[i]] <- gene_info
    
    # Save progress
    saveRDS(i, "progresso_genes.rds")
    saveRDS(gene_info_list, "gene_info_list.rds")
  }, silent = TRUE)
}

# Ensure all entries have the same column structure
columns <- c("uid", "name","otheraliases", "chrstart", "chrstop", "maplocation")
gene_info_list <- lapply(gene_info_list, function(df) {
  df[setdiff(columns, names(df))] <- NA
  return(df[columns])
})

# Combine results into a data frame
resultados <- do.call(rbind, gene_info_list)
write.csv(resultados, "gene_info_results.tsv", row.names = FALSE)
export(resultados, "gene_info_results.xlsx")
