setwd("")

library(rentrez)
library(rio)

# Define your API key
api_key <- "YOUR API KEY"

genes <- import("search_results.tsv")

# Define the RIS file name to save intermediate results
ris_file <- "resultados_pubmed.ris"
rds_file <- "processed_pmids.rds"

# Define your PMIDs
pmids <- genes$PMID

# Check previously processed PMIDs
if (file.exists(rds_file)) {
  processed_pmids <- readRDS(rds_file)
  pmids <- pmids[!(pmids %in% processed_pmids)]
} else {
  processed_pmids <- c()
}

# Function to split vector into batches
split_into_batches <- function(x, batch_size) {
  split(x, ceiling(seq_along(x) / batch_size))
}

# Split PMIDs into smaller batches
pmid_batches <- split_into_batches(pmids, 100)  # Here we define 100 PMIDs per batch, you can adjust as needed

# Loop over each batch of PMIDs
for (i in seq_along(pmid_batches)) {
  # Print message indicating the start of the batch
  cat("Searching PMIDs", min(pmid_batches[[i]]), "to", max(pmid_batches[[i]]), "\n")
  
  # Use the entrez_fetch function to get article details
  # Set db = "pubmed" to search in PubMed
  # Set idtype = "pmid" to indicate we are using PMIDs
  # Set rettype = "ris" to return results in RIS format
  # Set api_key = api_key to use your API key
  # The results will be stored in a list called 'articles'
  articles <- entrez_fetch(db = "pubmed", id = pmid_batches[[i]], rettype = "ris", idtype = "pmid", api_key = api_key)
  
  # Write the results to the RIS file
  cat(articles, file = ris_file, append = TRUE)
  
  # Add processed PMIDs to the processed PMIDs vector
  processed_pmids <- c(processed_pmids, pmid_batches[[i]])
  
  # Save the processed PMIDs vector to an RDS file
  saveRDS(processed_pmids, file = rds_file)
  
  # Wait 1 second between calls to avoid excessive access in a short period
  if (i < length(pmid_batches)) {
    Sys.sleep(1)
  }
}

# Confirm completion
cat("All results have been saved in", ris_file, "\n")