# Load necessary libraries
library(rentrez)       # To access database APIs such as PubMed
library(tidyverse)     # For data manipulation and visualization
library(rio)           # For data import and export

# Set the working directory where files will be saved and read
setwd("")

# API key for use in PubMed
api_key <- "YOUR API KEY"

# Import data from an Excel file containing genes related to panoptosis
all_genes_five_optosis_com_vias <- import("genes_simplificado_model_optosis_excel.xlsx")
vias <- import("Integrated_Cell_Death_Definitions_cut_V2.xlsx")

my_genes <- all_genes_five_optosis_com_vias$name  # Extract only gene names

# Search terms used for PubMed queries
search_terms <- vias$`Cell Death Type`

retmax <- 3  # Maximum number of results returned per search

# Initialize variables to store results and track progress
search_results <- list()
batch_size <- 10
progress_file <- "progress.rds"
results_file <- "search_results.tsv"

# Check if a progress file exists and load previous progress
if (file.exists(progress_file)) {
  progress <- readRDS(progress_file)
  start_index <- progress$start_index
  current_gene <- progress$current_gene
} else {
  start_index <- 1
  current_gene <- 1
}

# Check if results were partially saved and load them
if (file.exists(results_file)) {
  df_results <- import(results_file)
  df_results <- df_results %>% mutate(across(everything(), as.character))
} else {
  df_results <- tibble(Gene_Term_Category = character(), Title = character(), PubDate = character(), DOI = character(), PMID = character())
}

# Main loop that processes genes in batches
for (i in seq(start_index, length(my_genes), batch_size)) {
  genes_batch <- my_genes[i:min(i + batch_size - 1, length(my_genes))]
  
  # Adjust the current gene index if resuming from a previous batch
  if (i == start_index && current_gene > 1) {
    genes_batch <- genes_batch[current_gene:length(genes_batch)]
  } else {
    current_gene <- 1
  }
  
  # Inner loop that processes each gene in the current batch
  for (gene_index in seq_along(genes_batch)) {
    gene <- as.character(genes_batch[gene_index])
    cat("Fetching information for gene:", gene, "\n")
    
    # Search each associated term in PubMed
    for (term in search_terms) {
      for (category in c("cancer", "tumor")) {
        query <- paste(gene, "AND", term, "-related genes AND", category, "[ALL] AND Homo sapiens[Organism] NOT Neanderthal NOT Denisova AND Free full text[Filter]")
        
        # Repeat the search until a valid response is obtained
        repeat {
          res <- tryCatch({
            entrez_search(db = "pubmed", term = query, retmax = retmax, api_key = api_key)
          }, error = function(e) {
            e
          })
          
          if (!inherits(res, "error")) {
            break
          } else {
            Sys.sleep(20)  # Pause to avoid API overload
          }
        }
        
        # Process article summaries if results are found
        if (res$count > 0) {
          summaries <- entrez_summary(db = "pubmed", id = res$ids, api_key = api_key)
          articles <- lapply(summaries, function(x) {
            if (!is.null(x)) {
              Title <- ifelse("title" %in% names(x), as.character(x$title), NA)
              PubDate <- ifelse("pubdate" %in% names(x), as.character(x$pubdate), NA)
              DOI <- ifelse("elocationid" %in% names(x), as.character(x$elocationid), NA)
              PMID <- ifelse("uid" %in% names(x), as.character(x$uid), NA)
              list(Title = Title, PubDate = PubDate, DOI = DOI, PMID = PMID)
            }
          })
          search_results[[paste(gene, term, category)]] <- articles
          
          # Append new results to the results DataFrame
          new_results <- bind_rows(lapply(articles, as_tibble), .id = "Gene_Term_Category")
          new_results <- new_results %>% mutate(across(everything(), as.character))
          new_results$Gene_Term_Category <- paste(gene, term, category)
          df_results <- bind_rows(df_results, new_results)
          
          # Export accumulated results
          export(df_results, results_file)
        }
      }
    }
    
    # Update progress after each gene
    current_gene <- gene_index + 1
    progress <- list(start_index = i, current_gene = current_gene)
    saveRDS(progress, progress_file)
  }
  
  # Pause between gene batches to avoid API overload
  Sys.sleep(20)
}

# Display and save final results
print(df_results)
export(df_results, "search_results.tsv")
