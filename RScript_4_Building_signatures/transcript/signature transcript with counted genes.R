# Load necessary packages
library(dplyr)
library(tidyr)
library(rio)

# Set the working directory
setwd("C:/Users/Emamnuell/Desktop/panoptosis/genes/genes_24_optosis/codigos/4°_signature/transcript")

# Load the data
signature_transcript <- import("signature_transcript.tsv")
gene_info_transcritos <- import("gene_info_trancritos.xlsx")

# Remove extra spaces in transcript IDs
signature_transcript$Transcript <- trimws(signature_transcript$Transcript)
gene_info_transcritos$Transcript_ID <- trimws(gene_info_transcritos$Transcript_ID)

# Function to map transcripts to genes and count unique transcripts
map_transcripts_to_genes <- function(transcript_string, gene_info) {
  # Display the signature being processed
  cat("Processing signature:", transcript_string, "\n")
  
  # Split the transcripts
  transcripts <- unlist(strsplit(gsub("[()]", "", transcript_string), " \\+ "))
  
  # Map each transcript to its corresponding gene
  gene_list <- sapply(transcripts, function(transcript) {
    gene <- gene_info %>% 
      filter(Transcript_ID == transcript) %>% 
      pull(Display_Name)
    
    # Check if the transcript was found
    if(length(gene) == 0) {
      cat("Transcript not found:", transcript, "\n")
      return("N/A")
    } else {
      cat("Transcript found:", transcript, "->", gene, "\n")
      return(gene)
    }
  })
  
  # Count genes and format the output
  gene_table <- table(gene_list)
  gene_counts <- sapply(names(gene_table), function(gene) {
    total_transcripts <- nrow(gene_info %>% filter(Display_Name == gene))
    paste0("(", gene, "(", gene_table[gene], "/", total_transcripts, ")", ")")
  })
  
  # Join the genes back into a string
  paste(gene_counts, collapse = " + ")
}

# Apply the function to each row of the signature table
signature_transcript$member_gene <- sapply(signature_transcript$Transcript, 
                                           map_transcripts_to_genes, 
                                           gene_info = gene_info_transcritos)

# Check the result
print(signature_transcript)

# Save the result in a new file
export(signature_transcript, "signature_transcript_with_counted_genes.tsv")
