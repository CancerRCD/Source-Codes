library(rio)
library(dplyr)
setwd("C:/Users/Emamnuell/Desktop/panoptosis/genes/genes_24_optosis/codigos/4°_signature/transcript")

# Import the transcript signature data
transcript_stemness <- import("signature_transcript.tsv")

# Filter rows where the "Phenotype" column contains the word "Stemness"
transcript_stemness_filtered <- transcript_stemness %>%
  filter(grepl("Stemness", Phenotype))

# Filter rows where the "Phenotype" column contains the word "MSI"
transcript_MSI_filtered <- transcript_stemness %>%
  filter(grepl("MSI", Phenotype))

# Filter rows where the "Phenotype" column contains the word "TMB"
transcript_TMB_filtered <- transcript_stemness %>%
  filter(grepl("TMB", Phenotype))

# Create data frames for each filtered transcript set
signature_Stemness <- data.frame(transcript_stemness_filtered$Transcript)
signature_MSI <- data.frame(transcript_MSI_filtered$Transcript)
signature_TMB <- data.frame(transcript_TMB_filtered$Transcript)

# Rename columns to "Genes"
colnames(signature_Stemness) <- "Genes"
colnames(signature_MSI) <- "Genes"
colnames(signature_TMB) <- "Genes"

# Set the working directory and export the filtered gene signatures
setwd("C:/Users/Emamnuell/Desktop/panoptosis/genes/genes_24_optosis/codigos/5°_Radar_siganature/Stemness_PanCan/Transcript_vs_Stemness/data")
export(signature_Stemness, "Target_genes.txt")

setwd("C:/Users/Emamnuell/Desktop/panoptosis/genes/genes_24_optosis/codigos/5°_Radar_siganature/MSI_PanCan/Transcript_vs_MSI/data")
export(signature_MSI, "Target_genes.txt")

setwd("C:/Users/Emamnuell/Desktop/panoptosis/genes/genes_24_optosis/codigos/5°_Radar_siganature/TMB_PanCan/Transcript_vs_TMB/data")
export(signature_TMB, "Target_genes.txt")
