library(rio)
library(dplyr)

setwd("C:/Users/higor/OneDrive/Área de Trabalho/TIL_results")

##### mRNA #####

data1 <- import("C:/Users/higor/OneDrive/Área de Trabalho/Anti_pro_tumoral_classification/mRNA/classification_result_mRNA.tsv")
data2 <- import("C:/Users/higor/OneDrive/Área de Trabalho/Immune_phenotype_classification/mRNA/TIL_classification.tsv")

# Add the Expression column
TIL_results_mRNA <- data1 %>%
  left_join(data2 %>% select(Cancer_type, Gene_signature, classification, score_details) %>%
              distinct(), by = c("Cancer_type", "Gene_signature"))

# Reorder the columns
New_column_order <- c("Gene_signature", "Cancer_type", "Expression", "classification", "score_details", "Classification", "Score_Details")

TIL_results_mRNA <- TIL_results_mRNA[New_column_order]

export(TIL_results_mRNA, "TIL_results_mRNA.tsv")

##### miRNA #####

data3 <- import("C:/Users/higor/OneDrive/Área de Trabalho/Anti_pro_tumoral_classification/miRNA/classification_result_miRNA.tsv")
data4 <- import("C:/Users/higor/OneDrive/Área de Trabalho/Immune_phenotype_classification/miRNA/TIL_classification.tsv")

# Add the Expression column
TIL_results_miRNA <- data3 %>%
  left_join(data4 %>% select(Cancer_type, miRNA_signature, classification, score_details) %>%
              distinct(), by = c("Cancer_type", "miRNA_signature"))

# Reorder the columns
New_column_order <- c("miRNA_signature", "Cancer_type", "Expression", "classification", "score_details", "Classification", "Score_Details")

TIL_results_miRNA <- TIL_results_miRNA[New_column_order]

export(TIL_results_miRNA, "TIL_results_miRNA.tsv")

##### Transcript #####

data5 <- import("C:/Users/higor/OneDrive/Área de Trabalho/Anti_pro_tumoral_classification/Transcript/classification_result_transcript.tsv")
data6 <- import("C:/Users/higor/OneDrive/Área de Trabalho/Immune_phenotype_classification/Transcript/TIL_classification.tsv")

# Add the Expression column
TIL_results_transcript <- data5 %>%
  left_join(data6 %>% select(Cancer_type, Transcript, classification, score_details) %>%
              distinct(), by = c("Cancer_type", "Transcript"))

# Reorder the columns
New_column_order <- c("Transcript", "Cancer_type", "Expression", "classification", "score_details", "Classification", "Score_Details")

TIL_results_transcript <- TIL_results_transcript[New_column_order]

export(TIL_results_transcript, "TIL_results_transcript.tsv")
