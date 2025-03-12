# Load necessary libraries
library(tidyverse)
library(rio)

setwd("C:/Users/higor/OneDrive/Área de Trabalho/Anti_pro_tumoral_classification/miRNA")

# Load the data
data <- import("TIL_filtered_results_merge.tsv")

# Check the distribution of expressions
expression_dist <- data %>% 
  group_by(Expression) %>% 
  summarize(Count = n())

# Define the immune cell classification table
cell_classification <- tibble(
  Cell_Type = c("B cell memory_CIBERSORT", "B cell naive_CIBERSORT", "B cell plasma_CIBERSORT",
                "Cancer associated fibroblast_XCELL", "Class-switched memory B cell_XCELL",
                "Common lymphoid progenitor_XCELL", "Endothelial cell_XCELL", "Eosinophil_CIBERSORT",
                "Granulocyte-monocyte progenitor_XCELL", "Hematopoietic stem cell_XCELL", 
                "Macrophage M0_CIBERSORT", "Macrophage M1_CIBERSORT", "Macrophage M2_CIBERSORT", 
                "Mast cell activated_CIBERSORT", "Monocyte_CIBERSORT", "Myeloid dendritic cell activated_CIBERSORT",
                "Myeloid dendritic cell resting_CIBERSORT", "Neutrophil_CIBERSORT", 
                "NK cell activated_CIBERSORT", "NK cell resting_CIBERSORT", 
                "T cell CD4+ memory activated_CIBERSORT", "T cell CD4+ memory resting_CIBERSORT",
                "T cell CD4+ naive_CIBERSORT", "T cell CD4+ Th1_XCELL", "T cell CD4+ Th2_XCELL", 
                "T cell CD8+_CIBERSORT", "T cell follicular helper_CIBERSORT", 
                "T cell gamma delta_CIBERSORT", "T cell regulatory (Tregs)_CIBERSORT"),
  Classification = c("dual", "dual", "dual", "pro-tumoral", "dual", "dual",
                     "pro-tumoral", "dual", "pro-tumoral", "pro-tumoral", "dual", "anti-tumoral",
                     "pro-tumoral", "dual", "dual", "dual", "dual", "dual", 
                     "anti-tumoral", "dual", "anti-tumoral", "dual", "dual", "anti-tumoral", 
                     "pro-tumoral", "anti-tumoral", "dual", "dual", "pro-tumoral")
)

# Merge correlation data with immune cell classification
merged_data <- data %>% 
  left_join(cell_classification, by = c("immune_cells" = "Cell_Type"))

# Function to classify signatures based on correlations and return scoring details
classify_signature <- function(data) {
  if (nrow(data) == 0 || all(is.na(data$cor))) {
    return(list(Classification = "NS", Score_Details = "No data"))
  }
  
  anti_tumoral_score <- 0
  pro_tumoral_score <- 0
  dualista_score <- 0
  
  anti_tumoral_cells <- c()
  pro_tumoral_cells <- c()
  dualista_cells <- c()
  
  for (i in 1:nrow(data)) {
    row <- data[i, ]
    if (is.na(row$cor) || is.na(row$Expression) || is.na(row$Classification)) {
      next
    }
    cell_info <- paste(row$immune_cells, "(", row$cor, ")", sep = "")
    correlation_value <- row$cor
    
    if (row$Expression == "Overexpression") {
      if (row$Classification == "anti-tumoral") {
        anti_tumoral_score <- anti_tumoral_score + correlation_value
        anti_tumoral_cells <- c(anti_tumoral_cells, cell_info)
      } else if (row$Classification == "pró-tumoral") {
        pro_tumoral_score <- pro_tumoral_score + correlation_value
        pro_tumoral_cells <- c(pro_tumoral_cells, cell_info)
      } else if (row$Classification == "dualista") {
        dualista_score <- dualista_score + correlation_value
        dualista_cells <- c(dualista_cells, cell_info)
      }
    } else if (row$Expression == "Underexpression") {
      if (row$Classification == "anti-tumoral") {
        anti_tumoral_score <- anti_tumoral_score - correlation_value
        anti_tumoral_cells <- c(anti_tumoral_cells, cell_info)
      } else if (row$Classification == "pró-tumoral") {
        pro_tumoral_score <- pro_tumoral_score - correlation_value
        pro_tumoral_cells <- c(pro_tumoral_cells, cell_info)
      } else if (row$Classification == "dualista") {
        dualista_score <- dualista_score - correlation_value
        dualista_cells <- c(dualista_cells, cell_info)
      }
    } else if (row$Expression == "Unchanged") {
      if (row$Classification == "anti-tumoral") {
        anti_tumoral_score <- anti_tumoral_score + correlation_value
        anti_tumoral_cells <- c(anti_tumoral_cells, cell_info)
      } else if (row$Classification == "pró-tumoral") {
        pro_tumoral_score <- pro_tumoral_score + correlation_value
        pro_tumoral_cells <- c(pro_tumoral_cells, cell_info)
      } else if (row$Classification == "dualista") {
        dualista_score <- dualista_score + correlation_value
        dualista_cells <- c(dualista_cells, cell_info)
      }
    }
  }
  
  classification <- if (anti_tumoral_score > pro_tumoral_score & anti_tumoral_score > dualista_score) {
    "anti-tumoral"
  } else if (pro_tumoral_score > anti_tumoral_score & pro_tumoral_score > dualista_score) {
    "pró-tumoral"
  } else {
    "dualista"
  }
  
  score_details <- paste("Anti-tumoral:", anti_tumoral_score, 
                         "(", paste(anti_tumoral_cells, collapse = ", "), ")",
                         "Pró-tumoral:", pro_tumoral_score, 
                         "(", paste(pro_tumoral_cells, collapse = ", "), ")",
                         "Dualista:", dualista_score, 
                         "(", paste(dualista_cells, collapse = ", "), ")")
  
  return(list(Classification = classification, Score_Details = score_details))
}

# Apply the classification function to data grouped by signature and cancer type
classification_result <- merged_data %>%
  group_by(miRNA_signature, Cancer_type) %>%
  summarize(Result = list(classify_signature(pick(everything()))), .groups = 'drop') %>%
  mutate(Classification = map_chr(Result, "Classification"),
         Score_Details = map_chr(Result, "Score_Details")) %>%
  select(-Result)

# Add the Expression column
classification_result <- classification_result %>%
  left_join(data %>% select(Cancer_type, miRNA_signature, Expression) %>%
              distinct(), by = c("Cancer_type", "miRNA_signature"))

# Export data
export(cell_classification, "cell_classification.tsv")
export(merged_data, "merged_data.tsv")
export(classification_result, "classification_result_miRNA.tsv")
