#### RSCRIPT TO CLASSIFY THE TUMOR IMMUNE PHENOTYPE ADAPTIVELY

# Load necessary packages
library(dplyr)
library(tidyr)
library(purrr)
library(rio)

setwd("C:/Users/higor/OneDrive/Área de Trabalho/Immune_phenotype_classification/mRNA")

# Import data
data <- import("TIL_filtered_results_merge.tsv")

# Filter rows with immune cells and correlation present
data <- data %>%
  filter(!is.na(immune_cells) & !is.na(cor))

# Calculate descriptive statistics for each immune cell
stats <- data %>%
  group_by(immune_cells) %>%
  summarize(
    Q1 = quantile(cor, 0.25, na.rm = TRUE),
    Q2 = quantile(cor, 0.5, na.rm = TRUE),  # Median
    Q3 = quantile(cor, 0.75, na.rm = TRUE)
  )

# Function to classify based on correlations and gene expression
classify_tumor <- function(cor_data, stats) {
  hot_criteria <- c("T cell CD8+_CIBERSORT", "NK cell activated_CIBERSORT", "Macrophage M1_CIBERSORT")
  cold_criteria <- c("T cell CD8+_CIBERSORT", "NK cell activated_CIBERSORT", "Macrophage M1_CIBERSORT", 
                     "Macrophage M2_CIBERSORT", "T cell regulatory (Tregs)_CIBERSORT")
  variable_criteria <- c("T cell CD8+_CIBERSORT", "NK cell activated_CIBERSORT", "Macrophage M1_CIBERSORT", 
                         "Macrophage M2_CIBERSORT", "T cell regulatory (Tregs)_CIBERSORT")
  
  # Weights for critical cells 
  cell_weights <- list("T cell CD8+_CIBERSORT" = 3, "NK cell activated_CIBERSORT" = 3, 
                       "Macrophage M1_CIBERSORT" = 1, "Macrophage M2_CIBERSORT" = 1, "T cell regulatory (Tregs)_CIBERSORT" = 1)
  
  # Weighted counters for each classification
  hot_score <- 0
  cold_score <- 0
  variable_score <- 0
  expression_type <- unique(cor_data$Expression)
  
  # Check if there are immune cells in the signature
  if (nrow(cor_data) == 0) {
    return(list(classification = "NS", score_details = data.frame()))
  }
  
  # Table to store weight details
  score_details <- data.frame(immune_cells = character(), rho = numeric(), weight = numeric(), score_type = character(), stringsAsFactors = FALSE)
  
  # Analyze each immune cell
  for (i in 1:nrow(cor_data)) {
    row <- cor_data[i, ]
    cell_type <- row$immune_cells
    rho <- row$cor
    weight <- cell_weights[[cell_type]]
    
    cell_stats <- stats %>% filter(immune_cells == cell_type)
    Q1 <- cell_stats$Q1
    Q2 <- cell_stats$Q2
    Q3 <- cell_stats$Q3
    
    if (expression_type == "Overexpression") {
      if (cell_type %in% hot_criteria && rho >= Q3) {
        hot_score <- hot_score + weight
        score_details <- rbind(score_details, data.frame(immune_cells = cell_type, rho = rho, weight = weight, score_type = "Hot"))
      }
      if (cell_type %in% cold_criteria && rho <= Q1) {
        cold_score <- cold_score + weight
        score_details <- rbind(score_details, data.frame(immune_cells = cell_type, rho = rho, weight = weight, score_type = "Cold"))
      }
      if (cell_type %in% variable_criteria && rho > Q1 && rho < Q3) {
        variable_score <- variable_score + weight
        score_details <- rbind(score_details, data.frame(immune_cells = cell_type, rho = rho, weight = weight, score_type = "Variable"))
      }
    } else if (expression_type == "Underexpression") {
      if (cell_type %in% hot_criteria && rho <= Q1) {
        hot_score <- hot_score + weight
        score_details <- rbind(score_details, data.frame(immune_cells = cell_type, rho = rho, weight = weight, score_type = "Hot"))
      }
      if (cell_type %in% cold_criteria && rho >= Q3) {
        cold_score <- cold_score + weight
        score_details <- rbind(score_details, data.frame(immune_cells = cell_type, rho = rho, weight = weight, score_type = "Cold"))
      }
      if (cell_type %in% variable_criteria && rho > Q1 && rho < Q3) {
        variable_score <- variable_score + weight
        score_details <- rbind(score_details, data.frame(immune_cells = cell_type, rho = rho, weight = weight, score_type = "Variable"))
      }
    } else if (expression_type == "Unchanged") {
      if (cell_type %in% hot_criteria && rho >= Q3) {
        hot_score <- hot_score + weight
        score_details <- rbind(score_details, data.frame(immune_cells = cell_type, rho = rho, weight = weight, score_type = "Hot"))
      }
      if (cell_type %in% cold_criteria && rho <= Q1) {
        cold_score <- cold_score + weight
        score_details <- rbind(score_details, data.frame(immune_cells = cell_type, rho = rho, weight = weight, score_type = "Cold"))
      }
      if (cell_type %in% variable_criteria && rho > Q1 && rho < Q3) {
        variable_score <- variable_score + weight
        score_details <- rbind(score_details, data.frame(immune_cells = cell_type, rho = rho, weight = weight, score_type = "Variable"))
      }
    }
  }
  
  # Classify based on weighted counters
  if (hot_score == cold_score && hot_score > 0) {
    classification <- "Variable"
  } else if (hot_score == variable_score && hot_score > 0) {
    classification <- "Hot"
  } else if (cold_score == variable_score && cold_score > 0) {
    classification <- "Cold"
  } else {
    classification <- if (hot_score > cold_score && hot_score > variable_score) {
      "Hot"
    } else if (cold_score > hot_score && cold_score > variable_score) {
      "Cold"
    } else if (variable_score > 0) {
      "Variable"
    } else {
      "NS"
    }
  }
  
  return(list(classification = classification, score_details = score_details))
}

# Group data by gene and cancer and apply classification
result <- data %>%
  group_by(Gene_signature, Cancer_type) %>%
  nest() %>%
  mutate(classification_result = map(data, ~classify_tumor(.x, stats))) %>%
  mutate(classification = map_chr(classification_result, "classification"),
         score_details = map(classification_result, "score_details")) %>%
  select(-data, -classification_result)

# Select columns of interest
columns_of_interest <- c("Cancer_type", "Gene_signature", "classification")
selected_results <- result[, columns_of_interest, drop = FALSE]

# Export to a TSV file
export(selected_results, "TIL_classification_selecionado.tsv")

# Convert the score_details column from list to character
result$score_details <- sapply(result$score_details, toString)

# Save the results in a CSV file
export(result, "TIL_classification.tsv")
