# Load necessary libraries
library(survival)
library(UCSCXenaShiny)
library(UCSCXenaTools)
library(dplyr)
library(rio)
library(readr)
# install.packages("rio", dependencies = TRUE)

# Set working directory and import the table
setwd("C:/Users/Emamnuell/Desktop/panoptosis/genes/genes_24_optosis/codigos/7°_survival")
signature_gene <- import("all_signature.tsv")

# Convert the Genotype column to lowercase only if it is not "mRNA" or "miRNA"
signature_gene$Genotype <- ifelse(
  signature_gene$Genotype %in% c("mRNA", "miRNA"),
  signature_gene$Genotype,          # Keep values as they are
  tolower(signature_gene$Genotype)  # Convert to lowercase
)

# Function to process simple and compound gene names
process_gene_name <- function(gene) {
  if (is.na(gene) || gene == "") {
    return(NA_character_)  # Return NA if the gene is NA or empty
  }
  
  cat("Processing gene name:", gene, "\n")  # Debugging message
  
  # Check if the gene contains a plus (+) sign
  if (grepl("\\+", gene)) {
    return(gene)  # Keep the gene as it is if there is a plus sign
  } else {
    cleaned_gene <- gsub("`|\\(|\\)", "", gene)  # Remove backticks and parentheses if there is no plus sign
    cat("Cleaned gene name:", cleaned_gene, "\n")  # Debugging message
    return(cleaned_gene)
  }
}


run_analysis <- function(gene_ID, gene_var, tumor) {
  gene_ID <- process_gene_name(gene_ID)  # Process the gene name
  
  Check if gene_ID is valid
  if (is.na(gene_ID) || gene_ID == "") {
    warning(paste("Invalid gene ID:", gene_ID, "tumor:", tumor, "gene_var:", gene_var))
    return(data.frame(
      Gene = gene_ID,
      tumor = tumor,
      gene_var = gene_var,
      status = "DFI",
      p_val = NA,
      log_rank = NA,
      worst_prognosis_group = "NS"
    ))
  }
  
  opt_pancan = .opt_pancan
  opt_pancan$toil_cnv$use_thresholded_data = TRUE
  
  # Retrieve clinical data
  data_surv <- tryCatch({
    data_surv <- tcga_surv_get(
      item = gene_ID,
      TCGA_cohort = tumor,
      profile = gene_var,
      TCGA_cli_data = dplyr::full_join(load_data("tcga_clinical"), load_data("tcga_surv"), by = "sample"),
      opt_pancan = opt_pancan
    )
  }, error = function(e) {
    warning(paste("Failed to get data for gene:", gene_ID, "tumor:", tumor, "gene_var:", gene_var, ":", e$message))
    return(NULL)
  })
  
  # Check if data_surv is not empty or null
  if (is.null(data_surv) || nrow(data_surv) == 0) {
    warning(paste("No data found for gene:", gene_ID, "tumor:", tumor, "gene_var:", gene_var))
    return(data.frame(
      Gene = gene_ID,
      tumor = tumor,
      gene_var = gene_var,
      status = "DFI",
      p_val = NA,
      log_rank = NA,
      worst_prognosis_group = "NS"
    ))
  }
  
  # Sort and add percentile rank for the variable of interest (gene_var)
  data_surv <- data_surv %>%
    arrange(.data$value) %>%
    mutate(per_rank = 100 / nrow(.) * (1:nrow(.)))
  
  # Analyze different types of data
  if (gene_var == "cnv") {
    # Case for 'cnv'
    data_surv <- data_surv %>%
      mutate(group = case_when(
        .data$value == 0 ~ "Normal",
        .data$value > 0 ~ "Duplicated",
        .data$value < 0 ~ "Deleted",
        TRUE ~ NA_character_
      ))
    
  } else if (gene_var == "mutation") {
    # Case for 'mutation'
    data_surv <- data_surv %>%
      mutate(group = case_when(
        .data$value == "0" ~ "WT",  # Tipo selvagem
        .data$value != "0" ~ "MT",  # Tipo mutante
        TRUE ~ NA_character_
      ))
    
  } else {
    # Define cutoff points using quantiles
    cutpoint <- quantile(data_surv$per_rank, probs = c(0.5)) # Median as cutoff point
    
    data_surv <- data_surv %>%
      mutate(group = case_when(
        .data$per_rank > cutpoint ~ "High",  # Values above the 50th percentile
        .data$per_rank <= cutpoint ~ "Low",  # Values up to the 50th percentile
        TRUE ~ NA_character_
      ))
  }
  
  # Filter data to remove NA in 'group'
  data_surv <- data_surv %>%
    filter(!is.na(group))
  
  # Check if there is enough data after filtering
  if (nrow(data_surv) == 0) {
    warning(paste("Insufficient data after filtering for gene:", gene_ID, "tumor:", tumor, "gene_var:", gene_var))
    return(data.frame(
      Gene = gene_ID,
      tumor = tumor,
      gene_var = gene_var,
      status = "DFI",
      p_val = NA,
      log_rank = NA,
      worst_prognosis_group = "NS"
    ))
  }
  
  plot_surv <- tryCatch({
    # Check if gene_var is different from "mutation" and "cnv"
    if (!(gene_var %in% c("mutation"))) {
      # Set the cutoff mode and custom cutoff point
      tcga_surv_plot(
        data_surv,
        time = "DFI.time",
        status = "DFI",
        cutoff_mode = "Custom",             # Use the custom cutoff mode
        profile = gene_var,
        palette = "aaas",
        custom_cutpoint = if (gene_var != "cnv") cutpoint else NULL           # Add the custom cutoff point
      )
    } else {
      # Set default configuration when gene_var is "mutation" or "cnv"
      tcga_surv_plot(
        data_surv,
        time = "DFI.time",
        status = "DFI",
        profile = gene_var,
        palette = "aaas"
      )
    }
  }, error = function(e) {
    # Capture and display error messages
    message("Erro ao gerar o gráfico de sobrevivência: ", e$message)
    NULL
  })
  
  
  # Verify if survival data was generated correctly
  if (is.null(plot_surv)) {
    return(data.frame(
      Gene = gene_ID,
      tumor = tumor,
      gene_var = gene_var,
      status = "DFI",
      p_val = NA,
      log_rank = NA,
      worst_prognosis_group = "NS"
    ))
  }
  
  # Extract survival data from the plot_surv object
  surv_data <- plot_surv[["plot"]][["data"]]
  
  # Verify if survival data was generated correctly
  if (is.null(surv_data)) {
    warning(paste("Failed to generate survival data for gene:", gene_ID, "tumor:", tumor, "gene_var:", gene_var))
    return(data.frame(
      Gene = gene_ID,
      tumor = tumor,
      gene_var = gene_var,
      status = "DFI",
      p_val = NA,
      log_rank = NA,
      worst_prognosis_group = "NS"
    ))
  }
  
  # Convert the list to a data frame
  surv_df <- as.data.frame(surv_data)
  
  # Check if the 'group' column exists
  if (!"group" %in% colnames(data_surv)) {
    warning(paste("Column 'group' not found for gene:", gene_ID, "tumor:", tumor, "gene_var:", gene_var))
    return(data.frame(
      Gene = gene_ID,
      tumor = tumor,
      gene_var = gene_var,
      status = "DFI",
      p_val = NA,
      log_rank = NA,
      worst_prognosis_group = "NS"
    ))
  }
  
  # Calculate survival differences using the log-rank test
  surv_diff <- tryCatch({
    survdiff(Surv(DFI.time, DFI) ~ group, data = data_surv)
  }, error = function(e) {
    warning(paste("Failed to calculate survival differences for gene:", gene_ID, "tumor:", tumor, "gene_var:", gene_var, ":", e$message))
    return(NULL)
  })
  
  # Check if surv_diff is null
  if (is.null(surv_diff)) {
    return(data.frame(
      Gene = gene_ID,
      tumor = tumor,
      gene_var = gene_var,
      status = "DFI",
      p_val = NA,
      log_rank = NA,
      worst_prognosis_group = "NS"
    ))
  }
  
  p_val <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
  chisq_log_rank <- surv_diff$chisq
  
  # Check when each group reaches a 50% survival probability
  group_surv_50 <- tryCatch({
    surv_df %>%
      group_by(group) %>%
      summarize(time_50 = min(time[surv <= 0.5], na.rm = TRUE))
  }, error = function(e) {
    warning(paste("Failed to summarize survival data for gene:", gene_ID, "tumor:", tumor, "gene_var:", gene_var, ":", e$message))
    return(NULL)
  })
  
  # Print intermediate values
  print(group_surv_50)
  
  # Check if group_surv_50 is null
  if (is.null(group_surv_50)) {
    return(data.frame(
      Gene = gene_ID,
      tumor = tumor,
      gene_var = gene_var,
      status = "DFI",
      p_val = NA,
      log_rank = NA,
      worst_prognosis_group = "NS"
    ))
  }
  
  # Auxiliary function for performing pairwise comparisons
  compare_groups <- function(data, group1, group2) {
    subset_data <- data %>% filter(group %in% c(group1, group2))
    if (nrow(subset_data) == 0 || length(unique(subset_data$group)) < 2) {
      return(NA)
    }
    surv_diff_pair <- tryCatch({
      survdiff(Surv(DFI.time, DFI) ~ group, data = subset_data)
    }, error = function(e) {
      warning(paste("Failed to calculate paired survival differences for groups:", group1, "vs", group2, ":", e$message))
      return(NA)
    })
    if (is.null(surv_diff_pair) || !is.list(surv_diff_pair)) {
      return(NA)
    }
    p_val_pair <- 1 - pchisq(surv_diff_pair$chisq, length(surv_diff_pair$n) - 1)
    return(p_val_pair)
  }
  
  # Perform pairwise comparisons
  p_val_deleted_vs_normal <- compare_groups(data_surv, "Deleted", "Normal")
  p_val_deleted_vs_duplicated <- compare_groups(data_surv, "Deleted", "Duplicated")
  p_val_normal_vs_duplicated <- compare_groups(data_surv, "Normal", "Duplicated")
  
  # Check if data exists for all three groups
  groups_present <- unique(data_surv$group)
  has_deleted <- "Deleted" %in% groups_present
  has_duplicated <- "Duplicated" %in% groups_present
  has_normal <- "Normal" %in% groups_present
  
  # Verify if p_val is greater than 0.05
  if (p_val > 0.05) {
    worst_prognosis_group <- "NS"
  } else {
    # Determine the group with the worst prognosis for CNV
    if (gene_var == "cnv") {
      if (has_deleted && has_duplicated && has_normal) {
        # Compare "Deleted", "Duplicated", and "Normal"
        if (!is.na(p_val_deleted_vs_normal) && p_val_deleted_vs_normal < 0.05 &&
            (is.na(p_val_deleted_vs_duplicated) || p_val_deleted_vs_duplicated > 0.05) &&
            (is.na(p_val_normal_vs_duplicated) || p_val_normal_vs_duplicated > 0.05)) {
          worst_prognosis_group <- "Deleted"
        } else if (!is.na(p_val_normal_vs_duplicated) && p_val_normal_vs_duplicated < 0.05 &&
                   (is.na(p_val_deleted_vs_duplicated) || p_val_deleted_vs_duplicated > 0.05) &&
                   (is.na(p_val_deleted_vs_normal) || p_val_deleted_vs_normal > 0.05)) {
          worst_prognosis_group <- "Duplicated"
        } else if (!is.na(p_val_deleted_vs_duplicated) && p_val_deleted_vs_duplicated < 0.05) {
          if (!is.na(p_val_normal_vs_duplicated) && p_val_normal_vs_duplicated < 0.05 &&
              p_val_normal_vs_duplicated < p_val_deleted_vs_normal) {
            worst_prognosis_group <- "Duplicated"
          } else {
            worst_prognosis_group <- "Deleted"
          }
        } else if (!is.na(p_val_deleted_vs_normal) && p_val_deleted_vs_normal < 0.05 &&
                   !is.na(p_val_normal_vs_duplicated) && p_val_normal_vs_duplicated < 0.05 &&
                   (is.na(p_val_deleted_vs_duplicated) || p_val_deleted_vs_duplicated > 0.05)) {
          worst_prognosis_group <- "Deleted/Duplicated"
        } else {
          worst_prognosis_group <- "NS"
        }
      } else if (has_deleted && has_normal && !has_duplicated) {
        # Compare "Normal" e "Deleted"
        worst_prognosis_group <- ifelse(
          group_surv_50$time_50[group_surv_50$group == "Normal"] < group_surv_50$time_50[group_surv_50$group == "Deleted"],
          "Normal", "Deleted"
        )
      } else if (has_duplicated && has_normal && !has_deleted) {
        # Compare "Normal" e "Duplicated"
        worst_prognosis_group <- ifelse(
          group_surv_50$time_50[group_surv_50$group == "Normal"] < group_surv_50$time_50[group_surv_50$group == "Duplicated"],
          "Normal", "Duplicated"
        )
      } else if (has_deleted && has_duplicated && !has_normal) {
        # Compare "Deleted" e "Duplicated"
        worst_prognosis_group <- ifelse(
          group_surv_50$time_50[group_surv_50$group == "Deleted"] < group_surv_50$time_50[group_surv_50$group == "Duplicated"],
          "Deleted", "Duplicated"
        )
      } else {
        worst_prognosis_group <- "NS"
      }
    } else {
      # Other comparisons depending on gene_var
      if (gene_var %in% c("mRNA", "transcript", "miRNA", "methylation", "protein")) {
        # VCheck if there is any surv value ≤ 0.5
        if (any(surv_df$surv <= 0.5)) {
          worst_prognosis_group <- ifelse(
            group_surv_50$time_50[group_surv_50$group == "Low"] < group_surv_50$time_50[group_surv_50$group == "High"],
            "Low", "High"
          )
        } else {
          # No surv values ≤ 0.5
          # Take the smallest surv value > 0.5
          min_surv_greater_than_0_5 <- surv_df %>%
            filter(surv > 0.5) %>%
            slice(which.min(surv))
          
          # Check if the p-value is less than 0.05
          if (p_val < 0.05) {
            worst_prognosis_group <- min_surv_greater_than_0_5$group
          } else {
            worst_prognosis_group <- "NS"
          }
        }
      } else if (gene_var == "mutation") {
        # Apply the same logic for mutation
        if (any(surv_df$surv <= 0.5)) {
          worst_prognosis_group <- ifelse(
            group_surv_50$time_50[group_surv_50$group == "MT"] < group_surv_50$time_50[group_surv_50$group == "WT"],
            "MT", "WT"
          )
        } else {
          # No surv values ≤ 0.5
          min_surv_greater_than_0_5 <- surv_df %>%
            filter(surv > 0.5) %>%
            slice(which.min(surv))
          
          # Check if the p-value is less than 0.05
          if (p_val < 0.05) {
            worst_prognosis_group <- min_surv_greater_than_0_5$group
          } else {
            worst_prognosis_group <- "NS"
          }
        }
      } else {
        # Default case
        worst_prognosis_group <- "NS"
      }
    }
  } 
  # Print worst prognosis values for debugging
  print(paste("Worst prognosis group: ", worst_prognosis_group))
  
  
  
  if (gene_var == "cnv") {
    return(data.frame(
      Gene = gene_ID,
      tumor = tumor,
      gene_var = gene_var,
      status = "DFI",
      p_val = p_val,
      log_rank = chisq_log_rank,
      worst_prognosis_group = worst_prognosis_group,
      p_val_deleted_vs_normal = p_val_deleted_vs_normal,
      p_val_deleted_vs_duplicated = p_val_deleted_vs_duplicated,
      p_val_normal_vs_duplicated = p_val_normal_vs_duplicated
    ))
  } else {
    return(data.frame(
      Gene = gene_ID,
      tumor = tumor,
      gene_var = gene_var,
      status = "DFI",
      p_val = p_val,
      log_rank = chisq_log_rank,
      worst_prognosis_group = worst_prognosis_group
    ))
  }
}

# Function to save progress to an RDS file
save_progress <- function(results, file_name) {
  saveRDS(results, file_name)
}

# Function to load progress from an RDS file
load_progress <- function(file_name) {
  if (file.exists(file_name)) {
    return(readRDS(file_name))
  } else {
    return(list(results_with_cnv = list(), results_without_cnv = list()))
  }
}

# File to save progress
checkpoint_file <- "progress_checkpoint_DFI.rds"

# Load saved progress
results_list <- load_progress(checkpoint_file)
results_with_cnv <- results_list$results_with_cnv
results_without_cnv <- results_list$results_without_cnv

# Check the start index based on saved progress
start_index <- length(results_with_cnv) + length(results_without_cnv) + 1

# Apply the function to each row of the signature_gene table
for (i in start_index:nrow(signature_gene)) {
  gene_ID <- signature_gene$signature[i]  # Use the signature column
  gene_var <- signature_gene$Genotype[i]
  tumor <- signature_gene$Cancer_type[i]
  
  # Check if any of the values are NA and skip the iteration if so
  if (is.na(gene_ID) || is.na(gene_var) || is.na(tumor)) {
    warning(paste("Skipping analysis for row", i, "due to NA values."))
    next
  }
  
  cat("Running analysis for gene_ID:", gene_ID, "gene_var:", gene_var, "tumor:", tumor, "\n")  # Debug message
  
  result <- run_analysis(gene_ID, gene_var, tumor)
  
  # Add the result to the correct list
  if (gene_var == "cnv") {
    results_with_cnv[[length(results_with_cnv) + 1]] <- result
  } else {
    results_without_cnv[[length(results_without_cnv) + 1]] <- result
  }
  
  # Save progress every 10 genes
  if (i %% 10 == 0) {
    save_progress(list(results_with_cnv = results_with_cnv, results_without_cnv = results_without_cnv), checkpoint_file)
  }
}

# Function to remove specified columns
remove_columns <- function(df) {
  columns_to_remove <- c("p_val_deleted_vs_normal", "p_val_deleted_vs_duplicated", "p_val_normal_vs_duplicated")
  
  # Check if columns exist and remove them
  df <- df[, !(names(df) %in% columns_to_remove), drop = FALSE]
  return(df)
}

# Apply the function to each element in the list
results_with_cnv_clean <- lapply(results_with_cnv, remove_columns)

# Combine the result lists into a single data frame
df_with_cnv <- do.call(rbind, results_with_cnv_clean)
df_without_cnv <- do.call(rbind, results_without_cnv)

# Remove extra columns from CNV results to match with the others
df_with_cnv_base <- df_with_cnv[, !(names(df_with_cnv) %in% c("p_val_deleted_vs_normal", "p_val_deleted_vs_duplicated", "p_val_normal_vs_duplicated"))]

# Combine the two data frames
results_df <- rbind(df_with_cnv_base, df_without_cnv)

# View the results
print(results_df)

# Export results to a TSV file
write.table(results_df, "survival_analysis_all_signature_DFI.tsv", sep="\t", row.names=FALSE, quote=FALSE)
