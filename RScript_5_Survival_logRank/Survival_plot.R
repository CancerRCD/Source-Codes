library(survival)
library(survminer)
library(rio)
library(survival)
library(survminer)
library(ggplot2)
library(cowplot)
library(dplyr)
library(cowplot)
library(UCSCXenaShiny)

setwd("C:/Users/Emamnuell/Desktop/UENF/laboratorio/panoptosis/genes/genes_24_optosis/codigos/plots/HNSC-1855.4.3.P.3.71.64.1.1.4")

# Define the desired profile
profile <- "cnv"  # Make sure this value is correct
metric <- "PFI"  # Change to "DSS", "DFI" or "PFI" as needed
gene <- "(CXCL10 + TNFRSF4)"
tumor <- "KIRP"

# Load and filter data only for male gender
data <- tcga_surv_get(
  item = gene,
  TCGA_cohort = tumor,
  profile = profile,  # Now using the profile variable here
  TCGA_cli_data = dplyr::full_join(load_data("tcga_clinical"), 
                                   load_data("tcga_surv"), by = "sample"),
  opt_pancan = .opt_pancan
)

# Apply categorization based on profile type
if (profile == "mutation") {
  data$value_group <- ifelse(data$value == 0, "WT", "MT")
} else if (profile %in% c("mRNA", "miRNA", "protein", "transcript", "methylation")) {
  data$value_group <- ifelse(data$value > median(data$value, na.rm = TRUE), "High", "Low")
} else if (profile == "cnv") {
  data$value_group <- ifelse(data$value == 0, "Normal",
                             ifelse(data$value < 0, "Deleted", "Duplicated"))
} else {
  stop("Unknown profile type: ", profile)
}

# Create the survival object based on the selected metric
surv_object <- switch(metric,
                      "OS" = Surv(time = data$OS.time, event = data$OS),
                      "DSS" = Surv(time = data$DSS.time, event = data$DSS),
                      "DFI" = Surv(time = data$DFI.time, event = data$DFI),
                      "PFI" = Surv(time = data$PFI.time, event = data$PFI),
                      stop("Invalid metric! Choose between 'OS', 'DSS', 'DFI', or 'PFI'.")
)


# Fit the Kaplan-Meier model
fit <- survfit(surv_object ~ value_group, data = data)

# Compute the log-rank test and extract the p-value
surv_diff <- tryCatch({
  survdiff(surv_object ~ value_group, data = data)
}, error = function(e) {
  warning("Error computing the log-rank test:", e$message)
  return(NULL)
})

# Check if the log-rank test was computed correctly
if (!is.null(surv_diff)) {
  p_val <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
  chisq_log_rank <- surv_diff$chisq
} else {
  p_val <- NA
  chisq_log_rank <- NA
}

# Create the Kaplan-Meier plot and add the calculated p-value
if (profile == "cnv") {
  plot_km <- ggsurvplot(
    fit, 
    data = data, 
    conf.int = TRUE,  
    conf.int.alpha = 0.1,  
    conf.int.style = "step",  
    risk.table = TRUE, 
    ggtheme = theme_minimal(), 
    palette = c("#E69F00", "#56B4E9", "#009E73"),  
    linetype = 1,  
    surv.median.line = "hv",
    pval = paste0("p = ", signif(p_val, 3)),  # Manually adds the p-value
    pval.coord = c(90, 0.2)  # Adjusts the position of the p-value in the plot
  )
} else {
  plot_km <- ggsurvplot(
    fit, 
    data = data, 
    conf.int = TRUE,  
    conf.int.alpha = 0.1,  
    conf.int.style = "step",  
    risk.table = TRUE, 
    ggtheme = theme_minimal(), 
    palette = c("#E69F00", "#56B4E9"),  
    linetype = 1,  
    surv.median.line = "hv",
    pval = paste0("p = ", signif(p_val, 3)),  # Manually adds the p-value
    pval.coord = c(90, 0.2)  # Adjusts the position of the p-value in the plot
  )
}


# Define the title based on the metric
title_text <- switch(metric,
                     "OS" = "Overall Survival",
                     "DSS" = "Disease Specific Survival",
                     "DFI" = "Disease Free Interval",
                     "PFI" = "Progression Free Interval",
                     "Overall Survival"  # Default value in case metric does not match any of the options above
)

# Adjust the aesthetics of the main plot
plot_km$plot <- plot_km$plot + 
  theme(panel.grid = element_blank()) +  
  geom_hline(yintercept = 0, color = "black", size = 0.5) +  
  geom_vline(xintercept = 0, color = "black", size = 0.5) +
  labs(title = title_text) +  
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 25, vjust = 0.1))  

# Remove the background grid from the risk table
plot_km$table <- plot_km$table + theme(panel.grid = element_blank())

# Define the output file name
output_file <- paste0("KaplanMeier_", metric, "_", profile, ".pdf")

# Combine the main plot and risk table into a single figure
combined_plot <- plot_grid(plot_km$plot, plot_km$table, ncol = 1, rel_heights = c(3, 1))

# Save the combined plot as a single PDF file
ggsave(filename = output_file, plot = combined_plot, width = 6, height = 6, dpi = 600)

# Message indicating that the file was saved
message("Plot saved as: ", output_file)
