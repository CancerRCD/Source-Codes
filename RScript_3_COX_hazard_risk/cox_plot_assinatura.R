# Set working directory
setwd("C:/Users/Emamnuell/Desktop/UENF/laboratorio/panoptosis/genes/genes_24_optosis/codigos/plots/KIRP-107.3.2.N.1.44.44.1.1.2")

# Load necessary libraries
library(ggplot2)
library(tidyverse)
library(rio)

# Load different datasets
dados_OS <- import("cox_OS.csv")
dados_DFI <- import("cox_DFI.csv")
dados_DSS <- import("cox_DSS.csv")
dados_PFI <- import("cox_PFI.csv")

# Add a column to identify the source of the data
dados_OS$tipo <- "Overall Survival"
dados_DFI$tipo <- "Disease Free Interval"
dados_DSS$tipo <- "Disease Specific Survival"
dados_PFI$tipo <- "Progression Free Interval"

# Combine all datasets into a single dataframe
dados <- bind_rows(dados_PFI, dados_DSS, dados_DFI, dados_OS)

# Filter the data to include only cancers of interest
dados <- dados %>% filter(cancer %in% c("READ", "KIRP", "ESCA", "MESO"))

# Define the order of factors for the "cancer" column in the desired order
dados$cancer <- factor(dados$cancer, levels = c("READ", "KIRP", "ESCA", "MESO"))

# Convert the "cancer" variable to factor
dados$cancer <- factor(dados$cancer)

# Convert HR_log, lower_95_log, and upper_95_log variables to numeric
dados$HR_log <- as.numeric(dados$HR_log)
dados$lower_95_log <- as.numeric(dados$lower_95_log)
dados$upper_95_log <- as.numeric(dados$upper_95_log)

# Create a new column to define point colors
dados$color <- ifelse(dados$lower_95_log > 0 | dados$upper_95_log < 0, 
                      ifelse(dados$HR_log > 0, "risky", "protective"), 
                      "ns")

# Function to create the plot with a title
criar_grafico <- function(dados, tipo_metric) {
  dados_filtrados <- dados %>% filter(tipo == tipo_metric)
  
  plot <- ggplot(dados_filtrados, aes(x = factor(cancer, levels = rev(levels(cancer))), y = HR_log)) +
    geom_point(aes(color = color), size = 2) +  # Add colored points
    geom_errorbar(aes(ymin = lower_95_log, ymax = upper_95_log, color = color), width = 0.2) +  # Add colored error bars
    geom_hline(yintercept = 0, linetype = "solid", size = 0.5, color = "#009E73") +  # Add central line
    theme_classic() +
    labs(x = "Cancer", y = "HR_log") +  # Add titles and labels
    ggtitle(tipo_metric) +  # Use the metric name as the plot title
    theme(
      axis.text.x = element_text(size = 20, hjust = 0.5),  # Adjust X-axis text and size
      axis.text.y = element_text(size = 20),  # Adjust Y-axis text size (cancer names)
      axis.title.x = element_text(size = 20, vjust = 1),  # Adjust X-axis title and size
      axis.title.y = element_text(size = 20),  # Adjust Y-axis title and size
      strip.text = element_text(size = 30, hjust = 0.5, face = "bold"),  # Adjust facet text size
      strip.background = element_blank(),  # Remove the rectangle around facets
      legend.text = element_text(size = 15),  # Adjust legend text size
      legend.title = element_text(size = 20),  # Adjust legend title size
      legend.position = "bottom",  # Position the legend below the plot
      plot.title = element_text(size = 32, face = "bold", hjust = 0.5, vjust = 0.5)  # Increase font size, bold, and center the title
    ) + 
    coord_flip() +  # Flip axes
    scale_color_manual(values = c("risky" = "#E69F00", "protective" = "#56B4E9", "ns" = "black"), 
                       name = "Type", 
                       labels = c("risky" = "Risky", "protective" = "Protective", "ns" = "ns"))  # Define colors and labels
  
  return(plot)
}

# Create and save plots for each metric with specific titles
plot_OS <- criar_grafico(dados, "Overall Survival")
plot_OS
ggsave("cox_plot_OS.pdf", plot = plot_OS, height = 8, width = 8)

plot_DFI <- criar_grafico(dados, "Disease Free Interval")
plot_DFI
ggsave("cox_plot_DFI.pdf", plot = plot_DFI, height = 8, width = 8)

plot_DSS <- criar_grafico(dados, "Disease Specific Survival")
plot_DSS
ggsave("cox_plot_DSS.pdf", plot = plot_DSS, height = 8, width = 8)

plot_PFI <- criar_grafico(dados, "Progression Free Interval")
plot_PFI
ggsave("cox_plot_PFI.pdf", plot = plot_PFI, height = 8, width = 8)
