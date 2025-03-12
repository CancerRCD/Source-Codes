# Import necessary libraries
library(fmsb)
library(rio)
library(ggpubr)
library(dplyr)
library(rstatix)
library(ggplot2)
library(cowplot)

# Set working directory
setwd("C:/Users/Emamnuell/Desktop/panoptosis/genes/genes_24_optosis/codigos/plots/KIRC-169.2.1.P.2.71.45.1.1.2")

# Import data
data <- import("radar.csv")

# Filter for cancers of interest: BRCA, KIRP, and BLCA
cancers_of_interest <- c("THYM", "SARC", "BRCA", "CHOL", "ESCA", "KIRC")

data_filtered <- data %>% filter(cancer %in% cancers_of_interest)

# Add rows for maximum and minimum values
max_min <- data.frame(matrix(c(rep(1, 6), rep(-1, 6)), 
                             ncol = 6, byrow = TRUE))
colnames(max_min) <- data_filtered$cancer

# Transpose the 'cor' column and add it to the max_min dataframe
values <- t(data_filtered$cor)
colnames(values) <- data_filtered$cancer

# Combine the data
data_for_radar <- rbind(max_min, values)

# Create the radar plot and save it as an object
radar_plot <- function() {
  radarchart(data_for_radar, axistype = 1,
             pcol = "#56B4E9",  
             pfcol = scales::alpha("#56B4E9", 0.1),
             plwd = 2,  
             cglcol = "grey",  
             cglty = 1,  
             axislabcol = "black",  
             caxislabels = seq(-1, 1, 0.5),  
             vlcex = 2,  # Increase font size for labels at vertices
             cex.main = 2  # Increase title size
  )
  title(main = "Tumor Mutation Burden", cex.main = 3)
}

# Save the radar plot as a PDF
radar_plot()  # Generate the plot

ggsave("radar_plot.pdf",  height = 9, width = 9, units = "in")
#################      save with inches at 9x9     ########################
