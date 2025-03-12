library(ggplot2)
library(tidyverse)
library(rio)

# Set working directory
setwd("C:/Users/Emamnuell/Desktop/UENF/laboratorio/panoptosis/genes/genes_24_optosis/codigos/plots/ZZ_til")

# Import data
data <- import("til_KIRP-107.3.2.N.1.44.44.1.1.2.xlsx")

# Remove the suffix after "_"
data$immune_cells <- sub("_.*", "", data$immune_cells)

# Check if 'cor' is numeric and remove NAs and infinite values
data$cor <- as.numeric(data$cor)
data <- data %>%
  filter(!is.na(immune_cells) & !is.na(cor) & is.finite(cor))

# Manually define the min and max values for the color scale
min_cor <- -1
max_cor <- 1

# Create the plot in horizontal format
plot <- ggplot(data, aes(x = immune_cells, y = Cancer_type, fill = cor)) +
  geom_tile(color = "white", size = 0.5) + # Reduce the border size of the cells
  scale_fill_gradient2(low = "#56B4E9", high = "#E69F00", mid = "white", midpoint = 0, 
                       limit = c(min_cor, max_cor), space = "Lab", 
                       name = "Correlation") +
  theme_minimal() +
  labs(x = "Immune Cells", y = "Cancer Type", title = "Immune Infiltrates") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), # Adjust X-axis label rotation
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(vjust = -0.5, size = 17), 
        axis.title.y = element_text(vjust = 5, size = 17),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 27, vjust = 4),
        legend.position = "right") +
  scale_y_discrete(limits = c("KIRP")) + 
  scale_x_discrete(expand = c(0.1, 0.1)) + # Give more space to X-axis labels
  coord_fixed(ratio = 1) + # Maintain the appropriate aspect ratio
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) # Adjust margins

# Display the plot
print(plot)

# Save the plot in horizontal format
ggsave("imune_horizontal.pdf", plot = plot, units = "cm", height = 16, width = 34, dpi = 600)
