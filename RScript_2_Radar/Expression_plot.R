# Load necessary libraries
library(fmsb)
library(rio)
library(ggpubr)
library(dplyr)
library(rstatix)
library(ggplot2)
library(cowplot)

###########  Save PDF in inches 5x5  #################

# Set working directory
setwd("C:/Users/Emamnuell/Desktop/panoptosis/genes/genes_24_optosis/codigos/plots/BRCA-1496.1.3.P.3.71.71.1.1.2")

# Import data
df <- import("expression.xlsx")

# Perform Wilcoxon test
stat.test <- wilcox.test(value ~ type2, data = df) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

# Create boxplot using ggplot2
p <- ggboxplot(df, x = "type2", y = "value", fill = "type2", legend = "none") +
  labs(x = "", y = "mRNA Expression (log2(tpm+0.001))") +
  # scale_fill_manual(values = c("#56B4E9", "#E69F00")) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
  scale_x_discrete(labels = c("Tumor", "Normal")) +
  ylim(-10, 10) +
  ggtitle("Expression") +  # Plot title
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20), 
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    plot.margin = margin(10, 10, 10, 10)  # Adjust margins to provide more space
  )

# Add p-value
my_comparisons <- list(c("tumor", "normal"))
p <- p + stat_compare_means(method = "wilcox.test", comparisons = my_comparisons,
                            symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), 
                                               symbols = c("****", "***", "**", "*", "ns")),
                            size = 6,  
                            tip.length = 0.02,  
                            bracket.size = 0.4,  
                            step.increase = 0.1)

# Display the plot
print(p)

# Save the plot as a PDF
ggsave("expression_plot.pdf",  height = 5, width = 5, units = "in")
