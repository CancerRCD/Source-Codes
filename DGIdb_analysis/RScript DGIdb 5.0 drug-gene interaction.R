##### 
##### 
##### 
##### DGIdb 5.0: drug-gene interaction analysis
##### Based on Gene_Symbols_top_signatures_df1156 top signatures per omic layer (compilation of signatures in Tables, 1, 2, 3 and 4)
#####
#####
# Define the list of required packages
required_packages <- c(
  "dplyr", "stringr", "tidyr", "writexl", "rio", "ggplot2", "boot", 
  "openxlsx", "gganimate", "png", "gifski",
  "geomtextpath", "ggdist", "gghighlight", "ggiraph", "ggpubr", 
  "ggrepel", "ggstatsplot", "ggtext", "patchwork", "ggforce","lemon", "png",
  "gganimate", "plotly", "patchwork", "correctR", "openxlsx", "ggpattern", "forcats"
)

# Function to check, install, and load missing packages
install_and_load_packages <- function(packages) {
  # Find packages that are not installed
  missing_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
  
  # Install any missing packages
  if (length(missing_packages) > 0) {
    message("Installing missing packages: ", paste(missing_packages, collapse = ", "))
    install.packages(missing_packages)
  }
  
  # Load all packages
  sapply(packages, require, character.only = TRUE)
}

# Call the function to install and load packages
install_and_load_packages(required_packages)

# For any GitHub-hosted packages (if CRAN installation fails)
if (!"geomtextpath" %in% installed.packages()[, "Package"]) {
  remotes::install_github("clauswilke/geomtextpath")
}

setwd("D:/Pré-artigo 5-optosis model/Dataframes all metrics/DGIdb 5 drug-gene interaction/DGIdb_05022025")

GDIdb <- import("interactions.tsv")

setwd("D:/Pré-artigo 5-optosis model/Dataframes all metrics/DGIdb 5 drug-gene interaction")

GDIdb_genes <- as.data.frame(unique(GDIdb$gene_name))

# Trim leading and trailing spaces in the "gene_name" column of GDIdb
GDIdb <- GDIdb %>%
  mutate(gene_name = trimws(gene_name))

Gene_Symbols_top_signatures_df1156 <- import("Gene_Symbols_top_signatures_df1156.tsv")

# Trim leading and trailing spaces in the "Gene_Symbol" column of Gene_Symbols_top_signatures_df1156
Gene_Symbols_top_signatures_df1156 <- Gene_Symbols_top_signatures_df1156 %>%
  mutate(Gene_Symbol = trimws(Gene_Symbol))

# Filter GDIdb for rows where "gene_name" matches "Gene_Symbol"
filtered_GDIdb <- GDIdb %>%
  filter(gene_name %in% Gene_Symbols_top_signatures_df1156$Gene_Symbol)

# Filter meaningful interactions (excluding null, other, and unknown interactions)
DatasetS1R <- filtered_GDIdb %>%
  filter(interaction_type != "NULL" & interaction_type != "other/unknown" & drug_name != "NULL") %>%
  distinct()

# Filter meaningful interactions (excluding null, other and unknown interactions)
filtered_GDIdb_meaningful <- filtered_GDIdb %>%
  filter(interaction_type != "NULL" & interaction_type != "other/unknown" & drug_name !="NULL") %>%
  distinct()

rio::export(DatasetS1R, "DatasetS1R.xlsx")

GDIdb_gene_distribution <- filtered_GDIdb_meaningful %>%
  group_by(interaction_type, gene_name) %>%
  summarise(count = n(), .groups = 'drop')

gene_counts <- filtered_GDIdb_meaningful %>%
  group_by(gene_name) %>%
  summarise(total_count = n(), .groups = 'drop')

# Filter genes with at least 5 counts
genes_with_min_counts <- gene_counts[gene_counts$total_count >= 5, ]

# Order by total_count in descending order
genes_with_min_counts <- genes_with_min_counts[order(-genes_with_min_counts$total_count), ]

# Extract the ordered gene names
ordered_gene_list <- genes_with_min_counts$gene_name

print(ordered_gene_list)

# > print(ordered_gene_list)
# [1] "APBB1"   "NAT2"    "ITGB3"   "RHOB"    "TLR4"    "ATP5F1A" "TNFRSF4" "GATA3"   "PARP3"  
# [10] "RPL5" 
# 
set_genes_with_min_counts <- filtered_GDIdb_meaningful %>%
  group_by(gene_name) %>%
  filter(n() >= 5) %>%
  ungroup()

# Create the stacked histogram
ggplot(filtered_GDIdb_meaningful, aes(x = gene_name, fill = interaction_type)) +
  geom_histogram(stat = "count", position = "stack", color = "black") +
  labs(
    title = "Distribution of Gene Names by Interaction Type",
    x = "Gene Name",
    y = "Count",
    fill = "Interaction Type"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    plot.title = element_text(hjust = 0.5)
  )

# Load necessary libraries
library(ggplot2)
library(forcats)

# Define a color palette with sufficient unique colors
interaction_types <- unique(filtered_GDIdb_meaningful$interaction_type)
num_interaction_types <- length(interaction_types)

# Generate a color palette with the required number of unique colors
# You can use the 'hue_pal' function from the 'scales' package to generate distinct colors
library(scales)
color_palette <- hue_pal()(num_interaction_types)

# Create a named vector for the color palette
names(color_palette) <- interaction_types

# Create the stacked bar plot with the specified color palette
ggplot(filtered_GDIdb_meaningful, aes(x = fct_infreq(gene_name), fill = interaction_type)) +
  geom_bar(position = "stack", color = "black") +
  labs(
    title = "Distribution of Gene Names by Interaction Type",
    x = "Gene Name",
    y = "Count",
    fill = "Interaction Type"
  ) +
  scale_fill_manual(values = color_palette) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    plot.title = element_text(hjust = 0.5)
  )

# Install ggpattern if not already installed
if (!requireNamespace("ggpattern", quietly = TRUE)) {
  install.packages("ggpattern", repos = "https://cinc.rud.is")
}

#### Creating bar histogram using different colors and different patterns
# Load required libraries
library(ggplot2)
library(ggpattern)
library(forcats)

# Extract unique interaction types
interaction_types <- unique(filtered_GDIdb_meaningful$interaction_type)
num_interaction_types <- length(interaction_types)

# Define a color palette with sufficient unique colors
color_palette <- scales::hue_pal()(num_interaction_types)

# Define a set of patterns
available_patterns <- c("stripe", "crosshatch", "wave", "weave", "circle", "none")
pattern_palette <- rep(available_patterns, length.out = num_interaction_types)

# Create named vectors for colors and patterns
names(color_palette) <- interaction_types
names(pattern_palette) <- interaction_types

# Create the plot
plot_GDIdb_genes <- ggplot(filtered_GDIdb_meaningful, aes(x = fct_infreq(gene_name), fill = interaction_type, pattern = interaction_type)) +
  geom_bar_pattern(
    position = "stack",
    color = "black",
    pattern_density = 0.3,  # Adjusted for better clarity
    pattern_fill = "grey80",  # Adjusted for better contrast
    pattern_spacing = 0.02
  ) +
  scale_fill_manual(values = color_palette) +
  scale_pattern_manual(values = pattern_palette) +
  labs(
    title = "Distribution of Gene Names by Interaction Type",
    x = "Gene Name",
    y = "Count",
    fill = "Interaction Type",
    pattern = "Interaction Type"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  )

# Save the plot as PDF
ggsave(
  filename = "plot_GDIdb_genes.pdf",
  plot = plot_GDIdb_genes,
  device = "pdf",
  width = 11.69,  # A4 width in inches (landscape orientation)
  height = 8.27,  # A4 height in inches
  dpi = 600
)

# Save the plot as TIFF
ggsave(
  filename = "plot_GDIdb_genes.tiff",
  plot = plot_GDIdb_genes,
  device = "tiff",
  width = 11.69,
  height = 8.27,
  dpi = 600
)

#### 
####
####
#### 
#### Network analysis plotting, Part A. using the meaningful set of interactions
#### Network coloring by drug effect type - Revisited and corrected on 03/03/2025
#### 
####   Load required packages
library(igraph)
library(ggraph)
library(tidyverse)

# Filter meaningful interactions (excluding null, other, and unknown interactions)
filtered_GDIdb_network_meaningful <- filtered_GDIdb %>%
  filter(interaction_type != "NULL" & interaction_type != "other/unknown" & drug_name != "NULL") %>%
  distinct()

# Create an edge list (interactions between genes and drugs)
edges <- filtered_GDIdb_network_meaningful %>%
  select(gene_name, drug_name) %>%
  distinct()

# Create a graph object
g <- graph_from_data_frame(edges, directed = FALSE)

# Assign node types (Gene or Drug)
V(g)$type <- ifelse(V(g)$name %in% filtered_GDIdb_network_meaningful$gene_name, "Gene", "Drug")

# Assign interaction types to drugs (not genes)
drug_interaction_mapping <- filtered_GDIdb_network_meaningful %>%
  distinct(drug_name, interaction_type) %>%
  group_by(drug_name) %>%
  summarise(interaction_type = first(interaction_type), .groups = "drop") # Keep first type if multiple

# Assign interaction types to graph nodes (only for drugs)
V(g)$interaction_type <- ifelse(V(g)$type == "Drug",
                                drug_interaction_mapping$interaction_type[match(V(g)$name, drug_interaction_mapping$drug_name)],
                                NA)

# Generate a color palette for different interaction types (for drugs only)
interaction_types <- unique(na.omit(V(g)$interaction_type))
interaction_colors <- setNames(scales::hue_pal()(length(interaction_types)), interaction_types)

# Combine colors for genes and drugs
node_colors <- ifelse(V(g)$type == "Gene", "steelblue", interaction_colors[V(g)$interaction_type])

# Create a named color vector for the legend
legend_colors <- c("Gene" = "steelblue", interaction_colors)

# Plot using ggraph
Network_GDIdb_genes <- ggraph(g, layout = "fr") + 
  geom_edge_link(aes(edge_alpha = 0.5), color = "gray50") + 
  geom_node_point(aes(color = ifelse(V(g)$type == "Gene", "Gene", V(g)$interaction_type)), size = 5) + 
  geom_node_text(aes(label = name), repel = TRUE, size = 4) + 
  scale_color_manual(name = "Node Type", values = legend_colors) +  # Use scale_color_manual for legend
  theme_void() + 
  ggtitle("Drug-Gene Interaction Network: Colored by Interaction Type") + 
  theme(legend.position = "right")

# Save the Network plot as PDF
ggsave(
  filename = "Network_plot_GDIdb_genes.pdf",
  plot = Network_GDIdb_genes,
  device = "pdf",
  width = 14,  # A4 width in inches (landscape orientation)
  height = 10,  # A4 height in inches
  dpi = 600
)

# Save the Network plot as TIFF
ggsave(
  filename = "Network_plot_GDIdb_genes.tiff",
  plot = Network_GDIdb_genes,
  device = "tiff",
  width = 14,  # A4 width in inches (landscape orientation)
  height = 10,  # A4 height in inches
  dpi = 600
)

# Check for mismatches between graph nodes and gene/drug names in the dataframe

# Extract unique gene and drug names from the dataframe
genes_in_data <- unique(filtered_GDIdb_network_meaningful$gene_name)
drugs_in_data <- unique(filtered_GDIdb_network_meaningful$drug_name)

# Combine gene and drug names from the dataframe
all_names_in_data <- unique(c(genes_in_data, drugs_in_data))

# Extract node names from the graph
node_names_in_graph <- V(g)$name

# Find mismatches (nodes in the graph not found in the dataframe)
mismatches <- setdiff(node_names_in_graph, all_names_in_data)

# Print mismatches
if (length(mismatches) > 0) {
  cat("Mismatches found (nodes in the graph not found in the dataframe):\n")
  print(mismatches)
} else {
  cat("No mismatches found. All nodes in the graph are present in the dataframe.\n")
}
####
#### Interactive Network
####
# Load required packages
library(igraph)
library(visNetwork)
library(tidyverse)
library(htmlwidgets)  # Required to save HTML output

# Filter meaningful interactions (excluding null, other, and unknown interactions)
filtered_GDIdb_network_meaningful <- filtered_GDIdb %>%
  filter(interaction_type != "NULL" & interaction_type != "other/unknown" & drug_name != "NULL") %>%
  distinct()

# Create an edge list (interactions between genes and drugs)
edges <- filtered_GDIdb_network_meaningful %>%
  select(gene_name, drug_name) %>%
  distinct()

# Create a graph object
g <- graph_from_data_frame(edges, directed = FALSE)

# Check for mismatches between graph nodes and gene/drug names in the dataframe
genes_in_data <- unique(filtered_GDIdb_network_meaningful$gene_name)
drugs_in_data <- unique(filtered_GDIdb_network_meaningful$drug_name)
all_names_in_data <- unique(c(genes_in_data, drugs_in_data))
node_names_in_graph <- V(g)$name
mismatches <- setdiff(node_names_in_graph, all_names_in_data)

if (length(mismatches) > 0) {
  cat("Mismatches found (nodes in the graph not found in the dataframe):\n")
  print(mismatches)
} else {
  cat("No mismatches found. All nodes in the graph are present in the dataframe.\n")
}

# Assign node types (Gene or Drug)
V(g)$type <- ifelse(V(g)$name %in% filtered_GDIdb_network_meaningful$gene_name, "Gene", "Drug")

# Assign interaction types to drugs (not genes)
drug_interaction_mapping <- filtered_GDIdb_network_meaningful %>%
  distinct(drug_name, interaction_type) %>%
  group_by(drug_name) %>%
  summarise(interaction_type = first(interaction_type), .groups = "drop")

# Assign interaction types to graph nodes (only for drugs)
V(g)$interaction_type <- ifelse(V(g)$type == "Drug",
                                drug_interaction_mapping$interaction_type[match(V(g)$name, drug_interaction_mapping$drug_name)],
                                NA)

# Generate a distinct color palette for different interaction types
interaction_types <- unique(na.omit(V(g)$interaction_type))
interaction_colors <- setNames(scales::hue_pal()(length(interaction_types)), interaction_types)

# Assign colors: Genes in "steelblue", Drugs colored by interaction type
V(g)$color <- ifelse(V(g)$type == "Gene", "steelblue", interaction_colors[V(g)$interaction_type])

# Convert graph into a data frame format for visNetwork
nodes <- data.frame(
  id = V(g)$name,
  label = V(g)$name,
  color = V(g)$color,
  group = ifelse(V(g)$type == "Gene", "Gene", V(g)$interaction_type), # Use drug effect type as group
  title = ifelse(V(g)$type == "Gene",
                 paste0("<b>Gene:</b> ", V(g)$name),
                 paste0("<b>Drug:</b> ", V(g)$name, "<br><b>Effect:</b> ", V(g)$interaction_type))
)

edges <- data.frame(
  from = edges$gene_name,
  to = edges$drug_name
)

# Initialize the network
network <- visNetwork(nodes, edges) %>%
  visEdges(smooth = FALSE) %>%
  visNodes(size = 10, font = list(size = 20)) %>%
  visInteraction(navigationButtons = TRUE, tooltipDelay = 50)

# Add the "Gene" group
network <- network %>% visGroups(groupname = "Gene", color = "steelblue")

# Add each drug effect type separately in a loop
for (interaction in interaction_types) {
  network <- network %>% visGroups(groupname = interaction, color = interaction_colors[[interaction]])
}

# Add the correctly formatted legend
network <- network %>%
  visLegend(useGroups = TRUE, main = "Drug Effect Categorization") %>%
  visPhysics(enabled = TRUE) %>% 
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) 

# Save as HTML
saveWidget(network, "Drug_Gene_Interaction_Network.html", selfcontained = TRUE)

# Display the network in RStudio Viewer
network

####
####
####
#####
##### 
##### 
######################
######################
# Save the entire workspace to a file
save.image(file = "workspace.RData")
#####################
#####################
#####
#####
#####
#####


#####
#####
#####
# Load the saved workspace into the current session
load("workspace.RData")
#####
#####
#####
