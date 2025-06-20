library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)

setwd('/Users/zkarwowska/new_embl_folder/zeevi_dataset_v5/real_datasets/')

# Load data
counts <- read_csv("palleja_dataset/ready_palleja_counts.csv") %>%
  column_to_rownames(var = colnames(.)[1])

metadata <- read_csv("palleja_dataset/ready_palleja_metadata.csv") %>%
  column_to_rownames(var = colnames(.)[1]) %>%
  filter(Timepoint %in% c(0, 4))

# Create named vector for mapping sample ID to timepoint
timepoint_dict <- setNames(metadata$Timepoint, metadata$Sample_ID)

# Reindex rows of counts using timepoint_dict
counts <- counts[ , metadata$Sample_ID] %>% t() %>% as.data.frame()
counts$Timepoint <- timepoint_dict[rownames(counts)]

# Group by Timepoint and calculate mean (+ pseudocount)
counts_mean_df <- counts %>%
  group_by(Timepoint) %>%
  summarise(across(everything(), mean)) %>%
  column_to_rownames("Timepoint") %>%
  as.matrix() + 1e-3

# Calculate absolute log fold change between timepoints 0 and 4
logFC <- log2(counts_mean_df["0", ] / counts_mean_df["4", ])

# Plot using ggplot2
logFC_df <- data.frame(logFC = logFC)

ggplot(logFC_df, aes(y = logFC, x = 'feature')) + 
  geom_boxplot(
    width = 0.3,
    alpha = 0.7,
    fill = "#4CAF50",
    color = "#4CAF50",
    outlier.color = "#4CAF50",
    outlier.size = 2
  ) +
  geom_jitter(
    aes(y = logFC),
    width = 0.25,
    alpha = 0.6,
    color = "#2E7D32",
    size = 1.5
  ) +
  labs(
    title = "Distribution of Log Fold Changes",
    subtitle = "Boxplot with individual data points",
    y = "Log Fold Change (logFC)",
    x = NULL
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray60"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )
