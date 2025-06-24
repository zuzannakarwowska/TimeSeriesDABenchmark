library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)

setwd('/Users/zkarwowska/new_embl_folder/zeevi_dataset_v5/real_datasets/')

# Load data
counts <- read_csv('dietary_interventions_otu_table.csv') %>%
  column_to_rownames(var = colnames(.)[1])

metadata <- read_csv("dietary_interventions_metadata.csv") %>%
  column_to_rownames(var = colnames(.)[1]) %>% filter(study != 'baxter')

counts <- counts %>% select(metadata$Sample_ID)

# Get unique studies
studies <- unique(metadata$study)

# Initialize list to store results
logFC_results <- list()

# Calculate logFC for each study
for(study_name in studies) {
  cat("Processing study:", study_name, "\n")
  
  # Filter metadata for current study
  dataset_metadata <- metadata %>% 
    filter(study == study_name) %>% 
    dplyr::select(-c("study", 'treatment'))
  
  # Filter counts for current study samples
  dataset_counts <- counts %>% 
    dplyr::select(all_of(dataset_metadata$Sample_ID))
  
  # Create named vector for mapping sample ID to timepoint
  timepoint_dict <- setNames(dataset_metadata$Timepoint, dataset_metadata$Sample_ID)
  
  # Reindex rows of counts using timepoint_dict
  study_counts <- dataset_counts[, dataset_metadata$Sample_ID] %>% 
    t() %>% 
    as.data.frame()
  study_counts$Timepoint <- timepoint_dict[rownames(study_counts)]
  
  # Check if both timepoints 0 and 1 exist
  available_timepoints <- unique(study_counts$Timepoint)
  if(!all(c("0", "1") %in% available_timepoints)) {
    cat("Warning: Study", study_name, "does not have both timepoints 0 and 1. Available:", 
        paste(available_timepoints, collapse = ", "), "\n")
    next
  }
  
  # Group by Timepoint and calculate mean (+ pseudocount)
  counts_mean_df <- study_counts %>%
    group_by(Timepoint) %>%
    summarise(across(everything(), mean, .names = "{.col}")) %>%
    column_to_rownames("Timepoint") %>%
    as.matrix() + 1e-3
  
  # Calculate log fold change between timepoints 1 and 0
  logFC <- log2(counts_mean_df["1", ] / counts_mean_df["0", ])
  
  # Store results
  logFC_results[[study_name]] <- data.frame(
    study = study_name,
    feature = names(logFC),
    logFC = as.numeric(logFC),
    stringsAsFactors = FALSE
  )
}

# Combine all results
all_logFC <- do.call(rbind, logFC_results)

# Remove any infinite or NaN values
all_logFC <- all_logFC %>%
  filter(is.finite(logFC))

# Create combined data with individual studies + "All Studies" category
all_logFC_combined <- rbind(
  all_logFC,
  data.frame(
    study = "All Studies",
    feature = all_logFC$feature,
    logFC = all_logFC$logFC,
    stringsAsFactors = FALSE
  )
)

# Reorder factor levels to put "All Studies" at the end
all_logFC_combined$study <- factor(all_logFC_combined$study, 
                                   levels = c(unique(all_logFC$study), "All Studies"))

# Create boxplot with individual studies + combined "All Studies"
g <- ggplot(all_logFC_combined, aes(x = study, y = logFC, fill = study)) + 
  geom_boxplot(
    width = 0.6,
    alpha = 0.7,
    outlier.size = 1,
    outlier.alpha = 0.6
  ) +
  geom_jitter(
    width = 0.2,
    alpha = 0.4,
    size = 0.8
  ) +
  labs(
    title = "Distribution of Log Fold Changes by Study",
    subtitle = "Individual studies and combined data across all studies",
    y = "Log Fold Change (logFC)",
    x = "Study"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5, color = "gray60"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 12),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  scale_fill_viridis_d(option = "viridis")

print(g)

# Save the plot
ggsave(
  filename = "/Users/zkarwowska/Desktop/dietary_logFC_by_study.pdf",
  plot = g,
  width = 10,
  height = 7,
  units = "in",
  dpi = 300,
  device = "pdf"
)
