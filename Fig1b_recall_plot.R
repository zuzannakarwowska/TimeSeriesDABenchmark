library(dplyr)
library(ggplot2)
library(gridExtra)
library(purrr)
library(cowplot)

# Load required package
library(RColorBrewer)

# Define colors manually
color_palette <- c(
  "linda" = "#E41A1C",  # Red
  "maaslin.clr" = "#0358f1",  # Lighter Blue
  "maaslin.log" = "#70d6ff",  # Darker Blue1F5B99
  "metasplines" = "#4DAF4A",  # Green
  "trendyspliner.clr" = "#ff6d00",  # Lighter Yellow
  "trendyspliner.log" = "#ffd400",  # Darker Yellow
  "wilcoxon.clr" = "#9d4edd",  # Lighter Brown
  "wilcoxon.log_relab" = "#fb6f92",  # Darker Brown
  "nbmm" = "#A65628",  # Purple
  "nbzimm" = "#999999",  # Orange
  "zigmm" = "#F781BF"   # Pink
)



# Read all CSV files
results_dir <- '/Users/zkarwowska/new_embl_folder/zeevi_dataset_v5/results/results_files/one_arm/'
pattern <- ".*.csv"
files <- list.files(path = results_dir, pattern = pattern, full.names = TRUE)
all_data <- do.call(rbind, lapply(files, read.csv))
all_data[is.na(all_data)] <- 1

# Create a unique model identifier
all_data$model_name <- ifelse(is.na(all_data$transformation) | all_data$transformation == "None", 
                              all_data$model, 
                              paste0(all_data$model, ".", all_data$transformation))

efs_values <- c(125, 15, 2, 3, 5)
efs_labels <- c("EFS = 1.25", "EFS = 1.5", "EFS = 2", "EFS = 3", "EFS = 5")

sampling_values <- c(2, 4, 6)
sampling_labels <- c("M = 2", "M = 4", "M = 6")

process_effect_size <- function(efs_value, efs_label, sampling_value, sampling_label) {
  
  df1 <- all_data %>%
    filter(effect_size == efs_value, sampling == sampling_value) %>%
    group_by(sample_size, model_name) %>%
    summarise(Mean = mean(recall, na.rm = TRUE), .groups = "drop") %>%
    mutate(Metric = "Recall", EFS = efs_label, M=sampling_label)

  return(df1)
  
}

result_list <- list()
idx <- 1
for (i in seq_along(efs_values)) {
  for (j in seq_along(sampling_values)) {
    result_list[[idx]] <- process_effect_size(
      efs_values[i], 
      efs_labels[i],
      sampling_values[j],
      sampling_labels[j]
    )
    idx <- idx + 1
  }
}

df_combined <- bind_rows(result_list)
df_combined$M <- factor(df_combined$M, levels = c("M = 2", "M = 4", "M = 6"))

# Create the plot
plot <- ggplot(df_combined, aes(x = sample_size, y = Mean, color = model_name, group = model_name)) +
  geom_line(size = 1.2) +
  facet_grid(M ~ EFS, scales = "free_y") +  # Allow different y scales
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 1) # Add black borders
  ) +
  ylim(-.05, 1) + 
  scale_color_manual(values = color_palette) + 
  labs(x = "Sample Size", y = "Mean Recall", color = "Model") +
  theme(strip.text = element_text(size = 12, face = "bold"))

print(plot)
