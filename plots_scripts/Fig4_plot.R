library(dplyr)
library(ggplot2)
library(gridExtra)
library(glue)

# Install and load cowplot if needed
if (!requireNamespace("cowplot", quietly = TRUE)) {
  install.packages("cowplot")
}
library(cowplot)

# Define colors manually
color_palette <- c(
  "linda" = "#E41A1C",  # Red
  "maaslin.clr" = "#0358f1",  # Lighter Blue
  "maaslin.log_relab" = "#70d6ff",  # Darker Blue
  "metasplines" = "#4DAF4A",  # Green
  "spliner.clr" = "#ff6d00",  # Lighter Yellow
  "spliner.log_relab" = "#ffd400",  # Darker Yellow
  "wilcoxon.clr" = "#9d4edd",  # Lighter Brown
  "wilcoxon.log_relab" = "#fb6f92",  # Darker Brown
  "nbmm" = "#A65628",  # Purple
  "nbzimm" = "#999999",  # Orange
  "zigmm" = "#F781BF"   # Pink
)

# Define effect sizes
efs_values <- c(125, 15, 2, 3, 5)
efs_labels <- c("EFS = 1.25", "EFS = 1.5", "EFS = 2", "EFS = 3", "EFS = 5")

# Define sampling frequencies we want to compare
sampling_frequencies <- c(2, 4)

# Read all CSV files
results_dir <- '/Users/zkarwowska/new_embl_folder/zeevi_dataset_v5/results/results_files/one_arm/'
pattern <- ".*.csv"
files <- list.files(path = results_dir, pattern = pattern, full.names = TRUE)
all_data <- do.call(rbind, lapply(files, read.csv))
all_data[is.na(all_data)] <- 1

# Create a unique model identifier
all_data$transformation[all_data$transformation == 'none'] <- 'None'
all_data$transformation[all_data$transformation == 'log'] <- 'log_relab'
all_data$transformation[all_data$transformation == 'relab'] <- 'log_relab'
all_data$model[all_data$model == 'trendyspliner'] <- 'spliner'
all_data$model_name <- ifelse(is.na(all_data$transformation) | all_data$transformation == "None",
                              all_data$model,
                              paste0(all_data$model, ".", all_data$transformation))

# Function to process data for a given effect size and sampling frequency
process_effect_size <- function(efs_value, efs_label, sampling_freq) {
  df_recall <- all_data %>%
    filter(effect_size == efs_value, sampling == sampling_freq) %>%
    group_by(sample_size, model_name) %>%
    summarise(Mean = mean(recall, na.rm = TRUE), .groups = "drop") %>%
    mutate(Metric = "Recall", 
           EFS = efs_label, 
           Sampling = paste("F =", sampling_freq),
           EFS_Sampling = paste(efs_label, "F =", sampling_freq))
  
  df_precision <- all_data %>%
    filter(effect_size == efs_value, sampling == sampling_freq) %>%
    group_by(sample_size, model_name) %>%
    summarise(Mean = mean(precision, na.rm = TRUE), .groups = "drop") %>%
    mutate(Metric = "Precision", 
           EFS = efs_label, 
           Sampling = paste("F =", sampling_freq),
           EFS_Sampling = paste(efs_label, "F =", sampling_freq))
  
  df_auc <- all_data %>%
    filter(effect_size == efs_value, sampling == sampling_freq) %>%
    group_by(sample_size, model_name) %>%
    summarise(Mean = mean(auc, na.rm = TRUE), .groups = "drop") %>%
    mutate(Metric = "AUC", 
           EFS = efs_label, 
           Sampling = paste("F =", sampling_freq),
           EFS_Sampling = paste(efs_label, "F =", sampling_freq))
  
  return(bind_rows(df_recall, df_precision, df_auc))
}

# Apply function to all effect sizes and sampling frequencies, then combine results
df_combined <- bind_rows(
  mapply(process_effect_size, 
         rep(efs_values, each = length(sampling_frequencies)),
         rep(efs_labels, each = length(sampling_frequencies)),
         rep(sampling_frequencies, times = length(efs_values)),
         SIMPLIFY = FALSE)
)

# Reorder the metrics to Recall, Precision, then AUC  
df_combined$Metric <- factor(df_combined$Metric, levels = c("Recall", "Precision", "AUC"))
df_combined$model_name <- gsub(".none", "", df_combined$model_name)

# Apply the same filters as in the original code
df_combined <- df_combined %>% 
  filter(!model_name %in% c('baseline', 'baseline.clr', 'baseline.log', 'baseline.log_relab')) %>%
  filter((model_name == 'spliner.clr' & EFS == 'EFS = 3') | model_name != 'spliner.clr') %>%
  filter((model_name == 'spliner.log_relab' & EFS == 'EFS = 3') | model_name != 'spliner.log_relab')

df_combined <- df_combined %>%
  filter(!(EFS %in% c('EFS = 1.25', 'EFS = 1.5', 'EFS = 2', 'EFS = 5') & 
             model_name %in% c('spliner.clr', 'spliner.relab', 'metasplines')))

# Define the y-limits dynamically in the data
df_combined <- df_combined %>%
  mutate(y_min = ifelse(Metric == "AUC", 0.5, 0),
         y_max = ifelse(Metric == "AUC", 1.1, 1.1))

library(dplyr)
library(ggplot2)
library(gridExtra)
library(glue)

# Define colors manually
color_palette <- c(
  "linda" = "#E41A1C",  # Red
  "maaslin.clr" = "#0358f1",  # Lighter Blue
  "maaslin.log_relab" = "#70d6ff",  # Darker Blue
  "metasplines" = "#4DAF4A",  # Green
  "spliner.clr" = "#ff6d00",  # Lighter Yellow
  "spliner.log_relab" = "#ffd400",  # Darker Yellow
  "wilcoxon.clr" = "#9d4edd",  # Lighter Brown
  "wilcoxon.log_relab" = "#fb6f92",  # Darker Brown
  "nbmm" = "#A65628",  # Purple
  "nbzimm" = "#999999",  # Orange
  "zigmm" = "#F781BF"   # Pink
)

# Define effect sizes
efs_values <- c(125, 15, 2, 3, 5)
efs_labels <- c("EFS = 1.25", "EFS = 1.5", "EFS = 2", "EFS = 3", "EFS = 5")

# Define sampling frequencies we want to compare
sampling_frequencies <- c(2, 4)

# Read all CSV files
results_dir <- '/Users/zkarwowska/new_embl_folder/zeevi_dataset_v5/results/results_files/one_arm/'
pattern <- ".*.csv"
files <- list.files(path = results_dir, pattern = pattern, full.names = TRUE)
all_data <- do.call(rbind, lapply(files, read.csv))
all_data[is.na(all_data)] <- 1

# Create a unique model identifier
all_data$transformation[all_data$transformation == 'none'] <- 'None'
all_data$transformation[all_data$transformation == 'log'] <- 'log_relab'
all_data$transformation[all_data$transformation == 'relab'] <- 'log_relab'
all_data$model[all_data$model == 'trendyspliner'] <- 'spliner'
all_data$model_name <- ifelse(is.na(all_data$transformation) | all_data$transformation == "None",
                              all_data$model,
                              paste0(all_data$model, ".", all_data$transformation))

# Function to process data for a given effect size and sampling frequency
process_effect_size <- function(efs_value, efs_label, sampling_freq) {
  df_recall <- all_data %>%
    filter(effect_size == efs_value, sampling == sampling_freq) %>%
    group_by(sample_size, model_name) %>%
    summarise(Mean = mean(recall, na.rm = TRUE), .groups = "drop") %>%
    mutate(Metric = "Recall", 
           EFS = efs_label, 
           Sampling = paste("F =", sampling_freq))
  
  df_precision <- all_data %>%
    filter(effect_size == efs_value, sampling == sampling_freq) %>%
    group_by(sample_size, model_name) %>%
    summarise(Mean = mean(precision, na.rm = TRUE), .groups = "drop") %>%
    mutate(Metric = "Precision", 
           EFS = efs_label, 
           Sampling = paste("F =", sampling_freq))
  
  df_auc <- all_data %>%
    filter(effect_size == efs_value, sampling == sampling_freq) %>%
    group_by(sample_size, model_name) %>%
    summarise(Mean = mean(auc, na.rm = TRUE), .groups = "drop") %>%
    mutate(Metric = "AUC", 
           EFS = efs_label, 
           Sampling = paste("F =", sampling_freq))
  
  return(bind_rows(df_recall, df_precision, df_auc))
}

# Apply function to all effect sizes and sampling frequencies, then combine results
df_combined <- bind_rows(
  mapply(process_effect_size, 
         rep(efs_values, each = length(sampling_frequencies)),
         rep(efs_labels, each = length(sampling_frequencies)),
         rep(sampling_frequencies, times = length(efs_values)),
         SIMPLIFY = FALSE)
)

# Reorder the metrics to Recall, Precision, then AUC  
df_combined$Metric <- factor(df_combined$Metric, levels = c("Recall", "Precision", "AUC"))
df_combined$model_name <- gsub(".none", "", df_combined$model_name)

# Apply the same filters as in the original code
df_combined <- df_combined %>% 
  filter(!model_name %in% c('baseline', 'baseline.clr', 'baseline.log', 'baseline.log_relab')) %>%
  filter((model_name == 'spliner.clr' & EFS == 'EFS = 3') | model_name != 'spliner.clr') %>%
  filter((model_name == 'spliner.log_relab' & EFS == 'EFS = 3') | model_name != 'spliner.log_relab')

df_combined <- df_combined %>%
  filter(!(EFS %in% c('EFS = 1.25', 'EFS = 1.5', 'EFS = 2', 'EFS = 5') & 
             model_name %in% c('spliner.clr', 'spliner.relab', 'metasplines')))

# Define the y-limits dynamically in the data
df_combined <- df_combined %>%
  mutate(y_min = ifelse(Metric == "AUC", 0.5, 0),
         y_max = ifelse(Metric == "AUC", 1.1, 1.1))

# Create factors with the correct order
df_combined <- df_combined %>%
  mutate(
    EFS = factor(EFS, levels = efs_labels),
    Sampling = factor(Sampling, levels = paste("F =", sampling_frequencies))
  )

# Create a combined plot with boxed regions for each effect size
# Creating individual plots for each effect size and then combining them

plot_list <- list()
for (efs in efs_labels) {
  df_efs <- df_combined %>% filter(EFS == efs)
  
  # Create plot for this effect size
  p <- ggplot(df_efs, aes(x = sample_size, y = Mean, color = model_name, group = model_name)) +
    geom_line(size = 1.2) +
    facet_grid(Metric ~ Sampling, scales = "free_y") +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 10, face = "bold"),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      plot.background = element_rect(fill = "white", color = "black", size = 1.5),
      legend.position = "none",  # Remove legend from individual plots
      plot.margin = margin(15, 15, 15, 15),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    ) +
    labs(
      title = efs,
      x = NULL,  # Remove x label from individual plots
      y = NULL   # Remove y label from individual plots
    ) +
    scale_color_manual(values = color_palette)
  
  plot_list[[efs]] <- p
}

# Use grid.arrange to create the combined plot with a common legend
combined_plot <- do.call(grid.arrange, c(plot_list, ncol = 5))

# Create a legend separately
legend_plot <- ggplot(df_combined, aes(x = sample_size, y = Mean, color = model_name)) +
  geom_line() +
  scale_color_manual(values = color_palette) +
  theme(legend.position = "bottom")

# Extract the legend
legend <- cowplot::get_legend(legend_plot)

# Add the legend to the combined plot with grid.arrange
final_plot <- grid.arrange(
  combined_plot, 
  legend,
  heights = c(10, 1),
  ncol = 1
)

# Save the final plot
ggsave(
  final_plot,
  file = "/Users/zkarwowska/Desktop/EMBL_project/zeevi_dataset_v5/all_effect_sizes_sampling_comparison.pdf",
  width = 20,
  height = 5,
  dpi = 300,
  limitsize = FALSE
)

print(final_plot)
