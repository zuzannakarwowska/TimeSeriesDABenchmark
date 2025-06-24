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

#df_combined <- df_combined %>%
#  filter(!(EFS %in% c('EFS = 1.25', 'EFS = 1.5', 'EFS = 2', 'EFS = 5') & 
#             model_name %in% c('spliner.clr', 'spliner.relab', 'metasplines')))

df_combined <- df_combined %>%
  filter(!(EFS %in% c('EFS = 1.5', 'EFS = 5') & 
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
#efs_values <- c(125, 15, 2, 3, 5)
#efs_labels <- c("EFS = 1.25", "EFS = 1.5", "EFS = 2", "EFS = 3", "EFS = 5")

efs_values <- c(15, 3, 5)
efs_labels <- c("EFS = 1.5", "EFS = 3", "EFS = 5")

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

# Create a combined factor for Metric and Sampling to control row order
df_combined <- df_combined %>%
  mutate(Metric_Sampling = paste(Metric, " \n", Sampling, sep=""))

# Create the main plot
p <- ggplot(df_combined, aes(x = sample_size, y = Mean, color = model_name, group = model_name)) +
  geom_line(size = 1.2, alpha = 0.8) +
  facet_grid(Metric_Sampling ~ EFS, scales = "free_y") +  # Rows: Metric_Sampling, Columns: EFS
  theme_minimal() +
  theme(
    # Main plot styling
    plot.title = element_text(
      size = 18, 
      face = "bold", 
      hjust = 0.5,
      margin = margin(b = 10)
    ),
    plot.subtitle = element_text(
      size = 14, 
      hjust = 0.5, 
      color = "gray40",
      margin = margin(b = 20)
    ),
    plot.caption = element_text(
      size = 10, 
      color = "gray50",
      hjust = 0,
      margin = margin(t = 10)
    ),
    plot.background = element_rect(fill = "white", color = "gray20", linewidth = 1),
    plot.margin = margin(20, 20, 20, 20),
    
    # Facet styling
    strip.background = element_rect(
      fill = "gray95", 
      color = "gray60",
      linewidth = 0.5
    ),
    strip.text = element_text(
      size = 12, 
      face = "bold",
      color = "gray20",
      margin = margin(4, 4, 4, 4)
    ),
    strip.text.x = element_text(size = 12, face = "bold"),  # EFS labels on top
    strip.text.y = element_text(size = 12, face = "bold"),  # Metric-Sampling labels on side
    
    # Panel styling
    panel.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "gray40", fill = NA, linewidth = 0.5),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.5),
    panel.grid.minor = element_line(color = "gray95", linewidth = 0.3),
    
    # Axis styling
    axis.title.x = element_text(
      size = 14, 
      face = "bold",
      margin = margin(t = 10)
    ),
    axis.title.y = element_text(
      size = 14, 
      face = "bold",
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      angle = 45, 
      hjust = 1,
      size = 11,
      color = "gray20"
    ),
    axis.text.y = element_text(
      size = 11,
      color = "gray20"
    ),
    axis.ticks = element_line(color = "gray40", linewidth = 0.5),
    axis.ticks.length = unit(5, "pt"),
    
    # Legend styling
    legend.position = "bottom",
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 11),
    legend.background = element_rect(fill = "gray98", color = "gray60", linewidth = 0.5),
    legend.margin = margin(10, 10, 10, 10),
    legend.key.size = unit(1.2, "cm"),
    legend.key = element_rect(fill = "gray98", color = NA)
  ) +
  labs(
    title = "Model Performance Comparison Across Sample Sizes",
    subtitle = "Performance metrics by sampling method and feature selection approach",
    x = "Sample Size (n)",
    y = "Performance Value",
    color = "Model Type",
  ) +
  scale_color_manual(values = color_palette) +
  scale_x_continuous(
    breaks = scales::pretty_breaks(n = 6),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 5),
    expand = expansion(mult = c(0.02, 0.02))
  )

print(p)

# Save the plot
ggsave(
  p,
  file = "/Users/zkarwowska/Desktop/all_metrics_effect_sizes_sampling_comparison.pdf",
  width = 10,
  height = 15,  # Increased height to accommodate more rows
  dpi = 300,
  limitsize = FALSE
)

