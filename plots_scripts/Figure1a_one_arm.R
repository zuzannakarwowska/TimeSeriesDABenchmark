library(dplyr)
library(ggplot2)
library(gridExtra)
library(glue)
# Define colors manually
color_palette <- c(
  "linda" = "#E41A1C",  # Red
  "maaslin.clr" = "#0358f1",  # Lighter Blue
  "maaslin.log_relab" = "#70d6ff",  # Darker Blue1F5B99
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
sampling_frequency <- 6

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

# Function to process data for a given effect size
process_effect_size <- function(efs_value, efs_label) {
  df_recall <- all_data %>%
    filter(effect_size == efs_value, sampling == sampling_frequency) %>%
    group_by(sample_size, model_name) %>%
    summarise(Mean = mean(recall, na.rm = TRUE), .groups = "drop") %>%
    mutate(Metric = "Recall", EFS = efs_label)
  
  df_precision <- all_data %>%
    filter(effect_size == efs_value, sampling == sampling_frequency) %>%
    group_by(sample_size, model_name) %>%
    summarise(Mean = mean(precision, na.rm = TRUE), .groups = "drop") %>%
    mutate(Metric = "Precision", EFS = efs_label)
  
  df_auc <- all_data %>%
    filter(effect_size == efs_value, sampling == sampling_frequency) %>%
    group_by(sample_size, model_name) %>%
    summarise(Mean = mean(auc, na.rm = TRUE), .groups = "drop") %>%
    mutate(Metric = "AUC", EFS = efs_label)
  
  return(bind_rows(df_recall, df_precision, df_auc))
}

# Apply function to all effect sizes and combine results
df_combined <- bind_rows(mapply(process_effect_size, efs_values, efs_labels, SIMPLIFY = FALSE))

# Reorder the metrics to Recall, Precision, then AUC
df_combined$Metric <- factor(df_combined$Metric, levels = c("Recall", "Precision", "AUC"))

df_combined$model_name <- gsub(".none", "", df_combined$model_name)
df_combined <- df_combined %>% filter(!model_name %in% c('baseline', 'baseline.clr', 'baseline.log', 'baseline.log_relab')) %>%
  filter((model_name == 'spliner.clr' & EFS == 'EFS = 3') | model_name != 'spliner.clr') %>%
  filter((model_name == 'spliner.log_relab' & EFS == 'EFS = 3') | model_name != 'spliner.log_relab')# %>%
  #filter(sample_size <=50)

df_combined <- df_combined %>%
  filter(!(EFS %in% c('EFS = 1.25', 'EFS = 1.5', 'EFS = 2', 'EFS = 5') & model_name %in% c('spliner.clr', 'spliner.relab', 'metasplines')))

# Define the y-limits dynamically in the data
df_combined <- df_combined %>%
  mutate(y_min = ifelse(Metric == "AUC", 0.5, 0),
         y_max = ifelse(Metric == "AUC", 1.1, 1.1))

# Create the plot
plot <- ggplot(df_combined, aes(x = sample_size, y = Mean, color = model_name, group = model_name)) +
  geom_line(size = 1.2) +
  facet_grid(Metric ~ EFS, scales = "free_y") +  # Allow different y scales
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 1) # Add black borders
  ) +
  labs(x = "Sample Size", y = "Mean Score", color = "Model") +
  scale_color_manual(values = color_palette) + 
  theme(strip.text = element_text(size = 12, face = "bold"))

ggsave(
  plot, 
  file = glue("/Users/zkarwowska/Desktop/EMBL_project/zeevi_dataset_v5/effect_size_m_{sampling_frequency}.pdf"), 
  dpi = 300
)

print(plot)
