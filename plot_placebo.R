library(dplyr)
library(ggplot2)
library(patchwork)
library(glue)
library(purrr)


# Read one_arm files
one_arm_results_dir <- '/Users/zkarwowska/Desktop/EMBL_project/zeevi_dataset_v5/results/results_files/'
one_arm_pattern <- ".*one_arm.csv"

one_arm_files <- list.files(path = one_arm_results_dir, pattern = one_arm_pattern, full.names = TRUE)

one_arm_combined_data <- one_arm_files %>%
  map_df(read.csv) %>%  
  mutate(
    setup = 'one_arm',  
    model_name = paste(model, transformation, sep = "_") %>%  
      gsub("_none$|_NONE$", "", .)  
  )

# Read placebo files
placebo_results_dir <- '/Users/zkarwowska/Desktop/EMBL_project/zeevi_dataset_v5/results/results_files/'
placebo_pattern <- ".*placebo.csv"

# List all files matching the pattern
placebo_files <- list.files(path = placebo_results_dir, pattern = placebo_pattern, full.names = TRUE)

# Combine files and process the data
placebo_combined_data <- placebo_files %>%
  map_df(read.csv) %>%  # Combine all CSV files into one dataframe
  mutate(
    setup = 'placebo',  # Add the setup column
    model_name = paste(model, transformation, sep = "_") %>%  # Combine model and transformation
      gsub("_none$|_NONE$", "", .)  # Remove both "_none" and "_NONE" suffixes
  )
# Combine
combined_df <- bind_rows(one_arm_combined_data, placebo_combined_data)
combined_df[is.na(combined_df)] <- 1

#Summarise
summarized_data <- combined_df %>% 
  dplyr::filter(effect_size == 3) %>% 
  dplyr::filter(sampling %in% c(2, 4, 6)) %>% 
  group_by(sampling, sample_size, model_name, setup) %>% 
  dplyr::summarize(Mean = mean(precision, na.rm=TRUE)) 

options(repr.plot.width =9, repr.plot.height =3) 

# PRECISION
summarized_data %>% 
  ggplot(aes(x = sample_size, y = Mean, color = model_name, linetype = setup)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  facet_wrap(~sampling, ncol = 3, labeller = labeller(sampling = c('2' = 'F=2', '4' = 'F=4', '6' = 'F=6'))) +
  theme_minimal() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5)  # Center the title
  ) +
  labs(
    title = "Precision Across Different Sample Sizes for effect size = 3",
    x = "Sample Size", 
    y = "Mean Recall"
  ) +
  ylim(0, 1)

# RECALL  
summarized_data <- combined_df %>% 
  dplyr::filter(effect_size == 3) %>% 
  dplyr::filter(sampling %in% c(2, 4, 6)) %>% 
  group_by(sampling, sample_size, model_name, setup) %>% 
  dplyr::summarize(Mean = mean(recall, na.rm=TRUE)) 

options(repr.plot.width =9, repr.plot.height =3) 

summarized_data %>% 
  ggplot(aes(x = sample_size, y = Mean, color = model_name, linetype = setup)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  facet_wrap(~sampling, ncol = 3, labeller = labeller(sampling = c('2' = 'F=2', '4' = 'F=4', '6' = 'F=6'))) +
  theme_minimal() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5)  # Center the title
  ) +
  labs(
    title = "Recall Across Different Sample Sizes for effect size = 3",
    x = "Sample Size", 
    y = "Mean Recall"
  ) +
  ylim(0, 1)


## AUC
summarized_data <- combined_df %>% 
  dplyr::filter(effect_size == 3) %>% 
  dplyr::filter(sampling %in% c(2, 4, 6)) %>% 
  group_by(sampling, sample_size, model_name, setup) %>% 
  dplyr::summarize(Mean = mean(auc, na.rm=TRUE)) 

options(repr.plot.width =9, repr.plot.height =3) 

summarized_data %>% 
  ggplot(aes(x = sample_size, y = Mean, color = model_name, linetype = setup)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  facet_wrap(~sampling, ncol = 3, labeller = labeller(sampling = c('2' = 'F=2', '4' = 'F=4', '6' = 'F=6'))) +
  theme_minimal() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5)  # Center the title
  ) +
  labs(
    title = "AUC Across Different Sample Sizes for effect size = 3",
    x = "Sample Size", 
    y = "Mean Recall"
  ) +
  ylim(0, 1)
