library(dplyr)
library(pROC)
library(ggplot2)


dir <- '/Users/zkarwowska/Desktop/EMBL_project/zeevi_dataset_v5/results/results_files/'

file_list <- list.files(path = dir, pattern = "*.csv", full.names = TRUE)

all_data <- do.call(rbind, lapply(file_list, read.csv))

############ PLOT ############  
# PRECISION

summarized_data <- all_data %>% 
  dplyr::filter(effect_size == 3) %>% 
  dplyr::filter(sampling %in% c(2, 4, 6)) %>% 
  group_by(sampling, sample_size, model, effect_size) %>% 
  dplyr::summarize(Mean = mean(precision, na.rm=TRUE)) 

options(repr.plot.width =9, repr.plot.height =3) 

summarized_data %>% 
  ggplot(aes(x = sample_size, y = Mean, color = model)) +
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
    y = "Mean Precision"
  ) +
  ylim(0, 1)

# RECALL
summarized_data <- all_data %>% 
  dplyr::filter(effect_size == 3) %>% 
  dplyr::filter(sampling %in% c(2, 4, 6)) %>% 
  group_by(sampling, sample_size, model) %>% 
  dplyr::summarize(Mean = mean(recall, na.rm=TRUE)) 

options(repr.plot.width =9, repr.plot.height =3) 

summarized_data %>% 
  ggplot(aes(x = sample_size, y = Mean, color = model)) +
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
summarized_data <- all_data %>% 
  dplyr::filter(effect_size == 2) %>% 
  dplyr::filter(sampling %in% c(2, 4, 6)) %>% 
  group_by(sampling, sample_size, model) %>% 
  dplyr::summarize(Mean = mean(auc, na.rm=TRUE)) 

options(repr.plot.width =9, repr.plot.height =3) 

summarized_data %>% 
  ggplot(aes(x = sample_size, y = Mean, color = model)) +
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
    title = "AUC Across Different Sample Sizes for effect size = 2",
    x = "Sample Size", 
    y = "Mean AUC"
  ) +
  ylim(0, 1)
