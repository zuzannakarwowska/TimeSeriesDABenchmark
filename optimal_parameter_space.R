library(plotly)
library(reshape2)
library(dplyr)
library(tidyverse)
library(glue)

model <- 'LinDA'
df <- read.csv('/Users/zkarwowska/new_embl_folder/zeevi_dataset_v5/results/results_files/one_arm/linda.csv')
df <- df %>% filter(effect_size == 2)

df <- df %>% 
  rename(
    "N_samples"= sample_size,
    "sampling_density" = sampling,
    "Recall" = recall
  )

df_mean <- df %>%
  group_by(N_samples, sampling_density) %>%
  summarise(across(all_of("Recall"), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

recall_matrix <- df_mean %>%
  pivot_wider(
    names_from = sampling_density, 
    values_from = Recall,
    values_fill = 0
  ) %>%
  arrange(N_samples)  # Ensure correct order

color_palette <- colorRampPalette(c("white", "#1d73ff"))

z <- as.matrix(select(recall_matrix, -N_samples))
x <- recall_matrix$N_samples
y <- unique(df_mean$sampling_density)


# Create the plot
plot_obj <- persp(
  x = x, 
  y = y, 
  z = z, 
  col = color_palette(10),
  theta = 40,    
  phi = 30,      
  shade = 0.1,   
  main = model,
  xlab = "N-samples",
  ylab = "Sampling Density", 
  zlab = "Recall",
  border = "black",
  box = TRUE,
  ltheta = 120,  
  lphi = 100,     
  ticktype = "detailed"
)
