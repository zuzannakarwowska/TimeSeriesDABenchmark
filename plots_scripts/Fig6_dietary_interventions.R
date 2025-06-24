library(ggplot2)
library(ggrepel)
library(ragg)
.libPaths(c("/g/scb/zeller/karwowsk/R/4.4.1-gfbf-2023b/", .libPaths()))
library(dplyr)
library(glue)
library(tibble)
library(eulerr)
library(ggplot2)
library(tidyr)
library(ragg)

setwd('/Users/zkarwowska/new_embl_folder/zeevi_dataset_v5/real_datasets/')

n_palleja <- 80
palleja_res <- read.csv('dietary_interventions_results_per_dataset.csv', row.names = 'X') %>% 
  filter(study != 'baxter')

# Create sample_size_study combination for x-axis labels
palleja_res <- palleja_res %>%
  mutate(sample_size_study = paste0(study, "dataset \n N=", sample_size))

# Instead of averaging, we'll work with individual dataset entries
# But still normalize by n_palleja
new_df <- palleja_res %>% 
  mutate(across(where(is.numeric) & !matches("sample_size"), ~ .x / n_palleja))

results_long <- new_df %>%
  pivot_longer(
    cols = -c(sample_size, study, sample_size_study),
    names_to = "method",
    values_to = "n_features"
  ) %>%
  mutate(
    linetype = case_when(
      method %in% c("linda_n", "maaslin_n", "zigmm_n") ~ "dotted",
      TRUE ~ "solid"
    ),
    alpha = if_else(method == "all_three_overlap", 1, 0.7),
    size = if_else(method == "all_three_overlap", 1.5, 1),
    point_size = if_else(method == "all_three_overlap", 4, 3),
    stroke = if_else(method == "all_three_overlap", 1.2, 0.8)
  )

# Create label data for the last point of each method
label_data <- results_long %>%
  group_by(method) %>%
  filter(sample_size == max(sample_size)) %>%
  slice(1) %>%  # Take just one entry per method if there are multiple datasets with max sample size
  ungroup()

# Create custom x-axis positions
unique_combinations <- results_long %>%
  distinct(sample_size, study, sample_size_study) %>%
  arrange(sample_size, study) %>%
  mutate(x_position = row_number())

# Join back to get x_positions
results_long <- results_long %>%
  left_join(unique_combinations, by = c("sample_size", "study", "sample_size_study"))

label_data <- label_data %>%
  left_join(unique_combinations, by = c("sample_size", "study", "sample_size_study"))

g <- ggplot(results_long, aes(
  x = x_position,
  y = n_features,
  color = method,
  linetype = linetype,
  alpha = alpha,
  size = size)) +
  geom_line() +
  geom_point(aes(size = point_size, stroke = stroke), shape = 21, fill = "white") +
  ggrepel::geom_text_repel(
    data = label_data,
    aes(label = method, color = method),
    hjust = -0.5,
    vjust = -0.1,
    size = 5,
    fontface = "bold",
    show.legend = FALSE,
    direction = "y",
    nudge_x = -0.01,
    box.padding = 0.1,
    point.padding = 0.05,
    segment.color = "grey50",
    segment.size = 0.2,
    min.segment.length = 0.0,
    max.overlaps = Inf,
    force = 1,
    force_pull = 0.5
  ) +
  scale_linetype_identity() +
  scale_alpha_identity() +
  scale_size_identity() +
  scale_x_continuous(
    breaks = unique_combinations$x_position,
    labels = unique_combinations$sample_size_study,
    expand = expansion(mult = c(0.02, 0.8))
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal() +
  theme(
    plot.margin = margin(5.5, 60, 5.5, 5.5),
    plot.title = element_text(
      size = 18,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 5)
    ),
    plot.subtitle = element_text(
      size = 14,
      hjust = 0.5,
      color = "gray40",
      margin = margin(b = 15)
    ),
    plot.caption = element_text(
      size = 10,
      color = "gray50",
      hjust = 1,
      margin = margin(t = 10)
    ),
    axis.title.x = element_text(
      size = 16,
      face = "bold",
      margin = margin(t = 15)  # Increased margin to accommodate rotated labels
    ),
    axis.title.y = element_text(
      size = 16,
      face = "bold",
      margin = margin(r = 10)
    ),
    axis.line.x.bottom = element_line(color = "black", linewidth = 0.8),
    axis.line.y.left = element_line(color = "black", linewidth = 0.8),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # Rotated x-axis labels
    axis.text.y = element_text(size = 14),
    axis.ticks.length.x = unit(10, "pt"),
    axis.ticks.length.y = unit(5, "pt"),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  labs(
    title = "Feature Selection Performance Across Sample Sizes and Datasets",
    subtitle = "Comparison of different statistical methods for identifying significant features",
    x = "Sample Size_Dataset",
    y = "Number of Significant Features Detected",
    color = "Method"
  )

print(g)

# Save as PDF with adjusted dimensions to accommodate rotated labels
ggsave(
  filename = "/Users/zkarwowska/Desktop/dietary_lineplot_with_datasets.pdf",
  plot = g,
  width = 10,  # Increased width to accommodate more x-axis labels
  height = 5,  # Increased height for rotated labels
  units = "in",
  dpi = 300,
  device = "pdf"
)


full_dataset <- read.csv('dietary_interventions_results_full.csv', row.names = 'X')

# Calculate overlaps
fit <- euler(c(
  "LinDA" = full_dataset$linda_n,
  "MaAsLin2" = full_dataset$maaslin_n,
  "ZIGMM" = full_dataset$zigmm_n,
  "LinDA&MaAsLin2" = full_dataset$linda_maaslin_n,
  "LinDA&ZIGMM" = full_dataset$linda_zigmm_n,
  "MaAsLin2&ZIGMM" = full_dataset$maaslin_zigmm_n,
  "LinDA&MaAsLin2&ZIGMM" = full_dataset$all_three_overlap
))

# Plot

g <- plot(fit,
          fills = list(fill = c("skyblue", "pink", "lightgreen"), alpha = 0.6),
          edges = list(col = "black"),
          labels = list(font = 10, cex = 2),        # Adjust label size
          quantities = list(cex = 2))

print(g)


ggsave(
  filename = "/Users/zkarwowska/Desktop/dietary_interventions_venn.pdf",
  plot = g,
  width = 7,
  height = 7,
  units = "in",
  dpi = 300,
  device = "pdf"
)
