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

n_palleja <- 179
palleja_res <- read.csv('palleja_results.csv', row.names  ='X')

palleja_res_mean <- palleja_res %>%
  group_by(sample_size) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))

new_df <- palleja_res_mean %>% select(-c('sample_size', 'rep')) / n_palleja
new_df$sample_size <- palleja_res_mean$sample_size
new_df$rep <- palleja_res_mean$rep

results_long <- new_df %>%
  pivot_longer(
    cols = -c(sample_size, rep),
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
  ungroup()

g <- ggplot(results_long, aes(
  x = sample_size,
  y = n_features,
  color = method,
  linetype = linetype,
  alpha = alpha,
  size = size
)) +
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
    direction = "y",  # Allow vertical adjustment
    nudge_x = -0.01,   # Push labels slightly to the left
    box.padding = 0.1,  # Reduced padding for shorter segments
    point.padding = 0.05,  # Reduced padding for shorter segments
    segment.color = "grey50",  # Color for connecting segments
    segment.size = 0.2,  # Thickness of connecting segments
    min.segment.length = 0.0,  # Minimum segment length before hiding
    max.overlaps = Inf,  # Allow all labels to be shown
    force = 1,  # Increase repulsion force to keep labels closer
    force_pull = 0.5  # Add pull force to minimize segment length
  ) +
  scale_linetype_identity() +
  scale_alpha_identity() +
  scale_size_identity() +
  scale_x_continuous(
    expand = expansion(mult = c(0.02, 0.8))  # Added more space on the right for labels
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  coord_cartesian(clip = "off") +  # This allows labels to extend beyond plot area
  theme_minimal() +
  theme(
    plot.margin = margin(5.5, 60, 5.5, 5.5),  # Increased right margin from 40 to 60
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
      margin = margin(t = 10)
    ),
    axis.title.y = element_text(
      size = 16, 
      face = "bold",
      margin = margin(r = 10)
    ),
    axis.line.x.bottom = element_line(color = "black", linewidth = 0.8),
    axis.line.y.left = element_line(color = "black", linewidth = 0.8),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.ticks.length.x = unit(10, "pt"),
    axis.ticks.length.y = unit(5, "pt"),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  labs(
    title = "Feature Selection Performance Across Sample Sizes",
    subtitle = "Comparison of different statistical methods for identifying significant features",
    x = "Sample Size (n)",
    y = "Number of Significant Features Detected",
    color = "Method"
  )

print(g)

# Save as PDF with 10x5 inch dimensions
ggsave(
  filename = "/Users/zkarwowska/Desktop/palleja_lineplot.pdf",
  plot = g,
  width = 10,
  height = 5,
  units = "in",
  dpi = 300,
  device = "pdf"
)


full_dataset <- palleja_res_mean %>% filter(sample_size == 12) %>% round()

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
  filename = "/Users/zkarwowska/Desktop/palleja_venn.pdf",
  plot = g,
  width = 7,
  height = 7,
  units = "in",
  dpi = 300,
  device = "pdf"
)
