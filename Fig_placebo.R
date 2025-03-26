library(dplyr)
library(ggplot2)
library(gridExtra)

color_palette <- c(
  "linda." = "#E41A1C",  # Red
  "maaslin.clr" = "#0358f1",  # Lighter Blue
  "maaslin.log_relab" = "#70d6ff",  # Darker Blue1F5B99
  "metasplines" = "#4DAF4A",  # Green
  "trendyspliner.clr" = "#ff6d00",  # Lighter Yellow
  "trendyspliner.log_relab" = "#ffd400",  # Darker Yellow
  "wilcoxon.clr" = "#9d4edd",  # Lighter Brown
  "wilcoxon.log_relab" = "#fb6f92",  # Darker Brown
  "nbmm" = "#A65628",  # Purple
  "nbzimm" = "#999999",  # Orange
  "zigmm." = "#F781BF"   # Pink
)

effect_size <- 5
DF <- read.csv(glue('/Users/zkarwowska/Desktop/EMBL_project/zeevi_dataset_v5/one_arm_scripts/placebo_figure_input_efs{effect_size}.csv'))

plot_model <- function(m, legend = FALSE) {
  # Filter data for specific model and baseline
  df <- DF %>%
    filter(model_name %in% c(m, 'baseline.log_relab')) %>%
    mutate(
      hue = if_else(model_name == m, 'longitudinal', 'cross-sectional'),
      setup = factor(setup)
    ) %>%
    group_by(setup, hue, total_N) %>%
    summarise(`Recall` = mean(`Recall`), .groups = "drop")
  
  # Color dictionary with black colors
  model_color <- unname(color_palette[m])
  # Color dictionary
  color_dict <- c(
    'baseline' = 'black', 
    'one_arm' = model_color, 
    'placebo' = model_color, 
    m = model_color
  )
  
  # Create plot
  p <- ggplot(df, aes(x = total_N, y = `Recall`, color = setup, linetype = setup)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_color_manual(values = color_dict) +
    scale_linetype_manual(values = c('baseline' = 'dashed', 'one_arm' = 'solid', 'placebo' = 'dotted')) +
    ylim(0., 1) +
    labs(title = m, x = "N samples", y = "Recall") +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 14, face = "bold"),
      panel.grid = element_line(linewidth = 0.2),
      legend.position = "none",  
      plot.title = element_text(hjust = 0.5, face="bold"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    ) 
  
  # Conditionally add legend
  if (legend) {
    p <- p + theme(
      legend.position = "right",
      legend.box = "vertical",
      legend.margin = margin(0, 10, 10, 10)
    )
  }
  
  return(p)
}

# Create plots
p1 <- plot_model('linda.')
p2 <- plot_model('maaslin.clr')
p3 <- plot_model('zigmm.')
p4 <- plot_model('wilcoxon.clr', legend = TRUE)  # Only last plot gets legend

# Combine plots
library(patchwork)

final_plot <- p1 + p2 + p3 + p4 + 
  plot_layout(
    ncol = 4, 
    axis_titles = "collect",
    guides = "collect"  
  ) & 
  theme(legend.position = "right")


# Save the plot
ggsave(
  final_plot, 
  file = glue("/Users/zkarwowska/Desktop/EMBL_project/zeevi_dataset_v5/placebo_efs{effect_size}.pdf"), 
  width = 15, 
  height = 4, 
  dpi = 300
)

# Print the plot
print(final_plot)
