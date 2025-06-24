library(ggplot2)
library(dplyr)
library(gridExtra)
library(readr)

# Read data
grouped <- read_csv('/Users/zkarwowska/new_embl_folder/zeevi_dataset_v5/results/results_for_placebo_plots/grouped.csv')


# Function to create individual plots
plot_metric <- function(data, metric, title, show_legend = FALSE) {
  
  # Set y-axis limits
  y_limits <- if (metric == "recall") c(-0.1, 1.1) else c(0.5, 1.1)
  
  p <- ggplot(data, aes(x = factor(total), y = .data[[metric]], 
                        color = factor(N), linetype = setup, shape = factor(N))) +
    geom_point(size = 2) +
    geom_line(aes(group = interaction(N, setup)), linewidth = 0.8) +
    labs(title = title,
         x = "",
         y = stringr::str_to_title(metric)) +
    ylim(y_limits) +
    theme_minimal() +
    theme(plot.title = element_text(size = 18),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linewidth = 0.3),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14)) +
    scale_color_discrete(name = "N") +
    scale_linetype_discrete(name = "Setup") +
    scale_shape_discrete(name = "N")
  
  if (!show_legend) {
    p <- p + theme(legend.position = "none")
  } else {
    p <- p + theme(legend.position = "right",
                   legend.title = element_text(size = 14))
  }
  
  return(p)
}

# Filter data for each tool
linda_data <- grouped %>% filter(tool == "LINDA")
maaslin_data <- grouped %>% filter(tool == "MaAsLin2")
zigmm_data <- grouped %>% filter(tool == "ZIGMM")

# Create individual plots
p1 <- plot_metric(linda_data, "recall", "LINDA - Recall")
p2 <- plot_metric(linda_data, "precision", "LINDA - Precision", show_legend = TRUE)
p3 <- plot_metric(maaslin_data, "recall", "MaAsLin2 - Recall")
p4 <- plot_metric(maaslin_data, "precision", "MaAsLin2 - Precision")
p5 <- plot_metric(zigmm_data, "recall", "ZIGMM - Recall")
p6 <- plot_metric(zigmm_data, "precision", "ZIGMM - Precision")

# Arrange plots in a 3x2 grid
p <- grid.arrange(p1, p2, p3, p4, p5, p6, 
             ncol = 2, nrow = 3)
print(p)

ggsave(
  filename = "/Users/zkarwowska/Desktop/FIg5b.pdf",
  plot = p,
  width = 7,
  height = 10,
  units = "in",
  dpi = 300,
  device = "pdf"
)
