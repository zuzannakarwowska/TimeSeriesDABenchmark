library(ggplot2)
library(dplyr)
library(gridExtra)
library(readr)

# Read data
d <- '/Users/zkarwowska/new_embl_folder/zeevi_dataset_v5/results/results_for_placebo_plots/'
LINDA_100 <- read_csv(paste0(d, 'LINDA_100.csv'))
MAASLIN_100 <- read_csv(paste0(d, 'MAASLIN_100.csv'))
ZIGMM_100 <- read_csv(paste0(d, 'ZIGMM_100.csv'))
LINDA_200 <- read_csv(paste0(d, 'LINDA_200.csv'))
MAASLIN_200 <- read_csv(paste0(d, 'MAASLIN_200.csv'))
ZIGMM_200 <- read_csv(paste0(d, 'ZIGMM_200.csv'))

# Function to create bar plots with threshold
plot_bars_with_threshold <- function(data, title = "", add_legend = FALSE) {
  
  p <- ggplot(data, aes(x = factor(N))) +
    # Background bars (ypred) in grey
    geom_col(aes(y = ypred, fill = "Background"), 
             position = position_dodge(width = 0.8), width = 0.7, alpha = 0.7) +
    # Foreground bars (tp) with setup colors
    geom_col(aes(y = tp, fill = setup), 
             position = position_dodge(width = 0.8), width = 0.7) +
    # Threshold line
    geom_hline(yintercept = 40, color = "red", linetype = "dashed", linewidth = 0.5) +
    labs(title = title, x = "", y = "") +
    theme_minimal() +
    theme(panel.grid.minor = element_blank()) +
    scale_fill_manual(
      name = "Setup",
      values = c(
        "one arm" = "#1f77b4",
        "placebo" = "#ff7f0e", 
        "baseline" = "#2ca02c",
        "cross sectional" = "#d62728",
        "Background" = "grey"
      ),
      breaks = c("one arm", "placebo", "baseline", "cross sectional")  # Exclude background from legend
    )
  
  if (!add_legend) {
    p <- p + theme(legend.position = "none")
  } else {
    p <- p + theme(legend.position = "right")
  }
  
  return(p)
}

# Create individual plots
p1 <- plot_bars_with_threshold(LINDA_100, "LINDA 100")
p2 <- plot_bars_with_threshold(LINDA_200, "LINDA 200", add_legend = TRUE)
p3 <- plot_bars_with_threshold(MAASLIN_100, "MaAsLin2 100")
p4 <- plot_bars_with_threshold(MAASLIN_200, "MaAsLin2 200")
p5 <- plot_bars_with_threshold(ZIGMM_100, "ZIGMM 100")
p6 <- plot_bars_with_threshold(ZIGMM_200, "ZIGMM 200")

# Arrange plots in a 3x2 grid
final_plot <- grid.arrange(p1, p2, p3, p4, p5, p6, 
                           ncol = 2, nrow = 3)

# Save the plot
ggsave('/Users/zkarwowska/Desktop/placebo_metrics_barplot.png', 
       final_plot, width = 7, height = 10, dpi = 300)