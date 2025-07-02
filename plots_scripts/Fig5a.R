library(ggplot2)
library(dplyr)
library(readr)
library(cowplot)
library(ggplot2)
library(dplyr)

# Read data
d <- '/Users/zkarwowska/new_embl_folder/zeevi_dataset_v5/results/results_for_placebo_plots/'

LINDA_100 <- read_csv(paste0(d, 'LINDA_100.csv'))
MAASLIN_100 <- read_csv(paste0(d, 'MAASLIN_100.csv'))
ZIGMM_100 <- read_csv(paste0(d, 'ZIGMM_100.csv'))
LINDA_200 <- read_csv(paste0(d, 'LINDA_200.csv'))
MAASLIN_200 <- read_csv(paste0(d, 'MAASLIN_200.csv'))
ZIGMM_200 <- read_csv(paste0(d, 'ZIGMM_200.csv'))

plot_bars_with_threshold <- function(data, title = "", add_legend = FALSE) {
  library(ggplot2)
  library(dplyr)
  
  # Prepare M column as factor
  data <- data %>%
    mutate(
      M = ifelse(N == "baseline", "baseline", N),
      M = factor(M, levels = unique(M)),  # Preserve order
      setup = factor(setup, levels = c("baseline", "one arm", "placebo", "cross sectional"))
    )
  
  dodge <- position_dodge(width = 0.8)
  
  ggplot(data, aes(x = M, group = setup)) +
    # Grey background bars (no legend)
    geom_col(aes(y = ypred), fill = "grey30", width = 0.7, alpha = 0.5, position = dodge, show.legend = FALSE) +
    
    # Foreground TP colored bars (legend shown)
    geom_col(aes(y = tp, fill = setup), width = 0.7, position = dodge) +
    
    # Optional threshold
    geom_hline(yintercept = 40, color = "red", linetype = "dashed", linewidth = 0.5) +
    
    labs(title = title, x = "Sample size (M)", y = "") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 18),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      legend.position = ifelse(add_legend, "right", "none"),
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 12)
    ) +
    
    scale_fill_manual(
      name = "Setup",
      values = c(
        "baseline" = "#2ca02c",         # green
        "one arm" = "#1f77b4",          # blue
        "placebo" = "#ff7f0e",          # orange
        "cross sectional" = "#d62728"   # red
      )
    )
}

# Create all plots (legend only on one)
p1 <- plot_bars_with_threshold(LINDA_100, "LINDA 100")
p2 <- plot_bars_with_threshold(LINDA_200, "LINDA 200")
p3 <- plot_bars_with_threshold(MAASLIN_100, "MaAsLin2 100")
p4 <- plot_bars_with_threshold(MAASLIN_200, "MaAsLin2 200")
p5 <- plot_bars_with_threshold(ZIGMM_100, "ZIGMM 100")
p6 <- plot_bars_with_threshold(ZIGMM_200, "ZIGMM 200", add_legend = TRUE)

# Extract legend from one plot
legend <- get_legend(p6)

# Combine plots without legend
plot_grid_no_legend <- plot_grid(
  p1, p2, p3, p4, p5, p6 + theme(legend.position = "none"),
  ncol = 2, nrow = 3
)

# Add legend on the right
final_plot <- plot_grid(plot_grid_no_legend, legend, ncol = 2, rel_widths = c(0.85, 0.15))
print(final_plot)

# Save
ggsave('/Users/zkarwowska/Desktop/Fig5a.pdf',
       final_plot, width = 7, height = 10, dpi = 300)
