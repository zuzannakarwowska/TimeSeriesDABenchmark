library(ggplot2)
library(ggrepel)
library(ragg)
library(dplyr)
library(glue)
library(tibble)
library(eulerr)
library(ggplot2)
library(tidyr)
library(ragg)
library(patchwork)

# Define colors manually
color_palette <- c(
  "linda_n" = "#E41A1C",  # Red
  "maaslin_n" = "#0358f1",  # Lighter Blu
  "zigmm_n" = "#F781BF"   # Pink
)



# ------------------------------------------------------------------------------
# Line Plot
# ------------------------------------------------------------------------------
plot_lineplot <- function(df, dataset){
  

  new_df <- df  %>%
    select(c("linda_n", "maaslin_n", "zigmm_n"))

  new_df$sample_size <- df$sample_size
  new_df$rep <- df$rep
  
  results_long <- new_df %>%
    pivot_longer(
      cols = -c(sample_size, rep),
      names_to = "method",
      values_to = "n_features"
    ) %>%
    mutate(
      linetype = case_when(
        method %in% c("linda_n", "maaslin_n", "zigmm_n") ~ "solid",
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
    scale_color_manual(values = color_palette) + 
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
      panel.grid.major.y = element_line(color = "grey80", size = 0.5),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    labs(
      x = "Sample Size (n)",
      y = "Number of Significant Features Detected",
      color = "Method"
    )
  
  print(g)
  
  # Save as PDF with adjusted dimensions to accommodate rotated labels
  ggsave(
    filename = paste0("/Users/zkarwowska/Desktop/", dataset, "_lineplot.pdf"),
    plot = g,
    width = 10,  # Increased width to accommodate more x-axis labels
    height = 5,  # Increased height for rotated labels
    units = "in",
    dpi = 300,
    device = "pdf"
  )
}

# ------------------------------------------------------------------------------
# Venn Diagram
# ------------------------------------------------------------------------------

library(ggVennDiagram)

plot_venn <- function(full_dataset, title){
  
  # Define counts (assuming full_dataset is already defined)
  linda_only <- full_dataset$linda_n -
    full_dataset$linda_maaslin_n -
    full_dataset$linda_zigmm_n +
    full_dataset$all_three_overlap
  
  maaslin_only <- full_dataset$maaslin_n -
    full_dataset$linda_maaslin_n -
    full_dataset$maaslin_zigmm_n +
    full_dataset$all_three_overlap
  
  zigmm_only <- full_dataset$zigmm_n -
    full_dataset$linda_zigmm_n -
    full_dataset$maaslin_zigmm_n +
    full_dataset$all_three_overlap
  
  linda_maaslin <- full_dataset$linda_maaslin_n - full_dataset$all_three_overlap
  linda_zigmm <- full_dataset$linda_zigmm_n - full_dataset$all_three_overlap
  maaslin_zigmm <- full_dataset$maaslin_zigmm_n - full_dataset$all_three_overlap
  all_three <- full_dataset$all_three_overlap
  
  # Create a named list of sets based on overlaps
  # You'll create dummy element names to simulate the overlaps
  linda <- rep("a", linda_only)
  maaslin <- rep("b", maaslin_only)
  zigmm <- rep("c", zigmm_only)
  
  # Add overlaps
  linda_maaslin_set <- rep("ab", linda_maaslin)
  linda_zigmm_set <- rep("ac", linda_zigmm)
  maaslin_zigmm_set <- rep("bc", maaslin_zigmm)
  all_three_set <- rep("abc", all_three)
  
  # Combine all elements into one list per set
  # We'll convert this into sets of character vectors
  elements <- c(
    rep("linda", linda_only),
    rep("maaslin", maaslin_only),
    rep("zigmm", zigmm_only),
    rep("linda_maaslin", linda_maaslin),
    rep("linda_zigmm", linda_zigmm),
    rep("maaslin_zigmm", maaslin_zigmm),
    rep("all_three", all_three)
  )
  
  sets <- list(
    LINDA = which(elements %in% c("linda", "linda_maaslin", "linda_zigmm", "all_three")),
    MAASLIN = which(elements %in% c("maaslin", "linda_maaslin", "maaslin_zigmm", "all_three")),
    ZIGMM = which(elements %in% c("zigmm", "linda_zigmm", "maaslin_zigmm", "all_three"))
  )
  
  venn.plot <- ggVennDiagram(sets, label_alpha = 0) +
    scale_fill_gradient(low = "white", high = "grey") +
    coord_fixed(ratio = 0.6) +
    scale_x_continuous(expand = expansion(mult = 0.2))
  
  print(venn.plot)
  
  # Save to file
  pdf(paste0("/Users/zkarwowska/Desktop/", title, "_venn.pdf"), width = 7, height = 7)
      grid.draw(venn.plot)
      dev.off()
}


# ------------------------------------------------------------------------------
# LogFC
# ------------------------------------------------------------------------------
plot_logFC <- function(counts, metadata, title){

  # Create named vector for mapping sample ID to timepoint
  timepoint_dict <- setNames(metadata$Timepoint, metadata$Sample_ID)
  
  # Reindex rows of counts using timepoint_dict
  counts <- counts[ , metadata$Sample_ID] %>% t() %>% as.data.frame()
  counts$Timepoint <- timepoint_dict[rownames(counts)]
  
  # Group by Timepoint and calculate mean (+ pseudocount)
  counts_mean_df <- counts %>%
    group_by(Timepoint) %>%
    summarise(across(everything(), mean)) %>%
    column_to_rownames("Timepoint") %>%
    as.matrix() + 1e-3
  
  # Calculate absolute log fold change between timepoints 0 and 4
  logFC <- log2(counts_mean_df["1", ] / counts_mean_df["0", ])
  
  # Plot using ggplot2
  logFC_df <- data.frame(logFC = logFC)
  logFC_df$logFC <- abs(logFC_df$logFC)

  g <- ggplot(logFC_df, aes(y = logFC, x = 'feature')) + 
    geom_violin(
      width = 0.3,
      alpha = 0.7,
      fill = "#4CAF50",
      color = "#4CAF50",
      outlier.color = "#4CAF50",
      outlier.size = 2
    ) +
    geom_jitter(
      aes(y = logFC),
      width = 0.25,
      alpha = 0.6,
      color = "#2E7D32",
      size = 1.5
    ) +
    labs(
      title = "Distribution of Log Fold Changes",
      subtitle = "Boxplot with individual data points",
      y = "Log Fold Change (logFC)",
      x = NULL
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5, color = "gray60"),
      axis.title.y = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 14),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 14), 
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  print(g)
  
  ggsave(
    filename = paste0("/Users/zkarwowska/Desktop/", title, "_logFC.pdf"),
    plot = g,
    width = 5,
    height = 7,
    units = "in",
    dpi = 300,
    device = "pdf"
  )
}

plot_clrFC <- function(counts, metadata, title){
  
  # Create named vector for mapping sample ID to timepoint
  timepoint_dict <- setNames(metadata$Timepoint, metadata$Sample_ID)
  
  # Reindex rows of counts using timepoint_dict
  counts <- counts[ , metadata$Sample_ID] %>% t() %>% as.data.frame()
  counts$Timepoint <- timepoint_dict[rownames(counts)]

  clr_data <- clr(t(counts + 1)) 
  clr_data <- t(clr_data)
  
  # Reindex rows of counts using timepoint_dict
  counts <- counts[ , metadata$Sample_ID] %>% t() %>% as.data.frame()
  counts$Timepoint <- timepoint_dict[rownames(counts)]
  
  # Group by Timepoint and calculate mean (+ pseudocount)
  counts_mean_df <- counts %>%
    group_by(Timepoint) %>%
    summarise(across(everything(), mean)) %>%
    column_to_rownames("Timepoint") %>%
    as.matrix() 
  
  # Calculate absolute log fold change between timepoints 0 and 4
  clrFC <- counts_mean_df["1", ] - counts_mean_df["0", ]
  
  # Plot using ggplot2
  clrFC_df <- data.frame(clrFC = clrFC)
  clrFC_df$clrFC <- abs(clrFC_df$clrFC)
  
  ###
  g <- ggplot(clrFC_df, aes(y = clrFC, x = 'feature')) + 
    geom_violin(
      width = 0.3,
      alpha = 0.7,
      fill = "#4CAF50",
      color = "#4CAF50",
      outlier.color = "#4CAF50",
      outlier.size = 2
    ) +
    geom_jitter(
      aes(y = logFC),
      width = 0.25,
      alpha = 0.6,
      color = "#2E7D32",
      size = 1.5
    ) +
    labs(
      title = "Distribution of CLR Fold Changes",
      subtitle = "Boxplot with individual data points",
      y = "CLR Fold Change (clrFC)",
      x = NULL
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5, color = "gray60"),
      axis.title.y = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 14),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 14), 
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  print(g)
  
  ggsave(
    filename = paste0("/Users/zkarwowska/Desktop/", title, "_clrFC.pdf"),
    plot = g,
    width = 5,
    height = 7,
    units = "in",
    dpi = 300,
    device = "pdf"
  )
}

# ------------------------------------------------------------------------------
# venkataraman dataset
# ------------------------------------------------------------------------------

setwd('/Users/zkarwowska/new_embl_folder/zeevi_dataset_v5/real_datasets/')

# models results
venkataraman_results <- read.csv('venkataraman_results_subsampled.csv', row.names = 'X')

venkataraman_results_mean <- venkataraman_results  %>% 
  group_by(sample_size) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))


# Plot lineplot
plot_lineplot(venkataraman_results_mean, 'venkataraman')

# Plot venn
full_dataset <- venkataraman_results_mean %>%
  filter(sample_size == 10) %>% 
  round()
plot_venn(full_dataset, 'venkataraman')

# plot LogFC
venkataraman_metadata <- read_csv("baxter_metadata.csv") %>%
  column_to_rownames(var = colnames(.)[1]) %>% filter(study != 'baxter')

venkataraman_counts <- read_csv('baxter_otu_table.csv') %>%
  column_to_rownames(var = colnames(.)[1]) %>% select(venkataraman_metadata$Sample_ID)

plot_logFC(venkataraman_counts, venkataraman_metadata, 'venkataraman')


# Load necessary library
library(compositions)

plot_clrFC(venkataraman_counts, venkataraman_metadata, 'venkataraman')

# ------------------------------------------------------------------------------
# Mesnage
# ------------------------------------------------------------------------------

mesnage_counts <- read_csv("mesnage_dataset/ready_mesnage_counts.csv") %>%
  column_to_rownames(var = colnames(.)[1])

mesnage_metadata <- read_csv("mesnage_dataset/ready_mesnage_metadata.csv") %>%
  column_to_rownames(var = colnames(.)[1])

plot_logFC(mesnage_counts, mesnage_metadata, 'mesnage')

mesnage_res <- read.csv('mesnage_results.csv', row.names  ='X')

mesnage_res_mean <- mesnage_res  %>% 
  group_by(sample_size) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))

plot_lineplot(mesnage_res_mean, 'mesnage')

mesnage_res_full <- mesnage_res_mean %>%
  filter(sample_size == 70) %>% 
  round()
plot_venn(mesnage_res_full, 'mesnage')

# ------------------------------------------------------------------------------
# Palleja
# ------------------------------------------------------------------------------

palleja_counts <- read_csv("palleja_dataset/ready_palleja_counts.csv") %>%
  column_to_rownames(var = colnames(.)[1])

palleja_metadata <- read_csv("palleja_dataset/ready_palleja_metadata.csv") %>%
  column_to_rownames(var = colnames(.)[1])

palleja_metadata$Timepoint[palleja_metadata$Timepoint %in% c(4, 8)] <- 1

plot_logFC(palleja_counts, palleja_metadata, 'palleja')

palleja_res <- read.csv('palleja_results.csv', row.names  ='X')

palleja_res_mean <- palleja_res  %>% 
  group_by(sample_size) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))

plot_lineplot(palleja_res_mean, 'palleja')

palleja_res_full <- palleja_res_mean %>%
  filter(sample_size == 12) %>% 
  round()
plot_venn(palleja_res_full, 'palleja')

