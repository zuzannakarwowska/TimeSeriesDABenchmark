library(dplyr)
library(pROC)
library(ggplot2)
library(vegan)
library(patchwork)

calculate_tp <- function(dir, model_name, sample_size, sampling) {
  final_df <- data.frame()
  
  for (rep in 1:10) {
    pattern <- glue(".*_sampling_{sample_size}_efs_3_rep_{rep}_d_{sampling}.csv")
    
    files <- list.files(path = dir, pattern = pattern, full.names = TRUE)
    
    for (file_path in files) {
      df <- read.csv(file_path)
      
      if ('feature_assigned' %in% colnames(df)) {
        df$feature <- df$feature_assigned
        df$feature_assigned <- NULL  # Alternatively, set to NULL to remove it
      }
      
      if ('iteration' %in% colnames(df)) {
        # Column exists, proceed with your operations
      } else {
        df$iteration <- rep
      }
      
      if ('pvalue' %in% colnames(df)) {
        # Column exists, proceed with your operations
      } else {
        df$pvalue <- df$pval
      }
      
      df <- df %>% select(feature, pvalue, upregulated, iteration)
      
      df$q.value <- p.adjust(df$pvalue, method = "BH")
      df$ytrue <- abs(df$upregulated)
      
      df <- df %>%
        mutate(ypred = ifelse(q.value <= 0.05, 1, 0)) %>%
        mutate(!!glue("{model_name}_tp") := ifelse(ypred == 1 & ytrue == 1, 1, 0))
      
      final_df <- bind_rows(final_df, df) %>% select(feature, !!glue("{model_name}_tp"), iteration)
    }
  }
  
  return(final_df)
}

results_dir <- '/Users/zkarwowska/Desktop/EMBL_project/zeevi_dataset_v5/results/'

linda_dir <- glue("{results_dir}linda/one_arm/")
maaslin_clr_dir <- glue("{results_dir}maaslin/one_arm/clr/")
maaslin_log_dir <- glue("{results_dir}maaslin/one_arm/log/")
wilcoxon_clr_dir <- glue("{results_dir}wilcoxon/one_arm_all/clr/")
wilcoxon_log_dir <- glue("{results_dir}wilcoxon/one_arm_all/log_relab/")
nbmm_dir <- glue("{results_dir}nbmm/one_arm/")
nbzimm_dir <- glue("{results_dir}nbzimm/one_arm/")
zinb_dir <- glue("{results_dir}zigmm/one_arm/")
metasplines_dir <- glue("{results_dir}metasplines/one_arm/")
spliner_clr_dir <- glue("{results_dir}trendyspliner/one_arm/clr/")
spliner_log_dir <- glue("{results_dir}trendyspliner/one_arm/log/")



pcoa_plots <- list()

# Loop over different sample sizes
for (sample_size in c(10, 30, 60)){
  sampling <- 6
  linda_tp <- calculate_tp(linda_dir, 'linda', sample_size, sampling)
  maaslin_clr_tp <- calculate_tp(maaslin_clr_dir, 'maaslin_clr', sample_size, sampling)
  maaslin_log_tp <- calculate_tp(maaslin_log_dir, 'maaslin_log', sample_size, sampling)
  wilcoxon_clr_tp <- calculate_tp(wilcoxon_clr_dir, 'wilcoxon_clr', sample_size, sampling)
  wilcoxon_log_tp <- calculate_tp(wilcoxon_log_dir, 'wilcoxon_log', sample_size, sampling)
  nbmm_tp <- calculate_tp(nbmm_dir, 'nbmm', sample_size, sampling)
  nbzimm_tp <- calculate_tp(nbzimm_dir, 'nbzimm', sample_size, sampling)
  zinb_tp <- calculate_tp(zinb_dir, 'zinb', sample_size, sampling)
  
  # Combine all results
  dfs <- list(linda_tp, maaslin_clr_tp, maaslin_log_tp, wilcoxon_clr_tp, wilcoxon_log_tp, nbmm_tp, nbzimm_tp, zinb_tp)
  combined_data <- Reduce(function(x, y) merge(x, y, by = c("feature", "iteration")), dfs) %>% select(-c(feature, iteration))
  
  # Calculate Jaccard distance
  dist_matrix <- vegdist(t(combined_data), method = "jaccard", binary = TRUE)
  
  # Perform PCoA
  pcoa_result <- cmdscale(dist_matrix, k = 2, eig = TRUE)
  pcoa_scores <- data.frame(pcoa_result$points)
  colnames(pcoa_scores) <- c("PC1", "PC2")
  
  pcoa_scores$Sample <- rownames(pcoa_scores)
  
  options(repr.plot.width =9, repr.plot.height =3) 
  
  # Create the plot for the current sample size
  pcoa_plot <- ggplot(pcoa_scores, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = Sample), size = 6) +  # Increase point size
    labs(title = glue("Sample Size = {sample_size}"),
         x = "Principal Coordinate 1",
         y = "Principal Coordinate 2") +
    theme_minimal()

    if (sample_size != 60) {
    pcoa_plot <- pcoa_plot + theme(legend.position = "none")
  }
  
  # Store the plot in the list
  pcoa_plots[[as.character(sample_size)]] <- pcoa_plot
}

combined_plot <- pcoa_plots[["10"]] + pcoa_plots[["30"]] + pcoa_plots[["60"]] + plot_layout(ncol = 3)

combined_plot + 
  plot_layout(guides = 'collect') + 
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  theme(aspect.ratio = 1) 
