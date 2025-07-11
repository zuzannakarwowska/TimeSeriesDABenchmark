library(dplyr)
library(pROC)
library(ggplot2)
library(vegan)
library(patchwork)
library(glue)
library(ggrepel)
library(tidyverse)

# Define colors manually
color_palette <- c(
  "linda" = "#E41A1C",  # Red
  "maaslin.clr" = "#0358f1",  # Lighter Blue
  "maaslin.log" = "#70d6ff",  # Darker Blue1F5B99
  "metasplines" = "#4DAF4A",  # Green
  "trendyspliner.clr" = "#ff6d00",  # Lighter Yellow
  "trendyspliner.log" = "#ffd400",  # Darker Yellow
  "wilcoxon.clr" = "#9d4edd",  # Lighter Brown
  "wilcoxon.log" = "#fb6f92",  # Darker Brown
  "nbmm" = "#A65628",  # Purple
  "nbzimm" = "#999999",  # Orange
  "zigmm" = "#F781BF"   # Pink
)



wd <- '/Users/zkarwowska/new_embl_folder/zeevi_dataset_v5/simulation/'

treshold3 <- read.table(file = glue(wd,
                                    'efs3/rep3/gt_features.tsv'), sep = '\t',
                        header = TRUE) %>% 
  filter(upregulated != 0) %>% nrow()


efs <- '3'; sample_size= 50; sampling = 4
calculate_tp <- function(dir, model_name, sample_size, sampling) {
  
  final_df <- data.frame()
  
  for (rep in 1:10) {
    
    pattern <- glue(".*_sampling_{sample_size}_efs_{efs}_rep_{rep}_d_{sampling}.csv")
    
    files <- list.files(path = dir, pattern = pattern, full.names = TRUE)
    
    for (file_path in files) {
      df <- read.csv(file_path)

      if ('feature_assigned' %in% colnames(df)) {
        df$feature <- df$feature_assigned
        df$feature_assigned <- NULL  # Alternatively, set to NULL to remove it
      }
      if ('pval_results' %in% colnames(df)){
        
        df <- df %>% rename(pvalue = pval_results)
        
      }
      if ('pval' %in% colnames(df)){
        df <- df %>% rename(pvalue = pval)
      }
      if ('iteration' %in% colnames(df)) {
        # Column exists, proceed with your operations
      } else {
        df$iteration <- rep
        
      }

      
      df <- df %>% dplyr::select(feature, pvalue, upregulated, iteration)
      
      df$q.value <- p.adjust(df$pvalue, method = "BH")
      df$ytrue <- abs(df$upregulated)
      
      df <- df %>%
        mutate(ypred = case_when(
          (model_name %in% c('spliner_clr', 'spliner_log', 'metasplines') & q.value <= 0.05) ~ 1,#pvalue
          (!model_name %in% c('spliner_clr', 'spliner_log', 'metasplines') & q.value < 0.05) ~ 1,
          TRUE ~ 0
        ))
      
      #df <- df %>%
        #mutate(ypred = ifelse((q.value < 0.05), 1, 0))
      
      names(df)[names(df) == "ypred"] <- glue("{model_name}_ypred")
      
      df <- df %>%
        mutate(ypred = ifelse(q.value <= 0.05, 1, 0)) %>%
        mutate(!!glue("{model_name}_tp") := ifelse(ypred == 1 & ytrue == 1, 1, 0))
      
      final_df <- bind_rows(final_df, df) %>% 
        select(feature, !!glue("{model_name}_ypred"), !!glue("{model_name}_tp"), iteration)
    }
  }
  
  return(final_df)
}

results_dir <- '/Users/zkarwowska/new_embl_folder/zeevi_dataset_v5/results/'

linda_dir <- glue("{results_dir}linda/one_arm/")
maaslin_clr_dir <- glue("{results_dir}maaslin/one_arm/clr/")
maaslin_log_dir <- glue("{results_dir}maaslin/one_arm/log/")
nbmm_dir <- glue("{results_dir}nbmm/one_arm/")
nbzimm_dir <- glue("{results_dir}nbzimm/one_arm/")
zinb_dir <- glue("{results_dir}zigmm/one_arm/")
metasplines_dir <- glue("/Users/zkarwowska/new_embl_folder/zeevi_dataset_v5/results/metasplines/one_arm/")
spliner_clr_dir <- glue("/Users/zkarwowska/new_embl_folder/zeevi_dataset_v5/results/spliner/one_arm/clr/")
spliner_log_dir <- glue("/Users/zkarwowska/new_embl_folder/zeevi_dataset_v5/results/spliner/one_arm/relab/")
wilcoxon_clr_dir <- ("/Users/zkarwowska/Desktop/EMBL_project/zeevi_dataset_v5/results/wilcoxon/one_arm_all/clr/")
wilcoxon_log_dir <- ("/Users/zkarwowska/Desktop/EMBL_project/zeevi_dataset_v5/results/wilcoxon/one_arm_all/log_relab/")

linda_tp <- calculate_tp(linda_dir, 'linda', sample_size, sampling)
maaslin_clr_tp <- calculate_tp(maaslin_clr_dir, 'maaslin_clr', sample_size, sampling)
maaslin_log_tp <- calculate_tp(maaslin_log_dir, 'maaslin_log', sample_size, sampling)
wilcoxon_clr_tp <- calculate_tp(wilcoxon_clr_dir, 'wilcoxon_clr', sample_size, sampling)
wilcoxon_log_tp <- calculate_tp(wilcoxon_log_dir, 'wilcoxon_log', sample_size, sampling)
nbmm_tp <- calculate_tp(nbmm_dir, 'nbmm', sample_size, sampling)
nbzimm_tp <- calculate_tp(nbzimm_dir, 'nbzimm', sample_size, sampling)
zinb_tp <- calculate_tp(zinb_dir, 'zinb', sample_size, sampling)
#spliner_clr_tp <- calculate_tp(spliner_clr_dir, 'spliner_clr', sample_size, sampling)
#spliner_log_tp <- calculate_tp(spliner_log_dir, 'spliner_log', sample_size, sampling)
#metasplines_tp <- calculate_tp(metasplines_dir, 'metasplines', sample_size, sampling)


# Combine all results
dfs <- list(linda_tp, 
            maaslin_clr_tp, maaslin_log_tp, 
            zinb_tp, 
            wilcoxon_clr_tp, wilcoxon_log_tp, 
            #spliner_clr_tp, spliner_log_tp, 
            #metasplines_tp, 
            nbmm_tp, 
            nbzimm_tp)

combined_data <- Reduce(function(x, y) merge(x, y, by = c("feature", "iteration")), dfs) %>% select(-c(feature, iteration))

### ### ### ###  barplot

# Function to process each dataframe and extract mean TP
calculate_model_metrics <- function(df_input, model_name) {
  # Calculate TP metrics
  tp_results <- df_input %>% 
    select(matches(paste0(model_name, "_tp")), iteration) %>%
    filter(iteration %in% 1:10) %>%
    group_by(iteration) %>%
    summarise(Mean_TP = sum(get(paste0(model_name, "_tp")), na.rm = TRUE), .groups = "drop") %>%
    summarise(mean_tp = mean(Mean_TP))
  
  # Calculate YPred metrics
  ypred_results <- df_input %>% 
    select(matches(paste0(model_name, "_ypred")), iteration) %>%
    filter(iteration %in% 1:10) %>%
    group_by(iteration) %>%
    summarise(Mean_YPred = sum(get(paste0(model_name, "_ypred")), na.rm = TRUE), .groups = "drop") %>%
    summarise(mean_ypred = mean(Mean_YPred))
  
  # Combine results
  return(tibble(
    model = model_name,
    mean_tp = tp_results$mean_tp,
    mean_ypred = ypred_results$mean_ypred
  ))
}

# Apply to all your dataframes and combine results
results_df <- bind_rows(
  calculate_model_metrics(linda_tp, "linda"),
  calculate_model_metrics(maaslin_clr_tp, "maaslin_clr"),
  calculate_model_metrics(maaslin_log_tp, "maaslin_log"),
  calculate_model_metrics(wilcoxon_clr_tp, "wilcoxon_clr"),
  calculate_model_metrics(wilcoxon_log_tp, "wilcoxon_log"),
  calculate_model_metrics(nbmm_tp, "nbmm"),
  calculate_model_metrics(nbzimm_tp, "nbzimm"),
  calculate_model_metrics(zinb_tp, "zinb"),
  #calculate_model_metrics(spliner_clr_tp, "spliner_clr"),
  #calculate_model_metrics(spliner_log_tp, "spliner_log"),
  #calculate_model_metrics(metasplines_tp, "metasplines"),
)

results_df <- results_df %>%
  mutate(model = model %>%
           str_replace_all("_tp", "") %>%
           str_replace_all("_", ".") %>%
           str_replace_all("zinb", "zigmm") %>%
           str_replace_all("spliner", "trendyspliner"))


results_long <- results_df %>%
  tidyr::pivot_longer(
    cols = c(mean_ypred, mean_tp),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    fill_color = case_when(
      metric == "mean_tp" ~ model,  # Use model name as fill group for TP
      TRUE ~ "YPred"                # Use "YPred" as fill group for YPred
    )
  )


ggplot(results_long, aes(y = model, x = value, fill = fill_color)) +
  geom_bar(stat = "identity", position = "identity", width = 0.7, color = "white", alpha = 0.7) +
  scale_fill_manual(
    values = c("YPred" = "grey", color_palette),
    name = "Metric") +
  labs(x = "", y = "") +
  #xlim(0, 50) + 
  theme_minimal() +
  geom_vline(xintercept = treshold3, color = "black", size = 0.9, linetype = "dotted") + 
  theme(
    axis.text.y = element_text(face = "bold", size = 14),  # Change y-axis text size
    axis.text.x = element_text(size = 14),  # Change x-axis text size
    axis.title.x = element_text(size = 14),  # Change x-axis title size
    axis.title.y = element_text(size = 14),  # Change y-axis title size
    plot.title = element_text(hjust = 0.5, size = 14),  # Change title size
    legend.text = element_text(size = 14),  # Change legend text size
    legend.title = element_text(size = 14),  # Change legend title size
    legend.position = "none",
    panel.grid.major.y = element_blank()
  )
