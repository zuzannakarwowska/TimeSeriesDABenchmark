.libPaths(c("/g/scb/zeller/karwowsk/R/4.3.3", .libPaths()))
set.seed(10)

library(dplyr)
library(data.table)
library(tidyr)
library(foreach)
library(doParallel)
library(compositions)

cores <- 3
cl <- makeCluster(cores) 
registerDoParallel(cl)
clusterCall(cl, function(x) .libPaths(x), .libPaths())
pckgs <- c('dplyr', 'data.table', 'compositions', 'foreach', 'doParallel') 


Reps <- 1:10
efs=3
Sampling <- list(c(1, 19))#list(c(1, 19), c(1, 6, 15, 19), c(1, 5, 8, 12, 15, 19))
SampleSizes <- c(10, 20, 30, 40, 50, 60)
#efs=3; rep=1; sample_size=10; sampling <- c(1, 5, 8, 12, 15, 19)
simulation_path <- paste0('/g/scb/zeller/karwowsk//zeevi_dataset_v5/simulation/')

x <- foreach(sample_size=SampleSizes, .packages = pckgs, .verbose = TRUE) %:%
  foreach(sampling=Sampling, .packages = pckgs, .verbose = TRUE) %:%
  foreach(rep=Reps, .packages = pckgs, .verbose = TRUE) %dopar% {
    
    counts <- as.data.frame(readRDS(paste0(simulation_path, 'efs', efs, '/rep', rep, '/counts_n', sample_size, '.rds')))
    metadata <- as.data.frame(readRDS(paste0(simulation_path, 'efs', efs, '/rep', rep, '/metadata_n', sample_size, '.rds')))
    counts_transposed <- as.data.frame(t(counts))
    
    metadata <- metadata %>% 
      mutate(
        intervention = ifelse(delta.t >= 9 & delta.t < 20, '1', '0')
      ) %>%
      filter(delta.t %in% sampling)
    
    subsampled_counts <-  counts %>% dplyr::select(metadata$sample_name) 
    
    #transform with clr
    pseudocount <- 0.5
    filtered_subsampled_counts_clr <- apply(subsampled_counts, 2, function(column) {compositions::clr(column + pseudocount)})
    filtered_subsampled_counts_clr <- as.data.frame(filtered_subsampled_counts_clr)
    
    results_df <- data.frame(feature = character(), pvalue = numeric(), stringsAsFactors = FALSE)
    
    for (g in rownames(filtered_subsampled_counts_clr)) {
      
      df <- filtered_subsampled_counts_clr %>% as.data.frame()
      df$genus <- rownames(df)
      df <- df %>% tidyr::pivot_longer(-genus)
      
      colnames(df) <- c('genus', 'sample_name', 'value')
      
      df <- df %>% left_join(metadata, by = 'sample_name') %>% filter(genus == g)
      
      df <- df %>%
        group_by(host_subject_id, intervention) %>%
        summarize(mean = mean(value))
      
      df <- tidyr::pivot_wider(df, id_cols = host_subject_id, names_from = intervention, values_from = mean)
      pval <- wilcox.test(df$`1`, df$`0`, paired = TRUE, exact = FALSE)
      pvalue <- pval$p.value
      
      results_df <- rbind(results_df, data.frame(feature = g, pvalue = pvalue))
    }
    
    ground.truth <- read.table(paste0(simulation_path, 'efs', efs, '/rep', rep, '/gt_features.tsv'), sep = '\t', header=TRUE)
    wilcoxon_results_df <- merge(results_df, ground.truth, by='feature', all.x=TRUE)
    
    wilcoxon_results_df$effect_size <- efs
    wilcoxon_results_df$sample_size <- sample_size
    wilcoxon_results_df$rep <- rep
    wilcoxon_results_df$sampling_density <- length((sampling))
    
    # Save results
    output_path <- paste0('/g/scb/zeller/karwowsk/zeevi_dataset_v5/results/wilcoxon/one_arm/clr/')
    output.file <- paste0(output_path, 'wilcoxon_', 'sampling_', sample_size, '_efs_', efs, '_rep_', rep, '_d_', length(sampling), '.csv')
    write.csv(wilcoxon_results_df, output.file, quote=FALSE)
  }
    
    


