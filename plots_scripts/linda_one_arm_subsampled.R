.libPaths(c("/g/scb/zeller/karwowsk/R/4.2.2-foss-2022b", .libPaths()))

library(dplyr)
library(data.table)
library(LinDA)
library(foreach)
library(doParallel)

cores <- 4
cl <- makeCluster(cores) 
registerDoParallel(cl)
clusterCall(cl, function(x) .libPaths(x), .libPaths())
pckgs <- c('dplyr', 'data.table', 'LinDA', 'foreach', 'doParallel') 


EffectSizes <- c(3)#15, 2, 3, 5)
Reps <- 1:10
Datasets <- c('zeevi')
Sampling <- list(c(1, 19), c(1, 6, 15, 19), c(1, 5, 8, 12, 15, 19), seq(1:19))
SampleSizes <- seq(10, 60, 10)

x <- foreach(efs=EffectSizes, .packages = pckgs, .verbose = TRUE) %:%
  foreach(sample_size=SampleSizes, .packages = pckgs, .verbose = TRUE) %:%
  foreach(dataset=Datasets, .packages = pckgs, .verbose = TRUE) %:%
  foreach(sampling=Sampling, .packages = pckgs, .verbose = TRUE) %:%
  foreach(rep=Reps, .packages = pckgs, .verbose = TRUE) %dopar% {
    
    simulation_path <- paste0('/g/scb/zeller/karwowsk/zeevi_dataset_v5/simulation/')
    output_path <- paste0('/g/scb/zeller/karwowsk/zeevi_dataset_v5/results/linda/one_arm/')
    
    counts <- as.data.frame(readRDS(paste0(simulation_path, 'efs', efs, '/rep', rep, '/counts_n', sample_size, '.rds')))
    metadata <- as.data.frame(readRDS(paste0(simulation_path, 'efs', efs, '/rep', rep, '/metadata_n', sample_size, '.rds')))
    counts_transposed <- as.data.frame(t(counts))
    
    metadata <- metadata %>% 
      mutate(
        intervention = ifelse(delta.t >= 9 & delta.t < 20, '1', '0')
      ) %>%
      filter(delta.t %in% sampling)

    subsampled_counts <-  counts %>% dplyr::select(metadata$sample_name) #bacteria as rows
    
    linda.obj <- linda(subsampled_counts,
                       metadata,
                       prev.cut = 0,lib.cut=0, corr.cut= 0,
                       formula = '~ intervention + (1|host_subject_id)'
    )
    
    LinDA_results <- as.data.frame(linda.obj$output)
    LinDA_results$feature <- rownames(LinDA_results)
    LinDA_results$effect_size <- efs
    LinDA_results$iteration <- rep
    LinDA_results$sample_size <- sample_size
    LinDA_results$sampling_density <- length((sampling))
    
    
    ground.truth <- read.table(paste0(simulation_path, 'efs', efs, '/rep', rep, '/gt_features.tsv'), sep = '\t', header=TRUE)
    
    LinDA_results_df <- merge(LinDA_results, ground.truth, by='feature', all.x=TRUE)
    LinDA_results_df <- LinDA_results_df %>% dplyr::rename('pvalue'='intervention1.pvalue')
    
    output.file <- paste0(output_path, 'linda_', 'sampling_', sample_size, '_efs_', efs, '_rep_', rep, '_d_', length(sampling), '.csv')
    write.csv(LinDA_results_df, output.file, quote=FALSE)
  }

    
    
    



