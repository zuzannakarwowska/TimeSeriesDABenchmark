.libPaths(c("/g/scb/zeller/karwowsk/R/4.4.1-gfbf-2023b/", .libPaths()))
set.seed(10)

options(repos = c(CRAN = "https://cloud.r-project.org/"))
dyn.load("/g/scb/zeller/karwowsk/R/4.4.1-gfbf-2023b/nloptr/libs/nloptr.so")

library(dplyr)
library(data.table)
library(LinDA)
library(foreach)
library(doParallel)

Sys.setenv(LD_LIBRARY_PATH="/g/scb/zeller/karwowsk/R/4.4.1-gfbf-2023b/nloptr/libs")
library(nloptr)

cores <- 2
cl <- makeCluster(cores) 
registerDoParallel(cl)
clusterCall(cl, function(x) .libPaths(x), .libPaths())
pckgs <- c('dplyr', 'data.table', 'LinDA', 'foreach', 'doParallel') 

efs <- 3
Reps <- 1:10
Sampling <- list(c(10, 19), c(10, 12, 16, 18), c(10, 12, 14, 16, 18, 19))
SampleSizes <- seq(10, 60, 10)


x <- foreach(sample_size=SampleSizes, .packages = pckgs, .verbose = TRUE) %:%
  foreach(sampling=Sampling, .packages = pckgs, .verbose = TRUE) %:%
  foreach(rep=Reps, .packages = pckgs, .verbose = TRUE) %dopar% {
    
    simulation_path <- paste0('/g/scb/zeller/karwowsk/zeevi_dataset_v5/simulation/')
    output_path <- paste0('/g/scb/zeller/karwowsk/zeevi_dataset_v5//results/linda/placebo/')
    
    ground.truth <- read.table(paste0(simulation_path, 'efs', efs, '/rep', rep, '/gt_features.tsv'), sep = '\t', header=TRUE)
    
    intervention.counts <- as.data.frame(readRDS(paste0(simulation_path, 'efs', efs, '/rep', rep, '/counts_n', sample_size, '.rds')))
    intervention.metadata <- as.data.frame(readRDS(paste0(simulation_path, 'efs', efs, '/rep', rep, '/metadata_n', sample_size, '.rds')))
    
    placebo.counts <- as.data.frame(readRDS(paste0(simulation_path, 'efs', efs, '/rep', rep, '/placebo_counts_n', sample_size, '.rds')))
    placebo.metadata <- as.data.frame(readRDS(paste0(simulation_path, 'efs', efs, '/rep', rep, '/placebo_metadata_n', sample_size, '.rds')))
    
    #add Group variable to metadata
    intervention.metadata$Group <- 'intervention'
    placebo.metadata$Group <- 'placebo'
    
    master_metadata <- rbind(placebo.metadata, intervention.metadata)
    master_metadata <- master_metadata %>% 
      filter(delta.t %in% sampling) %>% 
      mutate(intervention = ifelse(delta.t >= 9 & delta.t <= 20, 1, 0)) 
    
    # Prepare master counts
    master_counts <- cbind(placebo.counts, intervention.counts)
    master_counts <- master_counts %>% dplyr::select(master_metadata$sample_name)
    
    
    linda.obj <- linda(master_counts,
                       master_metadata, 
                       prev.cut = 0, 
                       formula = '~ intervention:Group + (1|host_subject_id)')
    
    LinDA_results <- as.data.frame(linda.obj$output)
    LinDA_results$feature <- rownames(LinDA_results)
    
    results_df <- LinDA_results %>%
      dplyr::select(feature, intervention.Groupintervention.pvalue) %>%
      dplyr::rename(pvalue = intervention.Groupintervention.pvalue) %>%
      dplyr::left_join(ground.truth, by = "feature") %>%
      dplyr::mutate(
        effect_size_group = efs,
        iteration = rep,
        sample_size = sample_size,
        sampling_density = length(sampling)
      )
    
    output.file <- paste0(output_path, 'linda_', 'sampling_', sample_size, '_efs_', efs, '_rep_', rep, '_d_', length(sampling), '.csv')
    write.csv(results_df, output.file, quote=FALSE)
    
  }

