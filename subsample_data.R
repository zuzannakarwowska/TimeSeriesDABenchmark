#.libPaths(c("/g/scb/zeller/karwowsk/R/4.2.0-foss-2021b", .libPaths()))
.libPaths(c("/g/scb/zeller/karwowsk/R/4.3.3", .libPaths()))

library(dplyr)
library(tibble)
library(foreach)
library(doParallel)


#simulation_path <- '/g/scb/zeller/karwowsk/zeevi_dataset/simulation/efs3/rep1/'

filter_on_prevalence <- function(df){
  min_occurrence <- ncol(df) * 0.1
  occurrence_counts <- rowSums(df > 0)
  selected_features <- names(occurrence_counts[occurrence_counts > min_occurrence])
  return(selected_features)
}

subsample_data <- function(dataset, efs, rep, sample_size){
  
  simulation_path <- paste0('/g/scb/zeller/karwowsk/zeevi_dataset_v5/simulation/')

  counts <- as.data.frame(readRDS(paste0(simulation_path, 'efs', efs, '/rep', rep, '/counts.rds')))
  meatdata <- as.data.frame(readRDS(paste0(simulation_path, 'efs', efs, '/rep', rep, '/metadata.rds')))
  
  placebo_counts <- as.data.frame(readRDS(paste0(simulation_path, 'efs', efs, '/rep', rep, '/counts_placebo.rds')))
  placebo_meatdata <- as.data.frame(readRDS(paste0(simulation_path, 'efs', efs, '/rep', rep, '/metadata_placebo.rds')))
  
  ground.truth <- read.table(paste0(simulation_path, 'efs', efs, '/rep', rep, '/gt_features.tsv'), sep = '\t', header=TRUE)
  keep_bacteria <- ground.truth %>% filter(upregulated != 0) %>% pull(feature)
  
  df <- (rbind(t(counts), t(placebo_counts)))
  filtered_bacteria <- filter_on_prevalence(as.data.frame(t(df)))
  features <- unique(c(filtered_bacteria, keep_bacteria))
  
  filtered_counts <- filter(counts, rownames(counts) %in% features)
  filtered_counts_placebo <- filter(placebo_counts, rownames(placebo_counts) %in% features)
  
  #filter metadata only for before and after intervention, subsample to n subjects, add column with intervention
  subsampled_metadata <- meatdata %>% 
    dplyr::filter(host_subject_id %in% sample(unique(meatdata$host_subject_id), sample_size))
  
  subsampled_metadata_placebo <- placebo_meatdata %>% 
    dplyr::filter(host_subject_id %in% sample(unique(placebo_meatdata$host_subject_id), sample_size))
  
  #subsample and filter counts
  subsampled_counts <-  filtered_counts %>% dplyr::select(subsampled_metadata$sample_name)
  subsampled_counts_placebo <-  filtered_counts_placebo %>% dplyr::select(subsampled_metadata_placebo$sample_name)
  
  
  saveRDS(subsampled_counts, file = paste0(simulation_path, 'efs', efs, '/rep', rep, '/counts_n', sample_size, '.rds'))
  saveRDS(subsampled_metadata, file = paste0(simulation_path, 'efs', efs, '/rep', rep, '/metadata_n', sample_size, '.rds'))
  
  saveRDS(subsampled_counts_placebo, file = paste0(simulation_path, 'efs', efs, '/rep', rep, '/placebo_counts_n', sample_size, '.rds'))
  saveRDS(subsampled_metadata_placebo, file = paste0(simulation_path, 'efs', efs, '/rep', rep, '/placebo_metadata_n', sample_size, '.rds'))
}

#subsample_data('zeevi', 3, 1, 5)


cores <- 3
cl <- makeCluster(cores) 
registerDoParallel(cl)
clusterCall(cl, function(x) .libPaths(x), .libPaths())
pckgs <- c('dplyr', 'data.table', 'foreach', 'doParallel') 


EffectSizes <- c(2)#15, 2, 3, 5)
Reps <- seq(1, 10)
Datasets <- c('zeevi')#,'moore')
SampleSizes <- c(10, 20, 30, 40, 50, 60)

x <- foreach(efs=EffectSizes, .packages = pckgs, .verbose = TRUE) %:%
  foreach(sample_size=SampleSizes, .packages = pckgs, .verbose = TRUE) %:%
  foreach(dataset=Datasets, .packages = pckgs, .verbose = TRUE) %:%
  foreach(rep=Reps, .packages = pckgs, .verbose = TRUE) %dopar% {
    
    subsample_data(dataset, efs, rep, sample_size)
    
  }
