#.libPaths(c("/g/scb/zeller/karwowsk/R/4.2.2-foss-2022b", .libPaths()))
#.libPaths(c("/g/scb/zeller/karwowsk/R/4.3.3", .libPaths()))
.libPaths(c("/g/scb/zeller/karwowsk/R/4.4.1-gfbf-2023b", .libPaths()))

#install.packages("/g/scb/zeller/karwowsk/lobster_pckg5/LOBSTER-main/", repos = NULL, type = "source") 

library(LOBSTER)
library(dplyr)
library(readr)
library(foreach)
library(doParallel)

make_metadata <- function(metadata, sample_name){
  
  filtered_metadata <- metadata %>% dplyr::select(sample_id, subject_id) %>% filter(sample_id %in% sample_name)
  
  names(filtered_metadata)[names(filtered_metadata) == 'sample_id'] <- 'sample_name'
  names(filtered_metadata)[names(filtered_metadata) == 'subject_id'] <- 'host_subject_id'
  filtered_metadata$delta.t <- 0
  
  return(filtered_metadata)
}
#generate samples with lobster
generate_data <- function(feature_table, metadata, sample_size, efs){
  
  lobster <- lobster(training_counts = p.counts, 
                     training_meta = p.meta, 
                     initialization_counts = feature_table, 
                     initialization_meta = metadata)
  
  lobster <- train_models(lobster)
  lobster <- gen_mean_matrices(lobster, samp = 30)
  lobster <- gen_synth_meta(lobster)
  
  group_size <- as.integer(length(unique(metadata$host_subject_id)))
  
  lobster <- add_grouping(lobster, 
                          grouping_label = 'Group', 
                          n_in_group = sample_size)
  
  lobster <- select_shift_feat(lobster, 'upregulated', nbins = 10, effect_size=efs)
  
  lobster <- implant_effects(lobster = lobster, 
                             shift_feat = 'upregulated', 
                             grouping = 'Group', 
                             intervention_start = 10, 
                             intervention_end = 20, 
                             effect_size = efs, 
                             is_symmetric = TRUE, 
                             transition_duration_start = 3, 
                             transition_duration_end = 3, 
                             transition_type_start = 'linear',
                             transition_type_end = 'linear')
  
  lobster <- add_noise_and_sparsity(lobster)
  lobster <- summarize_synth_counts(lobster)
  
  
  synth_counts <- as.data.frame(synth_counts(lobster))
  synth_meta <- synth_meta(lobster)
  features <- as.data.frame(lobster@feature_selections$upregulated)
  
  simulated_counts <- lobster@synth_counts
  simulated_metadata <- lobster@synth_meta
  gt_features <- lobster@feature_selections$upregulated
  
  LobsterOutput <- list("simulated_counts" = simulated_counts, 
                        "simulated_metadata" = simulated_metadata, 
                        "gt_features" = gt_features)
  return(LobsterOutput)
  
}

cores <- 5
cl <- makeCluster(cores) 
registerDoParallel(cl)
clusterCall(cl, function(x) .libPaths(x), .libPaths())
pckgs <- c('dplyr', 'LOBSTER', 'readr', 'foreach', 'doParallel') 


EffectSizes <- c(2, 3)
Reps <-  seq(1:10)
efs <-2 
x <- foreach(rep=Reps, .packages = pckgs, .verbose = TRUE) %dopar% {

    counts.df <- as.data.frame(readRDS('/g/scb/zeller/karwowsk/zeevi_dataset_v4/input_data/zeevi_ready_profiles.rds'))
    unassigned_fraction <- tail(counts.df, 1)
    counts.df <- counts.df[!(row.names(counts.df) == 'unassigned'),]
    
    metadata <- as.data.frame(readRDS('/g/scb/zeller/karwowsk/zeevi_dataset_v4/input_data/zeevi_metadata.rds'))
    
    #randomly sample 100 subjects
    
    N <- 200
    sample_name <- sample(colnames(counts.df), N)
    intervention_sample_name <- sample(sample_name, N/2)
    placebo_sample_name <- setdiff(sample_name, intervention_sample_name)
    
    # Create metadata
    intervention_metadata <- make_metadata(metadata, intervention_sample_name)
    placebo_metadata <- make_metadata(metadata, placebo_sample_name)
    
    intervention_feature_table <- counts.df %>% dplyr::select(all_of(intervention_sample_name)) %>% as.matrix()
    placebo_feature_table <- counts.df %>% dplyr::select(all_of(placebo_sample_name)) %>% as.matrix()
    
    # Generate  placebo data
    placebo_generated_data <- generate_data(placebo_feature_table, placebo_metadata, 0, 0)
    
    placebo_metadata <- placebo_generated_data$simulated_metadata
    placebo_counts <- placebo_generated_data$simulated_counts
    placebo_gt <- placebo_generated_data$gt_features
    
    # Generate  intervention data
    intervention_generated_data <- generate_data(intervention_feature_table, intervention_metadata, N/2, efs)
    
    intervention_metadata <- intervention_generated_data$simulated_metadata
    intervention_counts <- intervention_generated_data$simulated_counts
    intervention_gt <- intervention_generated_data$gt_features
    
    if (efs == 1.5){
      output_path <- paste0('/g/scb/zeller/karwowsk/zeevi_dataset_v5/simulation/efs15/rep', rep, '/')
    } else {
      output_path <- paste0('/g/scb/zeller/karwowsk/zeevi_dataset_v5/simulation/efs', efs, '/rep', rep, '/')
    }
    
    write_rds(placebo_counts, paste0(output_path, 'counts_placebo.rds'))
    write_rds(placebo_metadata, paste0(output_path, 'metadata_placebo.rds'))
    write_tsv(placebo_gt, paste0(output_path, 'gt_features_placebo.tsv'))
    
    
    write_rds(intervention_counts, paste0(output_path, 'counts.rds'))
    write_rds(intervention_metadata, paste0(output_path, 'metadata.rds'))
    write_tsv(intervention_gt, paste0(output_path, 'gt_features.tsv'))
}

    
