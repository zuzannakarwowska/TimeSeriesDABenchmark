.libPaths(c("/g/scb/zeller/karwowsk/R/4.4.1-gfbf-2023b/", .libPaths()))
set.seed(10)
options(repos = c(CRAN = "https://cloud.r-project.org/"))

library(dplyr)
library(data.table)
library(NBZIMM)
library(foreach)
library(doParallel)

cores <- 1
cl <- makeCluster(cores) 
registerDoParallel(cl)
clusterCall(cl, function(x) .libPaths(x), .libPaths())
pckgs <- c('dplyr', 'data.table', 'NBZIMM', 'foreach', 'doParallel') 


prepare_zigmm_df <-function(counts){
  
  #prepare input for ZIGMM where columns are [sample_name, features{...}, offset]
  #in relative abundance
  
  relative_abundance_df <- scale(counts, center = FALSE, 
                                 scale = colSums(counts))
  
  rownames(relative_abundance_df) <- rownames(counts)
  relative_abundance_transposed <- as.data.frame(t(relative_abundance_df))
  relative_abundance_transposed <- tibble::rownames_to_column(relative_abundance_transposed, "sample_name")
  return(relative_abundance_transposed)
}

#dataset = 'moore'; sample_size=5; efs=15;rep=1; sampling = c(1, 19)
Reps <- 1:10
Sampling <- list(c(10, 19), c(10, 12, 16, 18), c(10, 12, 14, 16, 18, 19))
SampleSizes <- seq(10, 60, 10)
efs <-3

x <- foreach(sample_size=SampleSizes, .packages = pckgs, .verbose = TRUE) %:%
  foreach(sampling=Sampling, .packages = pckgs, .verbose = FALSE) %:%
  foreach(rep=Reps, .packages = pckgs, .verbose = FALSE) %dopar% {
    
    simulation_path <- paste0('/g/scb/zeller/karwowsk/zeevi_dataset_v5/simulation/')
    output_path <- paste0('/g/scb/zeller/karwowsk/zeevi_dataset_v5//results/zigmm/placebo/')

    intervention.counts <- as.data.frame(readRDS(paste0(simulation_path, 'efs', efs, '/rep', rep, '/counts_n', sample_size, '.rds')))
    intervention.metadata <- as.data.frame(readRDS(paste0(simulation_path, 'efs', efs, '/rep', rep, '/metadata_n', sample_size, '.rds')))
    intervention.metadata$Group <- 'intervention'
    
    placebo.counts <- as.data.frame(readRDS(paste0(simulation_path, 'efs', efs, '/rep', rep, '/placebo_counts_n', sample_size, '.rds')))
    placebo.metadata <- as.data.frame(readRDS(paste0(simulation_path, 'efs', efs, '/rep', rep, '/placebo_metadata_n', sample_size, '.rds')))
    placebo.metadata$Group <- 'placebo'
    
    ground.truth <- read.table(paste0(simulation_path, 'efs', efs, '/rep', rep, '/gt_features.tsv'), sep = '\t', header=TRUE)
    
    # Prepare metadata
    metadata <- rbind(placebo.metadata, intervention.metadata)
    metadata <- metadata %>% 
      filter(delta.t %in% sampling) %>% 
      mutate(intervention = ifelse(delta.t >= 9 & delta.t <= 20, 1, 0)) 
    
    metadata <-  metadata %>% dplyr::select(sample_name, everything()) 
    
    # Prepare counts
    master_counts <- cbind(placebo.counts, intervention.counts)
    master_counts <- master_counts %>% dplyr::select(metadata$sample_name)
    
    nbmm_counts <- prepare_zigmm_df(master_counts)
    
    nbmm.input.df <- left_join(metadata, nbmm_counts, by = 'sample_name') %>%
      relocate(sample_name, .before = host_subject_id)  
    
    
    nbmm_metadata <- nbmm.input.df %>% dplyr::select(sample_name:intervention)
    nbmm_counts <- nbmm.input.df %>% dplyr::select(-(sample_name:intervention))
    
    # Run model
    NBMM.RES <- mms(min.p = 0, 
                    y = nbmm_counts, 
                    fixed = ~  Group,
                    random = ~ 1 | host_subject_id,
                    data = nbmm_metadata,
                    zi_fixed = ~ host_subject_id,
                    method = "zig")
    
    feature_names <- colnames(nbmm_counts)
    
    fit_df <- NBMM.RES$fit 
    
    feature_res_list <- lapply(feature_names, function(feature) {
      
      feature_res <- summary(fit_df[[feature]])
      feature_res_table <- feature_res['tTable']
      feature_res_table_df <- as.data.frame(feature_res_table[[1]])
      
      if (nrow(feature_res_table_df) == 1) {
        feature_res_table_df$predictor <- 'error'
        feature_res_table_df$`p-value` <- 1
        feature_res_table_df$feature <- feature
      } else {
        feature_res_table_df <- feature_res_table_df %>% filter(row.names(feature_res_table_df) == 'Groupplacebo')
        feature_res_table_df$feature <- feature
        feature_res_table_df$predictor <- 'intervention'
      }
      
      return(feature_res_table_df)
    })
    
    # Combine all the data frames in the list into one
    NBZIMM.PREDICTION.DF <- bind_rows(feature_res_list)
    
    # Filter, select, and rename as before
    results_df <- NBZIMM.PREDICTION.DF %>%
      filter(predictor == 'intervention') %>%
      dplyr::select(feature, predictor, `t-value`, `p-value`) %>%
      rename(pval = `p-value`) %>% 
      dplyr::select('feature', 'pval')
    
    missing_features <- rownames(master_counts)[!(rownames(master_counts) %in% results_df$feature)]
    
    if (length(missing_features) > 0) {
      # If not empty, proceed to create 'missing_features_df' and update it
      missing_features_df <- as.data.frame(missing_features)
      missing_features_df$pval <- 1  # Assign a default p-value of 1 to missing features
      
      # Assuming 'feature' column is needed to match 'results_df' structure
      names(missing_features_df)[1] <- "feature"
      
      nbmm_results_df <- rbind(results_df, missing_features_df)
    } else {
      nbmm_results_df <- results_df
    }
    
    
    row.names(nbmm_results_df) <- NULL
    nbmm_results_df$effect_size <- efs
    nbmm_results_df$iteration <- rep
    nbmm_results_df$sample_size <- sample_size
    nbmm_results_df$sampling_density <- length((sampling))
    
    #merge with gt
    nbmm_results_df <- merge(nbmm_results_df, ground.truth, by='feature', all.x=TRUE)
    
    output.file <- paste0(output_path, 'zigmm_', 'sampling_', sample_size, '_efs_', efs, '_rep_', rep, '_d_', length(sampling), '.csv')
    write.csv(nbmm_results_df, output.file, quote=FALSE)
    
  }

