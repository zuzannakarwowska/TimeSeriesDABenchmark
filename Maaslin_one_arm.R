.libPaths(c("/g/scb/zeller/karwowsk/R/4.4.1-gfbf-2023b/", .libPaths()))

library(dplyr)
library(data.table)
library(Maaslin2)
library(foreach)
library(doParallel)


cores <- 3
cl <- makeCluster(cores) 
registerDoParallel(cl)
clusterCall(cl, function(x) .libPaths(x), .libPaths())
pckgs <- c('dplyr', 'data.table', 'Maaslin2', 'foreach', 'doParallel') 


rename_rows <-function(counts){
  cols <- rownames(counts)
  new_cols <- list()
  for (i in 1:length(cols)) {
    new_cols[[i]] <- paste0('feature', as.character(i))
  }
  cols_dict <- cols
  names(cols_dict) <- new_cols
  rename_counts <- counts
  rownames(rename_counts) <- new_cols
  
  return(list('counts' = rename_counts, 'dict' = cols_dict))
}


#dataset='zeevi'
EffectSizes <- c(2)
Reps <- 1:10
Sampling <- list(c(1, 19), c(1, 6, 15, 19), c(1, 5, 8, 12, 15, 19))
SampleSizes <- seq(10, 60, 10)
Datasets <- c('zeevi')

x <- foreach(efs=EffectSizes, .packages = pckgs, .verbose = TRUE) %:%
  foreach(sample_size=SampleSizes, .packages = pckgs, .verbose = TRUE) %:%
  foreach(dataset=Datasets, .packages = pckgs, .verbose = TRUE) %:%
  foreach(sampling=Sampling, .packages = pckgs, .verbose = TRUE) %:%
  foreach(rep=Reps, .packages = pckgs, .verbose = TRUE) %dopar% {
    
    simulation_path <- paste0('/g/scb/zeller/karwowsk/zeevi_dataset_v5/simulation/')
    output_path <- paste0('/g/scb/zeller/karwowsk/zeevi_dataset_v5//results/maaslin/one_arm/log/')
    
    ground.truth <- read.table(paste0(simulation_path, 'efs', efs, '/rep', rep, '/gt_features.tsv'), sep = '\t', header=TRUE)
    counts <- as.data.frame(readRDS(paste0(simulation_path, 'efs', efs, '/rep', rep, '/counts_n', sample_size, '.rds')))
    
    #rename columns as Maaslin changes names
    mapping_obj <- rename_rows(counts)
    counts <- mapping_obj$counts
    cols_dict <- mapping_obj$dict
    
    #prepare metadata
    metadata <- as.data.frame(readRDS(paste0(simulation_path, 'efs', efs, '/rep', rep, '/metadata_n', sample_size, '.rds')))
    
    metadata <- metadata %>% 
      filter(delta.t %in% sampling) %>% 
      mutate(
        intervention = ifelse(delta.t >= 9 & delta.t <= 20, '1', '0')
      )
    
    counts <- counts %>% dplyr::select(metadata$sample_name)
    maaslin_counts <- as.data.frame(t(counts)) #transpose counts, sample_name, features{...}
    maaslin_counts <- data.frame(sample_name = c(metadata$sample_name), maaslin_counts)
    rownames(maaslin_counts) <- c()
    rownames(maaslin_counts) <- maaslin_counts$sample_name
    
    metadata <- as.data.frame(metadata)
    rownames(metadata) <- metadata$sample_name
    
    maaslin_results <- Maaslin2(input_data = maaslin_counts,
                                input_metadata = metadata,
                                output ='None',
                                fixed_effects= c('intervention'),
                                min_prevalence = 0.0,
                                random_effects=c("host_subject_id"), 
                                plot_heatmap = FALSE,
                                plot_scatter = FALSE,
                                normalization="TSS", #TSS
                                transform = 'LOG') #LOG
    
    results_df <- maaslin_results$results %>% dplyr::select('feature', 'pval')
    
    #for samples that were not calculated add pvalue = 1
    feature <- rownames(counts)[!(rownames(counts) %in% results_df$feature)]
    
    if (length(feature) > 0) {
      # If not empty, proceed to create 'missing_features_df' and update it
      missing_features_df <- as.data.frame(feature)
      missing_features_df$pval <- 1  
      
      names(missing_features_df)[1] <- "feature"
      
      Maaslin_results_df <- rbind(results_df, missing_features_df)
    } else {
      Maaslin_results_df <- results_df
    }
    
    Maaslin_results_df$effect_size <- efs
    Maaslin_results_df$iteration <- rep
    Maaslin_results_df$sample_size <- sample_size
    Maaslin_results_df$sampling_density <- length((sampling))
    
    Maaslin_results_df <- Maaslin_results_df %>% 
      dplyr::mutate(feature_assigned = dplyr::recode(feature, !!!cols_dict, .default = NA_character_))
    
    colnames(ground.truth) <- c('feature_assigned', 'bin', 'upregulated')
    Maaslin_results_df <- merge(Maaslin_results_df, ground.truth, 
                                by='feature_assigned', all.x=TRUE)
    names(Maaslin_results_df)[names(Maaslin_results_df) == 'pval'] <- 'pvalue'
    
    output.file <- paste0(output_path, 'maaslin_', 'sampling_', sample_size, '_efs_', efs, '_rep_', rep, '_d_', length(sampling), '.csv')
    write.csv(Maaslin_results_df, output.file, quote=FALSE)
    
  }

#efs=3; rep=1;sampling = c(1, 6, 15, 19); sample_size=10
