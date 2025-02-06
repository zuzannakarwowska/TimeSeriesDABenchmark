.libPaths(c("/g/scb/zeller/karwowsk/R/4.4.1-gfbf-2023b/", .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org/"))

set.seed(10)

library(dplyr)
library(data.table)
library(NBZIMM)
library(foreach)
library(doParallel)
library(glue)

cores <- 1
cl <- makeCluster(cores) 
registerDoParallel(cl)
clusterCall(cl, function(x) .libPaths(x), .libPaths())
pckgs <- c('dplyr', 'data.table', 'NBZIMM', 'foreach', 'doParallel', 'glue') 


prepare_zigmm_df <-function(counts){
  
  offset <- c(colSums(counts))
  
  relative_abundance_df <- scale(counts, center = FALSE, 
                                 scale = colSums(counts))
  
  rownames(relative_abundance_df) <- rownames(counts)
  relative_abundance_transposed <- as.data.frame(t(relative_abundance_df))
  relative_abundance_transposed <- tibble::rownames_to_column(relative_abundance_transposed, "sample_name")
  return(relative_abundance_transposed)
}

Reps <- 1:10
SampleSizes <- c(seq(10, 60, 20))
Sampling <- list(c(1, 19), c(1, 6, 15, 19), c(1, 5, 8, 12, 15, 19))
efs <- 3

x <- foreach(sample_size=SampleSizes, .packages = pckgs, .verbose = TRUE) %:%
  foreach(sampling=Sampling, .packages = pckgs, .verbose = TRUE) %:%
  foreach(rep=Reps, .packages = pckgs, .verbose = TRUE) %dopar% {
    
    simulation_path <- paste0('/g/scb/zeller/karwowsk/zeevi_dataset_v5/simulation/')
    output_path <- paste0('/g/scb/zeller/karwowsk/zeevi_dataset_v5/results/zigmm/one_arm/')
    
    gt_file <- glue("{simulation_path}efs{efs}/rep{rep}/gt_features.tsv")
    counts_file <- glue("{simulation_path}efs{efs}/rep{rep}/counts_n{sample_size}.rds")
    metadata_file <- glue("{simulation_path}efs{efs}/rep{rep}/metadata_n{sample_size}.rds")
    
    ground.truth <- read.table(gt_file, sep = '\t', header = TRUE)
    counts <- readRDS(counts_file) %>% as.data.frame()
    metadata <- readRDS(metadata_file) %>% as.data.frame()
    
    nbmm_counts <- prepare_zigmm_df(counts)
    
    metadata <- metadata %>% 
      mutate(
        intervention = ifelse(delta.t >= 9 & delta.t < 20, '1', '0')
      ) %>%
      filter(delta.t %in% sampling)
    
    
    #first column must be sample_name
    metadata <-  metadata %>% dplyr::select(sample_name, everything()) 
    
    nbmm.input.df <- left_join(metadata, nbmm_counts, by = 'sample_name') %>%
      relocate(sample_name, .before = host_subject_id)  
    
    nbmm_metadata <- nbmm.input.df %>% dplyr::select(sample_name:intervention)
    nbmm_counts <- nbmm.input.df %>% dplyr::select(-(sample_name:intervention))
    
    NBMM.RES <- mms(min.p = 0, 
                    y = nbmm_counts, 
                    fixed = ~  intervention,
                    random = ~ 1 | host_subject_id,
                    data = nbmm_metadata,
                    zi_fixed = ~ host_subject_id,
                    method = "zig")
    
    
    feature_names <- colnames(nbmm_counts)
    
    fit_df <- NBMM.RES$fit # Moved outside the loop to avoid repeated access
    # Use lapply to iterate over features and store results in a list
    feature_res_list <- lapply(feature_names, function(feature) {
      feature_res <- summary(fit_df[[feature]])
      feature_res_table <- feature_res['tTable']
      feature_res_table_df <- as.data.frame(feature_res_table[[1]])
      feature_res_table_df$feature <- feature
      
      if (nrow(feature_res_table_df) == 1) {
        feature_res_table_df$predictor <- 'error'
        feature_res_table_df$`p-value` <- 1
      } else {
        feature_res_table_df$predictor <- c('Intercept', 'intervention')
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
    
    missing_features <- rownames(counts)[!(rownames(counts) %in% results_df$feature)]
    
    # Create a data frame with these missing features and assign a default p-value of 1
    missing_features_df <- data.frame(
      feature = missing_features,
      pval = rep(1, length(missing_features))
    )
    
    nbmm_results_df <- bind_rows(results_df, missing_features_df) %>%
      mutate(
        effect_size = efs,
        iteration = rep,
        sample_size = sample_size,
        sampling_density = length(sampling)
      ) %>%
      left_join(ground.truth, by = "feature")  
    
    output.file <- paste0(output_path, 'zigmm_', 'sampling_', sample_size, '_efs_', efs, '_rep_', rep, '_d_', length(sampling), '.csv')
    write.csv(nbmm_results_df, output.file, quote=FALSE)
    
  }

stopCluster(cl)

