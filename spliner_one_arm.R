.libPaths(c("/g/scb/zeller/karwowsk/R/4.4.1-gfbf-2023b/", .libPaths()))

set.seed(10)

library(splinectomeR)
library(dplyr)
library(foreach)
library(doParallel)
library(compositions)
library(glue)


rename_rows <- function(counts) {
  new_cols <- paste0('feature', seq_along(rownames(counts)))
  cols_dict <- setNames(rownames(counts), new_cols)
  rownames(counts) <- new_cols
  
  return(list(counts = counts, dict = cols_dict))
}

cores <- 2
cl <- makeCluster(cores) 
registerDoParallel(cl)
clusterCall(cl, function(x) .libPaths(x), .libPaths())
pckgs <- c('dplyr', 'data.table', 'compositions', 'foreach', 'doParallel', 'splinectomeR', 'glue') 

#EffectSizes <- c(15, 2, 3)
Reps <- 1:3
SampleSizes <- c(seq(10, 60, 10))
Sampling <- list(c(1, 19), c(1, 6, 15, 19), c(1, 5, 8, 12, 15, 19))
efs <- 3
sampling <- c(1, 6, 15, 19)

x <- foreach(sample_size=SampleSizes, .packages = pckgs, .verbose = TRUE) %:%
  #foreach(sampling=Sampling, .packages = pckgs, .verbose = TRUE) %:%
  foreach(rep=Reps, .packages = pckgs, .verbose = TRUE) %dopar% {
    
    #efs=3; rep=1;sampling = c(1, 6, 15, 19); sample_size=10
    simulation_path <- paste0('/g/scb/zeller/karwowsk/zeevi_dataset_v5/simulation/')
    output_path <- paste0('/g/scb/zeller/karwowsk/zeevi_dataset_v5/results/trendyspliner/one_arm/')
    
    gt_file <- glue("{simulation_path}efs{efs}/rep{rep}/gt_features.tsv")
    counts_file <- glue("{simulation_path}efs{efs}/rep{rep}/counts_n{sample_size}.rds")
    metadata_file <- glue("{simulation_path}efs{efs}/rep{rep}/metadata_n{sample_size}.rds")
    
    ground.truth <- read.table(gt_file, sep = '\t', header = TRUE)
    counts <- readRDS(counts_file) %>% as.data.frame()
    metadata <- readRDS(metadata_file) %>% as.data.frame()
    
    # rename counts name:feature {1...N}
    
    mapping_obj <- rename_rows(counts)
    counts <- mapping_obj$counts
    cols_dict <- mapping_obj$dict
    
    metadata <- metadata %>% 
      mutate(
        intervention = ifelse(delta.t >= 9 & delta.t < 20, '1', '0')
      ) %>%
      filter(delta.t %in% sampling)
    
    # Subsample
    subsampled_counts <- counts %>%
      dplyr::select(all_of(metadata$sample_name)) + 1
    
    # Transform with CLR
    subsampled_counts_clr <- apply(subsampled_counts, 2, clr) %>% t()
    
    subsampled_counts_transposed <- as.data.frame(subsampled_counts_clr) %>%
      dplyr::mutate(sample_name = colnames(subsampled_counts))
    
    # Merge with metadata
    merged_df <- dplyr::left_join(subsampled_counts_transposed, metadata, by = 'sample_name')
    
    # Run splinectomeR
    feature_names <- rownames(subsampled_counts)
    RESULTS = list()
    for (i in seq(1:length(feature_names))){
      trendyspliner_input <- merged_df %>% 
        dplyr::select(paste0('feature', i), 'delta.t', 'host_subject_id') %>%
        mutate(across(c(paste0('feature', i), 'delta.t'), as.numeric))
      
      result_one_group <- splinectomeR::trendyspliner(data = trendyspliner_input,
                                                      xvar = 'delta.t',
                                                      yvar = paste0('feature', i),
                                                      cases = 'host_subject_id',
                                                      mean_center = T,
                                                      ints = length(sampling)-1,
                                                      quiet = T,
                                                      perms = 999)
      res <- list("feature" = paste0('feature', i), 
                  "pvalue" = result_one_group$pval)
      
      RESULTS[[i]] <- as.data.frame(res)
      
    }
    
    results_df <- bind_rows(RESULTS) %>%
      mutate(
        effect_size = efs,
        sample_size = sample_size,
        rep = rep,
        sampling_density = length(sampling),
        feature_assigned = recode(feature, !!!cols_dict, .default = NA_character_) 
      ) %>%
      mutate(feature = feature_assigned) %>%
      select(-feature_assigned) %>%
      left_join(ground.truth, by = "feature")
    
    
    # Save files
    output.file <- glue("{output_path}spliner_sampling_{sample_size}_efs_{efs}_rep_{rep}_d_{length(sampling)}.csv")
    write.csv(results_df, output.file, quote=FALSE)
  
}


