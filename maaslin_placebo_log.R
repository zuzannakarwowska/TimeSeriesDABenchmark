.libPaths(c("/g/scb/zeller/karwowsk/R/4.4.1-gfbf-2023b/", .libPaths()))
set.seed(10)
options(repos = c(CRAN = "https://cloud.r-project.org/"))


library(dplyr)
library(data.table)
library(Maaslin2)
library(foreach)
library(doParallel)

cores <- 1
cl <- makeCluster(cores) 
registerDoParallel(cl)
clusterCall(cl, function(x) .libPaths(x), .libPaths())
pckgs <- c('dplyr', 'data.table', 'Maaslin2', 'foreach', 'doParallel') 

#sample_size=10; efs=3;rep=1; sampling = c(10, 19)

rename_rows <- function(counts) {
  new_cols <- paste0("feature", seq_len(nrow(counts)))
  cols_dict <- setNames(rownames(counts), new_cols)
  rownames(counts) <- new_cols
  return(list(counts = counts, dict = cols_dict))
}

Reps <- 1:10
Sampling <- list(c(10, 19), c(10, 12, 16, 18), c(10, 12, 14, 16, 18, 19))
SampleSizes <- seq(10, 60, 10)
efs <-3

x <- foreach(sample_size=SampleSizes, .packages = pckgs, .verbose = TRUE) %:%
  foreach(sampling=Sampling, .packages = pckgs, .verbose = FALSE) %:%
  foreach(rep=Reps, .packages = pckgs, .verbose = FALSE) %dopar% {
    
    simulation_path <- paste0('/g/scb/zeller/karwowsk/zeevi_dataset_v5/simulation/')
    output_path <- paste0('/g/scb/zeller/karwowsk/zeevi_dataset_v5//results/maaslin/placebo/relab/')
    
    ground.truth <- read.table(paste0(simulation_path, 'efs', efs, '/rep', rep, '/gt_features.tsv'), sep = '\t', header=TRUE)
    
    intervention.counts <- as.data.frame(readRDS(paste0(simulation_path, 'efs', efs, '/rep', rep, '/counts_n', sample_size, '.rds')))
    intervention.metadata <- as.data.frame(readRDS(paste0(simulation_path, 'efs', efs, '/rep', rep, '/metadata_n', sample_size, '.rds')))
    
    placebo.counts <- as.data.frame(readRDS(paste0(simulation_path, 'efs', efs, '/rep', rep, '/placebo_counts_n', sample_size, '.rds')))
    placebo.metadata <- as.data.frame(readRDS(paste0(simulation_path, 'efs', efs, '/rep', rep, '/placebo_metadata_n', sample_size, '.rds')))
    
    # Prepare master metadata
    placebo.metadata$Group <- 'placebo'
    intervention.metadata$Group <- 'intervention'
    
    maaslin_metadata <- rbind(placebo.metadata, intervention.metadata)%>% 
      filter(delta.t %in% sampling) %>% 
      mutate(intervention = ifelse(delta.t >= 9 & delta.t <= 20, 1, 0)) 
    
    # add interaction variable
    maaslin_metadata$interaction_term <-  with(maaslin_metadata, interaction(intervention, Group))
    
    maaslin_metadata <- maaslin_metadata %>% relocate('sample_name', .before='host_subject_id')
    colnames(maaslin_metadata) <-  c('ID', "host_subject_id","delta.t", "Group", "intervention", 'interaction_term')
    
    # Prepare master counts
    master_counts <- cbind(placebo.counts, intervention.counts)
    master_counts <- master_counts %>% dplyr::select(maaslin_metadata$ID)
    
    # Rename columns
    mapping_obj <- rename_rows(master_counts)
    master_counts <- mapping_obj$counts
    cols_dict <- mapping_obj$dict
    
    # transpose counts table to have features as columns and sampels as rows
    maaslin_counts <- as.data.frame(t(master_counts))
    maaslin_counts <- data.frame(sample_name = c(maaslin_metadata$ID), maaslin_counts)
    colnames(maaslin_counts) <- c('ID', rownames(master_counts))
    maaslin_counts$ID <- c()
    
    rownames(maaslin_metadata) <- maaslin_metadata$ID
    maaslin_metadata$ID <- c()
    #maaslin_metadata <- maaslin_metadata %>% select(host_subject_id, Group)
    
    #RUN MASSLIN2
    maaslin_results <- Maaslin2(input_data = maaslin_counts,
                                input_metadata = maaslin_metadata,
                                output ="NONE",
                                fixed_effects= c('Group'),
                                min_prevalence = 0,
                                random_effects=c("host_subject_id"), 
                                normalization="TSS",
                                transform = 'LOG',
                                plot_heatmap = FALSE,
                                plot_scatter = FALSE)
    
    results_df <- maaslin_results$results
    
    # For some features Maaslin does not count - fill with pvalue =1 
    feature <- setdiff(rownames(master_counts), results_df$feature)
    missing_features_df <- if (length(feature) > 0) {
      data.frame(feature = feature, pval = 1, stringsAsFactors = FALSE)
    } else {
      NULL
    }
    
    Maaslin_results_df <- results_df %>%
      dplyr::filter(value == "placebo") %>%
      dplyr::select(feature, pval) %>%
      dplyr::bind_rows(missing_features_df) %>%
      dplyr::mutate(
        effect_size = efs,
        iteration = rep,
        sample_size = sample_size,
        sampling_density = length(sampling)
      )
    
    Maaslin_results_df <- Maaslin_results_df %>% 
      dplyr::mutate(
        feature_assigned = dplyr::recode(feature, !!!cols_dict, .default = NA_character_))
    
    Maaslin_results_df$feature <- Maaslin_results_df$feature_assigned
    Maaslin_results_df$feature_assigned <-c()
    Maaslin_results_df <- merge(Maaslin_results_df, ground.truth, by='feature', all.x=TRUE)
    
    output.file <- paste0(output_path, 'maaslin_', 'sampling_', sample_size, '_efs_', efs, '_rep_', rep, '_d_', length(sampling), '.csv')
    write.csv(Maaslin_results_df, output.file, quote=FALSE)
    
  }
