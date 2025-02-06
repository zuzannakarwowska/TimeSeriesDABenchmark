.libPaths(c("/g/scb/zeller/karwowsk/R/4.4.1-gfbf-2023b/", .libPaths()))

set.seed(10)

library(splinectomeR)
library(dplyr)
library(compositions)
library(foreach)
library(doParallel)


rename_rows <- function(counts) {
  new_cols <- paste0('feature', seq_along(rownames(counts)))
  cols_dict <- setNames(rownames(counts), new_cols)
  rownames(counts) <- new_cols
  
  return(list(counts = counts, dict = cols_dict))
}

sample_size=10; efs=2;rep=1; sampling = c(10, 19)

if (length(sampling)==2){
  ints <- 2
} else{
  ints <- length(sampling)/2-1
}

simulation_path <- paste0('/g/scb/zeller/karwowsk/zeevi_dataset_v5/simulation/')
output_path <- paste0('/g/scb/zeller/karwowsk/zeevi_dataset_v5//results/trendyspliner/placebo/clr/')

ground.truth <- read.table(paste0(simulation_path, 'efs', efs, '/rep', rep, '/gt_features.tsv'), sep = '\t', header=TRUE)

intervention.counts <- as.data.frame(readRDS(paste0(simulation_path, 'efs', efs, '/rep', rep, '/counts_n', sample_size, '.rds')))
intervention.metadata <- as.data.frame(readRDS(paste0(simulation_path, 'efs', efs, '/rep', rep, '/metadata_n', sample_size, '.rds')))

placebo.counts <- as.data.frame(readRDS(paste0(simulation_path, 'efs', efs, '/rep', rep, '/placebo_counts_n', sample_size, '.rds')))
placebo.metadata <- as.data.frame(readRDS(paste0(simulation_path, 'efs', efs, '/rep', rep, '/placebo_metadata_n', sample_size, '.rds')))


# Prepare master metadata
placebo.metadata$Group <- 'placebo'
intervention.metadata$Group <- 'intervention'

master_metadata <- rbind(placebo.metadata, intervention.metadata)
master_metadata <- master_metadata %>% 
  filter(delta.t %in% sampling) %>% 
  mutate(intervention = ifelse(delta.t >= 9 & delta.t <= 20, 1, 0)) 

# Prepare master counts
master_counts <- cbind(placebo.counts, intervention.counts)

# Rename columns
mapping_obj <- rename_rows(master_counts)
master_counts <- mapping_obj$counts
cols_dict <- mapping_obj$dict

# Subsample count table

subsampled_counts <- master_counts %>%
  dplyr::select(all_of(master_metadata$sample_name)) + 1
# Transform with CLR
master_counts_clr_df <- apply(subsampled_counts, 2, clr) %>% t() %>% as.data.frame()
master_counts_clr_df$sample_name <- colnames(subsampled_counts)
merged_df <- merge(master_counts_clr_df, master_metadata, by='sample_name', all.x=TRUE)

RESULTS = list()
for (i in seq(1:length(rownames(intervention.counts)))){

  trendyspliner_input <- merged_df %>% 
    dplyr::select(paste0('feature', i), 'delta.t', 'host_subject_id', 'intervention', 'Group') %>%
    mutate(across(c(paste0('feature', i), 'delta.t'), as.numeric))
  
  result_one_group <- splinectomeR::trendyspliner(data = trendyspliner_input,
                                                  xvar = 'delta.t',
                                                  yvar = paste0('feature', i),
                                                  cases = 'host_subject_id',
                                                  group = 'Group',
                                                  mean_center = F,
                                                  ints = ints,
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
    sampling_density = length(sampling)
  )

results_df <- results_df %>% 
  dplyr::mutate(feature_assigned = dplyr::recode(feature, !!!cols_dict, .default = NA_character_))
results_df$feature <- results_df$feature_assigned
results_df$feature_assigned <- c()

results_df <- merge(results_df, ground.truth, by='feature', all.x=TRUE)

output.file <- paste0(output_path, 'spliner_', 'sampling_', sample_size, '_efs_', efs, '_rep_', rep, '_d_', length(sampling), '.csv')
write.csv(results_df, output.file, quote=FALSE)

