.libPaths(c("/g/scb/zeller/karwowsk/R/4.4.1-gfbf-2023b/", .libPaths()))
set.seed(10)

library(metagenomeSeq)
library(dplyr)
library(foreach)
library(doParallel)
library(purrr)
library(glue)
##Data should be preprocessed and prepared in tab-delimited files. 
## Measurements are stored in a matrix with samples along the columns and features along the rows.


efs=3; rep=1;sampling = c(10,19); sample_size=10

simulation_path <- paste0('/g/scb/zeller/karwowsk/zeevi_dataset_v5/simulation/')
output_path <- paste0('/g/scb/zeller/karwowsk/zeevi_dataset_v5//results/trendyspliner/placebo/')

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
subsampled_counts <-  master_counts %>% dplyr::select(master_metadata$sample_name) 

# filter out bacteria with mean == 0
counts_mean_df <- as.data.frame(rowSums(subsampled_counts))
keep_bacteria <-  counts_mean_df %>% filter(counts_mean_df[,1] > 0) %>% rownames()
counts_filtered <- subsampled_counts %>% 
  filter(row.names(subsampled_counts) %in% keep_bacteria) 

synth_meta_metasplines <- master_metadata
rownames(synth_meta_metasplines) <- synth_meta_metasplines$sample_name
phenotypeData <- AnnotatedDataFrame(synth_meta_metasplines)
obj <- newMRexperiment(counts=counts_filtered, phenoData=phenotypeData)

#Normalize
p = cumNormStatFast(obj)
obj = cumNorm(obj, p = p)

# Run Metasplines
classes <- unique(rownames(obj@assayData$counts))
timeSeriesFits = lapply(classes,function(i){
  fitTimeSeries(obj=obj,
                feature=i,
                class="intervention",
                id="host_subject_id",
                time="delta.t", 
                norm=FALSE,
                log = TRUE, 
                B=999) # B is the number of permutations 
})

names(timeSeriesFits) = classes

df <- lapply(timeSeriesFits,function(i){i[[1]]})[-grep("No statistically",timeSeriesFits)] 
results_df <- as.data.frame(do.call(rbind, df))
results_df$feature <- names(df)
names(results_df)[names(results_df) == 'p.value'] <- 'pvalue'
results_df <- results_df %>% dplyr::select(c(feature, pvalue))

missing_features <- setdiff(rownames(counts), results_df$feature)
missing_features_df <- purrr::map2_dfr(missing_features, rep(1, length(missing_features)), ~ tibble(feature = .x, pvalue = .y))

Metasplines_results_df <- rbind(results_df, missing_features_df)

Metasplines_results_df$effect_size <- efs
Metasplines_results_df$iteration <- rep
Metasplines_results_df$sample_size <- sample_size
Metasplines_results_df$sampling_density <- length((sampling))

Metasplines_results_df <- merge(Metasplines_results_df, ground.truth, by='feature', all.x=TRUE)

output.file <- paste0(output_path, 'metasplines_', 'sampling_', sample_size, '_efs_', efs, '_rep_', rep, '_d_', length(sampling), '.csv')
write.csv(Metasplines_results_df, output.file, quote=FALSE)

