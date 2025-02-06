.libPaths(c("/g/scb/zeller/karwowsk/R/4.4.1-gfbf-2023b/", .libPaths()))

library(metagenomeSeq)
library(dplyr)
library(foreach)
library(doParallel)
library(purrr)

##Data should be preprocessed and prepared in tab-delimited files. 
## Measurements are stored in a matrix with samples along the columns and features along the rows.

cores <- 1
cl <- makeCluster(cores) 
registerDoParallel(cl)
clusterCall(cl, function(x) .libPaths(x), .libPaths())
pckgs <- c('dplyr', 'metagenomeSeq', 'purrr','foreach', 'doParallel') 


EffectSizes <- c(3)
Reps <- 1:10
SampleSizes <- c(seq(10, 60, 10))
Sampling <- list(c(1, 6, 15, 19), c(1, 5, 8, 12, 15, 19), c(1, 5, 6, 8, 12, 14, 16, 19))
efs <-3 

x <- foreach(sample_size=SampleSizes, .packages = pckgs, .verbose = TRUE) %:%
  foreach(sampling=Sampling, .packages = pckgs, .verbose = TRUE) %:%
  foreach(rep=Reps, .packages = pckgs, .verbose = TRUE) %dopar% {
    
    simulation_path <- paste0('/g/scb/zeller/karwowsk/zeevi_dataset_v5/simulation/')
    output_path <- paste0('/g/scb/zeller/karwowsk/zeevi_dataset_v5/results/metasplines/one_arm/')
    
    gt_file <- glue("{simulation_path}efs{efs}/rep{rep}/gt_features.tsv")
    counts_file <- glue("{simulation_path}efs{efs}/rep{rep}/counts_n{sample_size}.rds")
    metadata_file <- glue("{simulation_path}efs{efs}/rep{rep}/metadata_n{sample_size}.rds")
    
    ground.truth <- read.table(gt_file, sep = '\t', header = TRUE)
    counts <- readRDS(counts_file) %>% as.data.frame()
    metadata <- readRDS(metadata_file) %>% as.data.frame()
      
      metadata <- metadata %>% 
        mutate(intervention = ifelse(delta.t >= 9 & delta.t < 20, 1, 0)) %>%
        filter(delta.t %in% sampling) 
      
      subsampled_counts <-  counts %>% dplyr::select(metadata$sample_name) 
      
      # filter out bacteria with mean == 0
      counts_mean_df <- as.data.frame(rowSums(subsampled_counts))
      keep_bacteria <-  counts_mean_df %>% filter(counts_mean_df[,1] > 0) %>% rownames()
      counts_filtered <- subsampled_counts %>% filter(row.names(subsampled_counts) %in% keep_bacteria)
      
      # create input object
      synth_meta_metasplines <- metadata
      rownames(synth_meta_metasplines) <- synth_meta_metasplines$sample_name
      phenotypeData <- AnnotatedDataFrame(synth_meta_metasplines)
      obj <- newMRexperiment(counts=counts_filtered, phenoData=phenotypeData)
      
      #Normalize
      p = cumNormStatFast(obj)
      obj = cumNorm(obj, p = p)
      
      set.seed(1000)
      
      # Run Metasplines
      classes <- unique(rownames(obj@assayData$counts))
      timeSeriesFits = lapply(classes,function(i){
        fitTimeSeries(obj=obj,
                      feature=i,
                      class="intervention",
                      id="host_subject_id",
                      time="delta.t",
                      B=1) # B is the number of permutations 
      })
      
      names(timeSeriesFits) = classes
      
      df <- lapply(timeSeriesFits,function(i){i[[1]]})[-grep("No statistically",timeSeriesFits)] 
      results_df <- as.data.frame(do.call(rbind, df))
      results_df$feature <- names(df)
      names(results_df)[names(results_df) == 'p.value'] <- 'pvalue'
      results_df <- results_df %>% dplyr::select(c(feature, pvalue))
      
      missing_features <- setdiff(rownames(counts), results_df$feature)
      missing_features_df <- purrr::map2_dfr(missing_features, rep(1, length(missing_features)), ~ tibble(feature = .x, pvalue = .y))
      
      Metasplines_results_df <- bind_rows(results_df, missing_features_df) %>%
        mutate(
          effect_size      = efs,
          iteration        = rep,
          sample_size      = sample_size,
          sampling_density = length(sampling)
        ) %>%
        left_join(ground.truth, by = "feature")
      
      output.file <- paste0(output_path, 'metasplines_', 'sampling_', sample_size, '_efs_', efs, '_rep_', rep, '_d_', length(sampling), '.csv')
      write.csv(Metasplines_results_df, output.file, quote=FALSE)
    }
    
  }

