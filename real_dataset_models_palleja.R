.libPaths(c("/g/scb/zeller/karwowsk/R/4.4.1-gfbf-2023b/", .libPaths()))

library(dplyr)
library(glue)
library(LinDA)
library(Maaslin2)
library(tibble)
library(ggplot2)
library(tidyr)

# ---- pvalue adjust function ----

adjust_pvalues <- function(pvalues){
  
  # Adjust only the non-NaN p-values
  non_nan_index <- !is.nan(pvalues)
  adjusted_pvalues <- rep(NA, length(pvalues))
  adjusted_pvalues[non_nan_index] <- p.adjust(pvalues[non_nan_index], method = 'BY')
  q.value <- round(adjusted_pvalues, 4)
  
}


setwd('/g/scb/zeller/karwowsk/zeevi_dataset_v5/real_datasets/')


run_linda <- function(counts, metadata){
  
  # counts - bacteria as rows, sample as columns
  linda.obj <- linda(counts,
                     metadata,
                     prev.cut = 0,
                     lib.cut = 0,
                     corr.cut = 0,
                     formula = '~ Timepoint + (1|Individual_ID)'
  )
  
  LinDA_results <- as.data.frame(linda.obj$output)
  LinDA_results$feature <- rownames(LinDA_results)
  LinDA_results$q.value <- adjust_pvalues(LinDA_results$Timepoint.pvalue)
  LinDA_results <- LinDA_results %>%
    mutate(ypred = ifelse((q.value <= 0.05), 1, 0))
  
  linda_predictions <- LinDA_results %>% 
    filter(ypred == 1) %>% 
    rownames() 
  return(linda_predictions)
}


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



## input Maaslin dataframe with sample_ID as rownames and bacteria as columns
## metadata has sample_ID as rownames
run_maaslin <- function(maaslin_counts, maaslin_metadata){
  # maaslin_counts: sample id in columns, bacteria as rows
  mapping_obj <- rename_rows(maaslin_counts)
  maaslin_counts <- mapping_obj$counts
  cols_dict <- mapping_obj$dict
  
  maaslin_counts <- data.frame(sample_name = c(maaslin_metadata$Sample_ID), t(maaslin_counts))
  maaslin_counts <- maaslin_counts %>% dplyr::select(-sample_name)
  rownames(maaslin_metadata) <- maaslin_metadata$Sample_ID
  
  maaslin_results <- Maaslin2(input_data = maaslin_counts,
                              input_metadata = maaslin_metadata,
                              output ='None',
                              fixed_effects= c('Timepoint'),
                              min_prevalence = 0.0,
                              random_effects=c("Individual_ID"), 
                              plot_heatmap = FALSE,
                              plot_scatter = FALSE,
                              normalization="TSS", #TSS
                              transform = 'LOG') #LOG
  
  maaslin_results_df <- maaslin_results$results %>% dplyr::select('feature', 'pval')
  maaslin_results_df <- maaslin_results_df %>% 
    dplyr::mutate(feature_assigned = dplyr::recode(feature, !!!cols_dict, .default = NA_character_))
  
  maaslin_results_df$q.value <- adjust_pvalues(maaslin_results_df$pval)
  maaslin_results_df <- maaslin_results_df %>%
    mutate(ypred = ifelse((q.value <= 0.05), 1, 0))
  
  maaslin_results_df$feature <- maaslin_results_df$feature_assigned
  maaslin_results_df$feature_assigned <- c()
  maaslin_predictions <- maaslin_results_df %>% filter(ypred == 1) %>%  pull(feature)
  return(maaslin_predictions)
  
}


### ZIGMM

library(NBZIMM)

prepare_zigmm_df <-function(counts){
  
  offset <- c(colSums(counts))
  
  relative_abundance_df <- scale(counts, center = FALSE, 
                                 scale = colSums(counts))
  
  rownames(relative_abundance_df) <- rownames(counts)
  relative_abundance_transposed <- as.data.frame(t(relative_abundance_df))
  relative_abundance_transposed <- tibble::rownames_to_column(relative_abundance_transposed, "Sample_ID")
  return(relative_abundance_transposed)
}


run_zigmm <- function(counts, metadata){
  
  zinb_prepared_counts <- prepare_zigmm_df(counts)
  
  nbmm.input.df <- left_join(metadata, zinb_prepared_counts, by = 'Sample_ID') %>%
    relocate(Sample_ID, .before = Individual_ID)  
  
  nbmm_metadata <- nbmm.input.df %>% dplyr::select(Sample_ID:Timepoint)
  nbmm_counts <- nbmm.input.df %>% dplyr::select(-(Sample_ID:Timepoint)) %>% dplyr::select(-c(Library_size))
  
  NBMM.RES <- mms(min.p = 0, 
                  y = nbmm_counts, 
                  fixed = ~  Timepoint,
                  random = ~ 1 | Individual_ID,
                  data = nbmm_metadata,
                  zi_fixed = ~ Individual_ID,
                  method = "zig")
  
  
  feature_names <- colnames(nbmm_counts)
  
  fit_df <- NBMM.RES$fit 
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
  
  nbmm_results_df <- bind_rows(results_df, missing_features_df) 
  nbmm_results_df$q.value <- adjust_pvalues(nbmm_results_df$pval)
  nbmm_results_df <- nbmm_results_df %>%
    mutate(ypred = ifelse((q.value <= 0.05), 1, 0))
  
  zigmm_predictions <- nbmm_results_df %>% 
    filter(ypred == 1) %>% 
    pull(feature)
  
  return(zigmm_predictions)
}


otu_table <- read.csv('ready_palleja_counts.csv', row.names = 'X') # sample id in columns
metadata <- read.csv('ready_palleja_metadata.csv', row.names = 'X')
absolute_counts <- read.csv('ready_palleja_absolute.csv', row.names = 1) #bacteria as columns
rownames(metadata) <- metadata$Sample_ID

results <- data.frame()

for (rep in seq(1, 5)){
  for (sample_size in c(5, 10, 12)) {
    
    # Sample metadata and subset OTU tables
    small_metadata <- metadata %>%
      filter(Individual_ID %in% sample(unique(metadata$Individual_ID), sample_size)) %>%
      arrange(desc(Individual_ID))
    
    small_otu_table <- otu_table %>% dplyr::select(small_metadata$Sample_ID)
    absolute_counts_small <- absolute_counts %>% dplyr::select(small_metadata$Sample_ID)
    
    # Run models
    linda_predictions   <- run_linda(small_otu_table[1:120, ], small_metadata)
    maaslin_predictions <- run_maaslin(small_otu_table, small_metadata)
    zigmm_predictions   <- run_zigmm(absolute_counts_small, small_metadata)
    
    # Define sets
    set1 <- linda_predictions
    set2 <- maaslin_predictions
    set3 <- zigmm_predictions
    
    # Calculate lengths and overlaps
    linda_len         <- length(set1)
    maaslin_len       <- length(set2)
    zigmm_len         <- length(set3)
    linda_maaslin_len <- length(intersect(set1, set2))
    linda_zigmm_len   <- length(intersect(set1, set3))
    maaslin_zigmm_len <- length(intersect(set2, set3))
    all_len           <- length(Reduce(intersect, list(set1, set2, set3)))
    
    # Combine into a row
    row <- data.frame(
      sample_size       = sample_size,
      rep               = rep,
      linda_n           = linda_len,
      maaslin_n         = maaslin_len,
      zigmm_n           = zigmm_len,
      linda_maaslin_n   = linda_maaslin_len,
      linda_zigmm_n     = linda_zigmm_len,
      maaslin_zigmm_n   = maaslin_zigmm_len,
      all_three_overlap = all_len
    )
    
    # Append to results
    results <- rbind(results, row)
  }
}

write.csv(results, 'palleja_results.csv')


## plot Venn 

# Load the library
library(eulerr)

## Define sets
set1 <- linda_predictions
set2 <- maaslin_predictions
set3 <- zigmm_predictions

linda_prediction_len <- length(set1)
maaslin_prediction_len <- length(set2)
zigmm_prediction_len <- length(set3)
linda_maslin_len <- length(intersect(set1, set2))
linda_zigmm_len <- length(intersect(set1, set3))
maaslin_zigmm_len <- length(intersect(set2, set3))
all_len <- length(Reduce(intersect, list(set1, set2, set3)))



# Calculate overlaps
fit <- euler(c(
  "LinDA" = length(set1),
  "MaAsLin2" = length(set2),
  "ZIGMM" = length(set3),
  "LinDA&MaAsLin2" = length(intersect(set1, set2)),
  "LinDA&ZIGMM" = length(intersect(set1, set3)),
  "MaAsLin2&ZIGMM" = length(intersect(set2, set3)),
  "LinDA&MaAsLin2&ZIGMM" = length(Reduce(intersect, list(set1, set2, set3)))
))

# Plot
plot(fit,
     fills = list(fill = c("skyblue", "pink", "lightgreen"), alpha = 0.6),
     edges = list(col = "black"),
     labels = list(font = 2),
     quantities = TRUE,
     main = "Overlap Between LinDA, MaAsLin2, and ZIGMM Predictions")

