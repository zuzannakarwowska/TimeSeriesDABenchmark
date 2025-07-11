.libPaths(c("/g/scb/zeller/karwowsk/R/4.4.1-gfbf-2023b/", .libPaths()))

library(dplyr)
library(glue)
library(LinDA)
library(Maaslin2)
library(tibble)
library(ggplot2)
library(tidyr)

library(NBZIMM)

# ---- pvalue adjust function ----

adjust_pvalues <- function(pvalues){
  
  # Adjust only the non-NaN p-values
  non_nan_index <- !is.nan(pvalues)
  adjusted_pvalues <- rep(NA, length(pvalues))
  adjusted_pvalues[non_nan_index] <- p.adjust(pvalues[non_nan_index], method = 'BY')
  q.value <- round(adjusted_pvalues, 4)
  
}


## -----------------------------------------------------------------------------
# LinDA
## -----------------------------------------------------------------------------

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

## -----------------------------------------------------------------------------
# Maaslin2
## -----------------------------------------------------------------------------


## input Maaslin dataframe with sample_ID as rownames and bacteria as columns
## metadata has sample_ID as rownames

run_maaslin <- function(counts, metadata) {
  maaslin_counts <- as.data.frame(t(counts))  # transpose: samples as rows
  maaslin_results <- Maaslin2(
    input_data = maaslin_counts,
    input_metadata = metadata,
    output = 'None',
    fixed_effects = "Timepoint",
    random_effects = "Individual_ID",
    min_prevalence = 0.0,
    normalization = "CLR",
    transform = "NONE",
    plot_heatmap = FALSE,
    plot_scatter = FALSE
  )
  
  predictions <- maaslin_results$results %>%
    dplyr::select(feature, pval) %>%
    mutate(
      q.value = adjust_pvalues(pval),
      ypred = as.integer(q.value <= 0.05)
    ) %>%
    filter(ypred == 1) %>%
    pull(feature)
  
  return(predictions)
}


## -----------------------------------------------------------------------------
# ZIGMM
## -----------------------------------------------------------------------------


run_zigmm <- function(otu_table, metadata){
  
  counts <- as.data.frame(t(otu_table)) #bacteria as columns
  
  NBMM.RES <- mms(
    min.p = 0,
    y = log(counts + 1),
    fixed = ~ Timepoint,
    random = ~ 1 | Individual_ID,
    data = metadata,
    zi_fixed = ~ Individual_ID,
    method = "zig"
  )
  
  pvals <- sapply(names(counts), function(f) {
    fit <- NBMM.RES$fit[[f]]
    if (length(fit) > 0) summary(fit)$tTable["Timepoint", "p-value"] else 1
  })
  
  results <- tibble(
    feature = names(pvals),
    pval = pvals,
    q.value = adjust_pvalues(pval)
  ) %>%
    mutate(ypred = as.integer(q.value <= 0.05))
  
  zigmm_predictions <- filter(results, ypred == 1) %>% pull(feature)
  return(zigmm_predictions)
  
}

## -----------------------------------------------------------------------------
# Run models
## -----------------------------------------------------------------------------

run_subsampled_models <- function(otu_table, metadata, sample_sizes){
  
  results <- data.frame()
  
  for (rep in seq(1, 3)){
    for (sample_size in sample_sizes) {
      
      small_metadata <- metadata %>%
        filter(Individual_ID %in% sample(unique(metadata$Individual_ID), sample_size)) %>%
        arrange(desc(Individual_ID))
      
      small_otu_table <- otu_table %>% dplyr::select(small_metadata$Sample_ID)
      
      # Run models
      linda_predictions   <- run_linda(small_otu_table[1:(nrow(otu_table)-2), ], small_metadata)
      maaslin_predictions <- run_maaslin(small_otu_table, small_metadata)
      zigmm_predictions   <- run_zigmm(small_otu_table, small_metadata)
      
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
  return(results)
}


setwd('/g/scb/zeller/karwowsk/zeevi_dataset_v5/real_datasets/')

palleja_otu_table <- read.csv('ready_palleja_counts.csv', row.names = 'X') # sample id in columns
palleja_metadata <- read.csv('ready_palleja_metadata.csv', row.names = 'X') %>% mutate(Timepoint = scale(Timepoint)) # Scale timepoint
rownames(palleja_metadata) <- palleja_metadata$Sample_ID

palleja_sample_sizes <- c(5, 10, 12)
palleja_results <- run_subsampled_models(palleja_otu_table, palleja_metadata, palleja_sample_sizes)
write.csv(palleja_results, 'palleja_results.csv')

mesnage_otu_table <- read.csv('ready_mesnage_counts.csv', row.names = 'X') # sample id in columns
mesnage_metadata <- read.csv('ready_mesnage_metadata.csv', row.names = 'X') %>% mutate(Timepoint = scale(Timepoint)) # Scale timepoint
rownames(mesnage_metadata) <- mesnage_metadata$Sample_ID

mesnage_sample_sizes <- c(10, 20, 30, 40, 50, 60, 70)
mesnage_results <- run_subsampled_models(mesnage_otu_table, mesnage_metadata, mesnage_sample_sizes)
write.csv(mesnage_results, 'mesnage_results.csv')


venkataraman_otu_table <- read.csv('baxter_otu_table.csv', row.names = 'X') # sample id in columns
venkataraman_metadata <- read.csv('baxter_metadata.csv', row.names = 'X') %>% 
  dplyr::select(c('Sample_ID', 'Individual_ID','Timepoint', 'study')) %>%
  mutate(Timepoint = scale(Timepoint)) %>%
  filter(study == 'venkataraman') %>% dplyr::select(-c(study))
venkataraman_otu_table <- venkataraman_otu_table %>% dplyr::select(venkataraman_metadata$Sample_ID)
rownames(venkataraman_metadata) <- venkataraman_metadata$Sample_ID

venkataraman_sample_sizes <- c(5, 10)
venkataraman_results <- run_subsampled_models(venkataraman_otu_table, venkataraman_metadata, venkataraman_sample_sizes)
write.csv(venkataraman_results, 'venkataraman_results.csv')


#zigmm_predictions   <- run_zigmm(mesnage_otu_table, mesnage_metadata)
#maaslin_predictions <- run_maaslin(mesnage_otu_table, mesnage_metadata)
#linda_predictions <- run_linda(mesnage_otu_table[1:mesnage_n, ], mesnage_metadata)
#



#write.csv(results, 'palleja_results.csv')


## plot Venn 

# Load the library
library(eulerr)

## Define sets
set1 <- mesnage_results$linda_predictions
set2 <- mesnage_results$maaslin_predictions
set3 <- mesnage_results$zigmm_predictions

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
