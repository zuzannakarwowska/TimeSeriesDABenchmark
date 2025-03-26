library(dplyr)
library(pROC)
library(glue)

compute_auc <- function(df) {
  controls <- df %>% filter(upregulated == 0) %>% pull(pvalue)
  cases <- df %>% filter(upregulated != 0) %>% pull(pvalue)
  
  if (length(controls) == 0 || length(cases) == 0) {
    return(NaN)  # Return NaN if there's no data for controls or cases
  }
  
  auc_result <- tryCatch({
    roc_obj <- roc(controls = controls, cases = cases, direction = '>')
    auc_value <- roc_obj$auc
    if (auc_value < 0 | auc_value > 1) NaN else auc_value
  }, error = function(e) NaN)  # If an error occurs, return NaN
  
  return(auc_result)
}

model_name <- 'baseline'
MEASURES <- data.frame()
transformation <- 'log'

for (N in c(10, 20, 40)) {
  for (efs in c(125, 15, 2, 3, 5)){
    for (rep in seq(10)) { 
      
      file <- glue('/Users/zkarwowska/Desktop/EMBL_project/zeevi_dataset_v5/results/baseline/{transformation}/baseline_sampling_{N}_efs_{efs}_rep_{rep}.csv')
      
      
      if(file.exists(file)) {
        
        df <- read.csv(file)
        names(df)[names(df) == 'pval'] <- 'pvalue'
        df$q.value <- p.adjust(df$pvalue, method = "BH")
        
        df$ytrue <- abs(df$upregulated)
        
        df <- df %>%
          mutate(ypred = ifelse((q.value <= 0.05), 1, 0))
        
        AUC <- roc(controls = df %>% filter(upregulated == 0) %>% pull(pvalue),
                   cases = df %>% filter(upregulated != 0) %>% pull(pvalue),
                   direction = '>')$auc
        
        FDR <- df %>%
          filter(pvalue < 0.05) %>%
          mutate(false_positive = ifelse((pvalue < 0.05 & upregulated == 0), TRUE, FALSE)) %>%
          summarise(fdr = sum(false_positive) / n())
        
        precision <- sum(df$ypred & df$ytrue) / sum(df$ypred)
        recall <- sum(df$ypred & df$ytrue) / sum(df$ytrue)
        
        sum_ypred <- sum(df$ypred)
        sum_ytrue <- sum(df$ytrue)
        
        # Create a row as a dataframe for the current iteration
        current_measure <- data.frame(auc = AUC,
                                      FDR = FDR$fdr,
                                      precision = precision,
                                      recall = recall,
                                      ytrue = sum_ytrue,
                                      ypred = sum_ypred,
                                      effect_size = efs,
                                      sample_size = N,
                                      sampling=sampling,
                                      model = model_name,
                                      transformation = transformation,
                                      version = 'v5',
                                      iteration = rep)
        
        # Bind the row to the main MEASURES dataframe
        MEASURES <- rbind(MEASURES, current_measure)
      } else {
        # File does not exist, so skip to the next file
        next
      }
    }
  }
}

out_dir <- '/Users/zkarwowska/new_embl_folder/zeevi_dataset_v5/results/results_files/one_arm/'
if (transformation != 'none'){
  write.csv(MEASURES, glue('{out_dir}/{model_name}_{transformation}_one_arm.csv'))
} else {
  write.csv(MEASURES, glue('{out_dir}/{model_name}_one_arm.csv'))
}
#write.csv(MEASURES, glue('{out_dir}/{model_name}_one_arm.csv'))

