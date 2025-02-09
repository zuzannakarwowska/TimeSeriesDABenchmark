library(dplyr)
library(pROC)
library(glue)

model_name <- 'trendyspliner'
MEASURES <- data.frame()
transformation <- 'clr'

efs <-3 
for (N in c(10, 20, 30, 40, 50, 60)) {
  for (rep in seq(10)) { 
    for (sampling in c(2, 4, 6)){
      
      res_dir <- glue('/Users/zkarwowska/Desktop/EMBL_project/zeevi_dataset_v5/results/{model_name}/placebo/')
      
      if (transformation == 'NONE'){
        file <- glue('{res_dir}/{model_name}_sampling_{N}_efs_3_rep_{rep}_d_{sampling}.csv')
      } else {
        #file <- glue('{res_dir}/{transformation}/{model_name}_sampling_{N}_efs_3_rep_{rep}_d_{sampling}.csv')
        file <- glue('{res_dir}/{transformation}/spliner_sampling_{N}_efs_3_rep_{rep}_d_{sampling}.csv')
        
      }
      
      
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
                                      setup = 'placebo',
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

out_dir <- '/Users/zkarwowska/Desktop/EMBL_project/zeevi_dataset_v5/results/results_files/'

if (transformation == 'NONE'){
  write.csv(MEASURES, glue('{out_dir}/{model_name}_placebo.csv'))
} else {
  write.csv(MEASURES, glue('{out_dir}/{model_name}_{transformation}_placebo.csv'))
}
  
  

