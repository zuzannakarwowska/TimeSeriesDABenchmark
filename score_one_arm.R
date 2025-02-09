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
      
      #file <-  paste0('/g/scb/zeller/karwowsk/zeevi_dataset_v5/results/wilcoxon/one_arm/log/wilcoxon_sampling_', N, '_efs_3_rep_', rep, '_d_4.csv')
      #file <-  paste0('/g/scb/zeller/karwowsk/zeevi_dataset_v5/results/maaslin/one_arm/', transformation, '/maaslin_sampling_', N, '_efs_', efs, '_rep_', rep, '_d_', sampling, '.csv')
      res_dir <- '/Users/zkarwowska/Desktop/EMBL_project/zeevi_dataset_v5/zeevi_dataset_v5/results/'
      file <-  paste0('{res_dir}/', model_name, 'one_arm/log/spliner', '_sampling_', N, '_efs_3_rep_', rep, '_d_', sampling, '.csv')
      
      if(file.exists(file)) {
        
        df <- read.csv(file)
        names(df)[names(df) == 'pval'] <- 'pvalue'
        df[is.na(df)] <- 1
        #df$q.value <- round(p.adjust(df %>% pull(pvalue), method = 'BY'), 4)
        
        non_nan_index <- !is.nan(df$pvalue)
        
        # Adjust only the non-NaN p-values
        adjusted_pvalues <- rep(NA, length(df$pvalue))  # Initialize with NA
        adjusted_pvalues[non_nan_index] <- p.adjust(df$pvalue[non_nan_index], method = 'BH')
        
        # Round the adjusted p-values
        df$q.value <- round(adjusted_pvalues, 4)
        df$ytrue <- abs(df$upregulated)
        
        df <- df %>%
          mutate(ypred = ifelse((q.value < 0.05), 1, 0))
        print(sum(df$ypred))
        
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
                                      model = MODEL,
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

out_dir <- '/Users/zkarwowska/Desktop/EMBL_project/zeevi_dataset_v5/zeevi_dataset_v5/results/results_files/'
write.csv(MEASURES, paste0('{out_dir}/{model_name}_{transformation}_one_arm.csv'))
