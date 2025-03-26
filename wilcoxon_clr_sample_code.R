.libPaths(c("/g/scb/zeller/karwowsk/R/4.3.3", .libPaths()))
set.seed(10)

library(dplyr)
library(data.table)
library(tidyr)
library(compositions)

efs=3; rep=1; sample_size=10; sampling <- c(1, 5, 8, 12, 15, 19)
simulation_path <- paste0('/g/scb/zeller/karwowsk/zeevi_dataset/simulation/')

counts <- as.data.frame(readRDS(paste0(simulation_path, 'efs', efs, '/rep', rep, '/counts_n', sample_size, '.rds')))
metadata <- as.data.frame(readRDS(paste0(simulation_path, 'efs', efs, '/rep', rep, '/metadata_n', sample_size, '.rds')))

metadata <- metadata %>%
  mutate(
    intervention = ifelse(delta.t >= 9 & delta.t < 20, '1', '0')
  ) %>%
  filter(delta.t %in% sampling)

subsampled_counts <-  counts %>% dplyr::select(metadata$sample_name)

#transform with clr
pseudocount <- 1
filtered_subsampled_counts_clr <- apply(subsampled_counts, 2, function(column) {compositions::clr(column + pseudocount)})
filtered_subsampled_counts_clr <- as.data.frame(filtered_subsampled_counts_clr)

results_df <- data.frame(feature = character(), pvalue = numeric(), stringsAsFactors = FALSE)

g <- rownames(filtered_subsampled_counts_clr)[1]

df <- filtered_subsampled_counts_clr %>% as.data.frame()
df$genus <- rownames(df)
df <- df %>% tidyr::pivot_longer(-genus)
colnames(df) <- c('genus', 'sample_name', 'value')

df <- df %>% left_join(metadata, by = 'sample_name') %>% filter(genus == g)

df <- df %>%
  group_by(host_subject_id, intervention) %>%
  summarize(mean = mean(value))

df <- tidyr::pivot_wider(df, id_cols = host_subject_id, names_from = intervention, values_from = mean)
pval <- wilcox.test(df$`1`, df$`0`, paired = TRUE, exact = FALSE)
pvalue <- pval$p.value

