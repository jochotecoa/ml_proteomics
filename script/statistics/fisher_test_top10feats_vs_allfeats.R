
source('script/functions/functions_JOA.R')
forceLibrary(c('mlbench', 'caret', 'doParallel', 'dplyr', 'RANN'))

feature_names = readRDS(file = 'data/whole_data/preds_na_omit_cardiac.rds') %>% 
  colnames()

feature_names_table = feature_names %>% 
  strsplit('_') %>% 
  sapply('[[', 1) %>% 
  table %>% 
  data.frame(all_features = .)


top10_features = c('circ_sum', 'circ_mean', 'circ_mean_log2', 'linear_density', 
                   'TPM_value_mean_log2', 'mirna_mean_stringent', 
                   'mirna_mean_stringent_log2', 'mirna_sum_stringent', 
                   'Aa_M_prop', 'Aa_D_prop')

top10_features_table = top10_features %>% 
  strsplit('_') %>% 
  sapply('[[', 1) %>% 
  table %>% 
  data.frame(top10 = .)

table_fisher = merge(x = top10_features_table, y = feature_names_table, 
                     by.x= "top10..", by.y = 'all_features..', all = T) %>% 
  naToZero() %>% 
  column_to_rownames('top10..')

rna_expr = colSums(table_fisher[grepl(pattern = 'circ|mirna|TPM', x = rownames(table_fisher)), ])
non_rna_expr = colSums(table_fisher[!grepl(pattern = 'circ|mirna|TPM', x = rownames(table_fisher)), ])

table_fisher_rnaexpr = rbind(rna_expr, non_rna_expr)
rownames(table_fisher_rnaexpr) = c('rna_expr', 'non_rna_expr')

fisher.test(table_fisher_rnaexpr)

# Fisher's Exact Test for Count Data
# 
# data:  table_fisher_rnaexpr
# p-value = 0.002736
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   1.749321 59.744689
# sample estimates:
# odds ratio
#   8.829389
