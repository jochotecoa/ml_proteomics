if (all(file.exists('data/circRNA/score_feats.rds'), file.exists('data/circRNA/circ_feats_dcast.rds'))) {
  score_feats = readRDS('data/circRNA/score_feats.rds')
  circ_feats_dcast = readRDS('data/circRNA/circ_feats_dcast.rds')
  seq_depth_trx = readRDS('data/seq_depth_trx.rds')
} else {
  source('script/data_cleaning/circrna_data_cleaning.R')
}

seq_depth_trx$sample_name = seq_depth_trx$sample_name %>% 
  gsub(pattern = 'Con', replacement = '')

mrna_prot_df$uniprotswissprot = mrna_prot_df %>% 
  rownames() %>% 
  strsplit('--') %>% 
  sapply('[[', 1) %>% 
  as.character()
mrna_prot_df$sample_name = mrna_prot_df %>% 
  rownames() %>% 
  strsplit('--') %>% 
  sapply('[[', 2) %>% 
  as.character()

i_df = mrna_prot_df
mrna_prot_df = i_df

mrna_prot_df = mrna_prot_df %>% 
  rownames_to_column('uniprot_sample') %>% 
  merge.data.frame(score_feats, 'uniprotswissprot', all.x = T) %>% 
  merge.data.frame(circ_feats_dcast, 'uniprot_sample', all.x = T) %>% 
  merge.data.frame(seq_depth_trx, 'sample_name', all.x = T) %>% 
  column_to_rownames('uniprot_sample')

circrna_abs_cols = colnames(mrna_prot_df) %>% 
  subset(., grepl('circ', .)) %>% 
  subset(., !grepl('log', .))

circrna_log_cols = colnames(mrna_prot_df) %>% 
  subset(., grepl('circ', .)) %>% 
  subset(., grepl('log', .))

circrna_score_cols = colnames(mrna_prot_df) %>% 
  subset(., grepl('score', .)) 

mrna_prot_df[, circrna_abs_cols][is.na(mrna_prot_df[, circrna_abs_cols])] = 0
mrna_prot_df[, circrna_log_cols][is.na(mrna_prot_df[, circrna_log_cols])] = log2(2e-06)
mrna_prot_df[, circrna_score_cols][is.na(mrna_prot_df[, circrna_score_cols])] = 0

mrna_prot_df = mrna_prot_df %>%  
  dplyr::select(-c(uniprotswissprot, sample_name))
