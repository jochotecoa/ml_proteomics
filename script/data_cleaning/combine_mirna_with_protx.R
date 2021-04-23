if (all(file.exists('data/miRNA/score_feats.rds'), file.exists('data/miRNA/mirna_feats_dcast.rds'))) {
  score_feats = readRDS('data/miRNA/score_feats.rds')
  score_feats_all = readRDS('data/miRNA/score_feats_all.rds')
  mirna_feats_dcast = readRDS('data/miRNA/mirna_feats_dcast.rds')
  mirna_feats_dcast_all = readRDS('data/miRNA/mirna_feats_dcast_all.rds')
  seq_depth_mir = readRDS('data/miRNA/seq_depth_mir.rds')
} else {
  source('script/data_cleaning/mirna_data_cleaning.R')
}

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
  merge.data.frame(score_feats_all, 'uniprotswissprot', all.x = T) %>% 
  merge.data.frame(mirna_feats_dcast, 'uniprot_sample', all.x = T) %>% 
  merge.data.frame(mirna_feats_dcast_all, 'uniprot_sample', all.x = T) %>% 
  merge.data.frame(seq_depth_mir, 'sample_name', all.x = T) %>% 
  column_to_rownames('uniprot_sample')

colnames(mrna_prot_df)[grep('seq_depth', colnames(mrna_prot_df))] = 'seq_depth_mirna'

mirna_abs_cols = colnames(mrna_prot_df) %>% 
  subset(., grepl('mirna', .)) %>% 
  subset(., !grepl('log', .))

mirna_log_cols = colnames(mrna_prot_df) %>% 
  subset(., grepl('mirna', .)) %>% 
  subset(., grepl('log', .))

mirna_score_cols = colnames(mrna_prot_df) %>% 
  subset(., grepl('score', .)) 
  
mrna_prot_df[, mirna_abs_cols][is.na(mrna_prot_df[, mirna_abs_cols])] = 0
mrna_prot_df[, mirna_log_cols][is.na(mrna_prot_df[, mirna_log_cols])] = log2(2e-06)
mrna_prot_df[, mirna_score_cols][is.na(mrna_prot_df[, mirna_score_cols])] = 0


# forceLibrary('pbmcapply')
# pb_2 = progressBar(max = length(unique(mrna_prot_df$uniprotswissprot)))
# 
# for (swissprot_i in unique(mrna_prot_df$uniprotswissprot)) {
#   rows_prot = rownames(mrna_prot_df)[mrna_prot_df$uniprotswissprot == swissprot_i]
#   cols_na = is.na(mrna_prot_df[rows_prot,]) %>% colSums() %>% as.logical() %>% colnames(mrna_prot_df)[.]
#   for (col_na in cols_na) {
#     if (class(mrna_prot_df[, col_na]) != 'numeric') {next}
#     rows_na = rows_prot[is.na(mrna_prot_df[rows_prot, col_na])]
#     mrna_prot_df[rows_na, col_na] = median(mrna_prot_df[rows_prot, col_na], na.rm = T)
#   }
#   setTxtProgressBar(pb_2, grep(swissprot_i, unique(mrna_prot_df$uniprotswissprot)))
# }
# close(pb_2)

mrna_prot_df = mrna_prot_df %>%  
  dplyr::select(-c(uniprotswissprot, sample_name))
