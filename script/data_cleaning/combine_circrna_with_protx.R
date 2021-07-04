
if (tissue = 'hepatic') {
  if (all(file.exists('data/circRNA/score_feats.rds'), file.exists('data/circRNA/circ_feats_dcast.rds'))) {
    score_feats = readRDS('data/circRNA/score_feats.rds')
    circ_feats_dcast = readRDS('data/circRNA/circ_feats_dcast.rds')
    seq_depth_trx = readRDS('data/seq_depth_trx.rds')
  } else {
    source('script/data_cleaning/circrna_data_cleaning.R')
  }
  
} else {
  if (all(file.exists(paste0('data/circRNA/score_feats', tissue, '.rds')), file.exists(paste0('data/circRNA/circ_feats_dcast', tissue, '.rds')))) {
    score_feats = readRDS(paste0('data/circRNA/score_feats', tissue, '.rds'))
    circ_feats_dcast = readRDS(paste0('data/circRNA/circ_feats_dcast', tissue, '.rds'))
    seq_depth_trx = readRDS(paste0('data/transcriptomics/seq_depth_trx', tissue, '.rds'))
    } else {
      source('script/data_cleaning/circrna_data_cleaning.R')
    }
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


circ_feats_dcast$sample_name = circ_feats_dcast$uniprot_sample %>% 
  strsplit('--') %>% 
  sapply('[[', 2) %>% 
  as.character()


noprotsamples = unique(mrna_prot_df$sample_name)[!(unique(mrna_prot_df$sample_name) %in% unique(circ_feats_dcast$sample_name))]
if (length(noprotsamples) > 0) {
  warning(paste(paste0(noprotsamples, collapse = ', '), 'did not have proteomics samples, only transcriptomics'))
  mrna_prot_df = mrna_prot_df %>% 
    dplyr::filter(!grepl(pattern = paste0(noprotsamples, collapse = '|'), sample_name))
}

notranscrsamples = unique(circ_feats_dcast$sample_name)[!(unique(circ_feats_dcast$sample_name) %in% unique(mrna_prot_df$sample_name))]
if (length(notranscrsamples) > 0) {
  warning(paste(paste0(notranscrsamples, collapse = ', '), 'did not have transcriptomics samples, only proteomics'))
  circ_feats_dcast = circ_feats_dcast %>% 
    dplyr::filter(!grepl(pattern = paste0(notranscrsamples, collapse = '|'), sample_name))
  
}

all(unique(circ_feats_dcast$sample_name) %in% unique(mrna_prot_df$sample_name)) %>% 
  stopifnot('sample names different between transcrx and protx'= .)



seq_depth_trx$sample_name = seq_depth_trx$sample_name %>% 
  gsub(pattern = 'Con', replacement = '')


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
