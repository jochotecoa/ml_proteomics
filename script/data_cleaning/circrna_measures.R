source('script/functions/functions_JOA.R')
forceLibrary(c('dplyr', 'tibble', 'biomaRt', 'reshape2', 'caret', 'caTools'))

# circ_scores_all = cir_mir2 

circ_df_stringent = readRDS('../ml_bigdata/circ_df_stringent.rds')

circ_feats_str = circ_df_stringent %>% 
  dplyr::select(!c(circBase_ID)) %>% 
  addVarsProt(fnc_list = c('mean', 'median', 'min', 'max', 'sum', 'sd'), by_str = 'uniprotswissprot')

# circ_feats_all = circ_df_all %>% 
#   dplyr::select(!c(circBase_ID, rowname, Var1, circ_score, miRBase_ID, ensembl_gene_id, refseq_mrna)) %>% 
#   addVarsProt(fnc_list = c('mean', 'median', 'min', 'max', 'sum', 'sd'), by_str = 'uniprotswissprot')

circ_feats_melt = circ_feats_str %>% 
  melt()

circ_feats_melt = circ_feats_melt[circ_feats_melt$uniprotswissprot != '', , F]
circ_feats_melt$measure = NA
circ_feats_melt$measure[grepl('mean', circ_feats_melt$variable)] = 'mean'
circ_feats_melt$measure[grepl('median', circ_feats_melt$variable)] = 'median'
circ_feats_melt$measure[grepl('min', circ_feats_melt$variable)] = 'min'
circ_feats_melt$measure[grepl('max', circ_feats_melt$variable)] = 'max'
circ_feats_melt$measure[grepl('sum', circ_feats_melt$variable)] = 'sum'
circ_feats_melt$measure[grepl('sd', circ_feats_melt$variable)] = 'sd'
circ_feats_melt$measure = circ_feats_melt$measure %>% 
  paste0('circ_', .)

circ_feats_melt$sample = circ_feats_melt$variable %>% 
  gsub(pattern = '_mean|_median|_min|_max|_sum|_sd', replacement = '')

circ_feats_melt$uniprot_sample = paste0(circ_feats_melt$uniprotswissprot, '--', circ_feats_melt$sample)

circ_feats_dcast = circ_feats_melt %>% 
  dplyr::select(-c(uniprotswissprot, variable, sample)) %>% 
  dcast(uniprot_sample ~ measure) 

circ_feats_dcast_log2 = circ_feats_dcast[, -grep('uniprot_sample', colnames(circ_feats_dcast))] %>% 
  log2()
colnames(circ_feats_dcast_log2) = colnames(circ_feats_dcast_log2) %>% 
  paste0('_log2')
circ_feats_dcast = circ_feats_dcast %>% 
  cbind.data.frame(circ_feats_dcast_log2)
circ_feats_dcast[circ_feats_dcast == -Inf] = log2(2e-06)

dir.create('data/circRNA')

saveRDS(circ_feats_dcast, 'data/circRNA/circ_feats_dcast.rds')


# Scores ------------------------------------------------------------------

cir_mir = readRDS('../ml_bigdata/circ_prot_scores_stringent.rds')  

score_per_prot = cir_mir %>% 
  dplyr::select(uniprotswissprot, circ_score)

score_feats = score_per_prot %>% 
  addVarsProt(fnc_list = c('mean', 'median', 'min', 'max', 'sum', 'sd'), by_str = 'uniprotswissprot')

score_feats = score_feats[score_feats$uniprotswissprot != '', , F]

saveRDS(score_feats, 'data/circRNA/score_feats.rds')
