

circ_df = mrna_df %>% 
  rownames_to_column() %>% 
  dplyr::filter(!grepl('ENST', rowname)) %>% 
  dplyr::rename(circBase_ID = rowname)

cir_mir = readRDS('../ciRmiR_unique_highsponging.rds')
cir_mir2 = readRDS('../ciRmiR_unique_nseedfrequency.rds')

cir_mir[, 'miRBase_ID'] = cir_mir[, 'Var1'] %>% 
  as.character() %>% 
  strsplit(' ') %>% 
  sapply('[[', 1) %>% 
  as.character()
cir_mir2['miRBase_ID'] = cir_mir2[, 'Var1'] %>% 
  as.character() %>% 
  strsplit(' ') %>% 
  sapply('[[', 1) %>% 
  as.character()

cir_mir['circBase_ID'] = cir_mir[, 'Var1'] %>% 
  as.character() %>% 
  strsplit(' ') %>% 
  sapply('[[', 2) %>% 
  as.character()
cir_mir2['circBase_ID'] = cir_mir2[, 'Var1'] %>% 
  as.character() %>% 
  strsplit(' ') %>% 
  sapply('[[', 2) %>% 
  as.character()

colnames(cir_mir)[2] = 'circ_score'
colnames(cir_mir2)[2] = 'circ_score'

circ_df['circBase_ID'] = circ_df[, 'rowname'] %>% 
  cleanCircNames()

circ_df_stringent = circ_df %>% 
  merge.data.frame(cir_mir, 'circBase_ID')
circ_df_all = circ_df %>% 
  merge.data.frame(cir_mir2, 'circBase_ID')

saveRDS(circ_df_stringent, '../circ_df_stringent.rds')
saveRDS(circ_df_all, '../circ_df_all.rds')


circ_df_stringent = circ_df_stringent %>% 
  merge.data.frame(swiss_mir, 'miRBase_ID') %>% 
  merge.data.frame(circ_df, 'circBase_ID') 
circ_df_all = circ_df_all %>% 
  merge.data.frame(swiss_mir_all, 'miRBase_ID') %>% 
  merge.data.frame(circ_df, 'circBase_ID') 


circ_feats = circ_df_stringent %>% 
  dplyr::select(!c(ensembl_gene_id, miRBase_ID, refseq_mrna, Var1, circ_score, circBase_ID)) %>% 
  addVarsProt(fnc_list = c('mean', 'median', 'min', 'max', 'sum', 'sd'), by_str = 'uniprotswissprot')





circ_df_stringent = circ_df_stringent %>% 
  dplyr::select(!c(Var1, circ_score))




circ_scores_stringent = cir_mir %>% 
  merge.data.frame(swiss_mir, 'miRBase_ID')
circ_scores_all = cir_mir2 %>% 
  merge.data.frame(swiss_mir_all, 'miRBase_ID')
