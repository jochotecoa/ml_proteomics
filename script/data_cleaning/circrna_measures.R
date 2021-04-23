source('script/functions/functions_JOA.R')
forceLibrary(c('dplyr', 'tibble', 'biomaRt', 'reshape', 'caret', 'caTools'))

cir_mir = readRDS('../ml_bigdata/circ_prot_scores_stringent.rds')  
# circ_scores_all = cir_mir2 

circ_df_stringent = readRDS('../ml_bigdata/circ_df_stringent.rds')

circ_feats_str = circ_df_stringent %>% 
  dplyr::select(!c(circBase_ID, rowname, Var1, circ_score, miRBase_ID, ensembl_gene_id, refseq_mrna)) %>% 
  addVarsProt(fnc_list = c('mean', 'median', 'min', 'max', 'sum', 'sd'), by_str = 'uniprotswissprot')

circ_feats_all = circ_df_all %>% 
  dplyr::select(!c(circBase_ID, rowname, Var1, circ_score, miRBase_ID, ensembl_gene_id, refseq_mrna)) %>% 
  addVarsProt(fnc_list = c('mean', 'median', 'min', 'max', 'sum', 'sd'), by_str = 'uniprotswissprot')
