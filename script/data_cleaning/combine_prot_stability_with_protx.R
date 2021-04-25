if (file.exists('data/stability/prot_stab_feats.rds')) {
  prot_stab_feats = readRDS('data/stability/prot_stab_feats.rds')
} else {
  source('script/data_cleaning/protein_stability_data_cleaning.R')
}

mrna_prot_df$uniprotswissprot = mrna_prot_df %>% 
  rownames() %>% 
  strsplit('--') %>% 
  sapply('[[', 1) %>% 
  as.character()

i_df = mrna_prot_df
mrna_prot_df = i_df

mrna_prot_df = mrna_prot_df %>% 
  rownames_to_column('uniprot_sample') %>% 
  merge.data.frame(prot_stab_feats, 'uniprotswissprot', all.x = T) %>% 
  column_to_rownames('uniprot_sample')

mrna_prot_df = mrna_prot_df %>%  
  dplyr::select(-c(uniprotswissprot, sample_name))
