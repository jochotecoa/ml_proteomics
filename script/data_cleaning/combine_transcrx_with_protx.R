# Combine mRNA with proteins ----------------------------------------------

mrna_prot_df = mrna_unip_df %>% 
  merge.data.frame(y = prot_df, by = 'uniprot_sample') # all.y = T if there's something else aside TPM

mrna_prot_df %>% dplyr::select(where(is.numeric)) %>% na.omit %>% cor

mrna_prot_df = mrna_prot_df %>% 
  dplyr::select(!matches('\\.y'))

colnames(mrna_prot_df) = colnames(mrna_prot_df) %>% 
  gsub(pattern = '\\.x', replacement = '')



mrna_prot_df = mrna_prot_df %>%
  dplyr::filter(!is.na(proteomics_value))

# mrna_prot_df = mrna_prot_df[mrna_prot_df$TPM_value != 0, ]

# mrna_prot_df$log2_TPM_value = log2(mrna_prot_df$TPM_value)

mrna_prot_df = mrna_prot_df %>% 
  dplyr::select(!c(uniprotswissprot, sample_name)) %>% 
  remove_rownames() %>% 
  column_to_rownames('uniprot_sample')

mrna_prot_df$strand_median = mrna_prot_df$strand_median %>% as.factor()
mrna_prot_df$strand_min = mrna_prot_df$strand_min %>% as.factor()
mrna_prot_df$strand_max = mrna_prot_df$strand_max %>% as.factor()

# saveRDS(object = mrna_prot_df, file = paste0('data/whole_raw_dataset', tissue, '.rds'))

