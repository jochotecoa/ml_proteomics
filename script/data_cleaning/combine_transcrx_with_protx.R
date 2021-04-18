# Combine mRNA with proteins ----------------------------------------------

mrna_prot_df = mrna_unip_df %>% 
  merge.data.frame(y = prot_df, by = 'uniprot_sample') # all.y = T if there's something else aside TPM

mrna_prot_df %>% dplyr::select(where(is.numeric)) %>% na.omit %>% cor

mrna_prot_df = mrna_prot_df %>% 
  dplyr::select(!contains('.y'))

colnames(mrna_prot_df) = colnames(mrna_prot_df) %>% 
  gsub(pattern = '.x', replacement = '')



mrna_prot_df = mrna_prot_df %>%
  dplyr::filter(!is.na(proteomics_value))

mrna_prot_df = mrna_prot_df[mrna_prot_df$TPM_value != 0, ]

mrna_prot_df$log2_TPM_value = log2(mrna_prot_df$TPM_value)

mrna_prot_df = mrna_prot_df %>% 
  dplyr::select(!c(uniprotswissprot, sample_name)) %>% 
  remove_rownames() %>% 
  column_to_rownames('uniprot_sample')

mrna_prot_df$strand = mrna_prot_df$strand %>% as.factor()

mrna_prot_df$noncds_length = mrna_prot_df$transcript_length - mrna_prot_df$cds_length
mrna_prot_df$proportion_noncds_length = mrna_prot_df$noncds_length / mrna_prot_df$transcript_length

saveRDS(object = mrna_prot_df, file = 'data/whole_raw_dataset.rds')

