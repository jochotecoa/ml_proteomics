# Check that we have the same samples in transcrx and protx ---------------



all(unique(prot_df$sample_name) %in% unique(mrna_unip_df$sample_name)) %>% 
  stopifnot('sample names different between transcrx and protx'= .)

# unique(prot_df$sample_name)[!(unique(prot_df$sample_name) %in% unique(mrna_unip_df$sample_name))]

noprotsamples = unique(mrna_unip_df$sample_name)[!(unique(mrna_unip_df$sample_name) %in% unique(prot_df$sample_name))]
warning(paste(paste0(noprotsamples, collapse = ', '), 'did not have proteomics samples, only transcriptomics'))
mrna_unip_df = mrna_unip_df %>% 
  dplyr::filter(!grepl(pattern = paste0(noprotsamples, collapse = '|'), sample_name))


mrna_unip_df_sum = mrna_unip_df[, 'TPM_value'] %>% 
  aggregate.data.frame(by = list(uniprot_sample = mrna_unip_df$uniprot_sample), FUN = sum, na.rm = T)

colnames(mrna_unip_df_sum)[2] = c('TPM_value')

mrna_unip_df_ids = mrna_unip_df[, c('uniprot_sample', 'uniprotswissprot')] %>% unique()

mrna_unip_df_median = mrna_unip_df[, c('strand', 'transcript_length', 'percentage_gene_gc_content', 'cds_length')] %>% 
  aggregate.data.frame(by = list(uniprotswissprot = mrna_unip_df$uniprotswissprot), FUN = median, na.rm = T)

mrna_unip_df = merge.data.frame(mrna_unip_df_ids, mrna_unip_df_sum, 'uniprot_sample')
mrna_unip_df = merge.data.frame(mrna_unip_df, mrna_unip_df_median, 'uniprotswissprot')
