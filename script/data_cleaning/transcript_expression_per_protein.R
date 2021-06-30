# Divide transcripts per protein (after connecting to protein expr --------

mrna_unip_df = mrna_unip_df[mrna_unip_df$uniprotswissprot %in% prot_df$uniprotswissprot, ]

dupl_mrna = mrna_unip_df$ensembl_transcript_id %>% 
  subset(., duplicated(.))

mrna_unip_df = mrna_unip_df[!(mrna_unip_df$ensembl_transcript_id %in% dupl_mrna), ]

mrna_unip_df = mrna_unip_df %>% 
  melt.data.frame(id = c('ensembl_transcript_id', 'uniprotswissprot', 'strand', 'transcript_length', 'percentage_gene_gc_content', 'cds_length'))

colnames(mrna_unip_df)[(ncol(mrna_unip_df)-1):ncol(mrna_unip_df)] = c('sample_name', 'TPM_value')

if (tissue == 'cardiac') {
  mrna_unip_df$sample_name = mrna_unip_df$sample_name %>% 
    gsub(pattern = '5FU', replacement = 'X5FU') %>% 
    gsub(pattern = 'con_DF2', replacement = 'DF2')
  
  mrna_unip_df = mrna_unip_df %>% 
    dplyr::filter(!grepl(pattern = 'DAU', x = sample_name))
  
  mrna_unip_df$sample_name = mrna_unip_df$sample_name %>% 
    gsub(pattern = 'Con', replacement = '')
  
}

if (tissue == 'hepatic') {
  mrna_unip_df$sample_name = mrna_unip_df$sample_name %>% 
    gsub(pattern = 'Con', replacement = '')
}


mrna_unip_df$uniprot_sample = paste(mrna_unip_df$uniprotswissprot, mrna_unip_df$sample_name, sep = '--')


