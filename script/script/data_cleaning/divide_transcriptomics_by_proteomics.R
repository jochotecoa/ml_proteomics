dupl_mrna = mrna_unip_df$ensembl_transcript_id[mrna_unip_df$ensembl_transcript_id %in% dupl_mrna] %>%
  table() %>%
  as.data.frame()

colnames(dupl_mrna) = c('ensembl_transcript_id', 'n_protein_per_transcript')

mrna_unip_df = mrna_unip_df %>%
  merge.data.frame(y = dupl_mrna, all = T)

mrna_unip_df$n_protein_per_transcript[is.na(mrna_unip_df$n_protein_per_transcript)] = 1

mrna_unip_df_ids = mrna_unip_df[, c(1, ncol(mrna_unip_df)-1, ncol(mrna_unip_df))]

mrna_unip_df = mrna_unip_df[, -c(1, ncol(mrna_unip_df)-1, ncol(mrna_unip_df))] / mrna_unip_df_ids$n_protein_per_transcript
mrna_unip_df = cbind.data.frame(mrna_unip_df, mrna_unip_df_ids)
mrna_unip_df = mrna_unip_df[, -ncol(mrna_unip_df)]
