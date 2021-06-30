addVarsProt <- function(x, fnc_list, by_str) {
  fnc_str = fnc_list[1]
  fnc = get(fnc_str)
  transf_cols = colnames(x)[!grepl(by_str, colnames(x))]
  by_lst = x[, by_str] %>% list()
  names(by_lst) = by_str
  df = x %>% 
    .[, !grepl(by_str, colnames(x)), F] %>% 
    aggregate.data.frame(by = by_lst, FUN = fnc, na.rm = TRUE) 
  colnames(df)[!grepl(by_str, colnames(df))] = paste0(transf_cols, '_', fnc_str)
  if (length(fnc_list) > 1) {
    for (fnc_str in fnc_list[-1]) {
      fnc = get(fnc_str)
      transf_cols = colnames(x)[!grepl(by_str, colnames(x))]
      by_lst = x[, by_str] %>% list()
      names(by_lst) = by_str
      x_2 = x %>% 
        .[, !grepl(by_str, colnames(x)), F] %>% 
        aggregate.data.frame(by = by_lst, FUN = fnc, na.rm = TRUE) 
      colnames(x_2)[!grepl(by_str, colnames(x_2))] = paste0(transf_cols, '_', fnc_str)
      df = x_2 %>% 
        merge.data.frame(x = df, by = by_str)
    }
  }
  return(df)
}

# Check that we have the same samples in transcrx and protx ---------------



all(unique(prot_df$sample_name) %in% unique(mrna_unip_df$sample_name)) %>% 
  stopifnot('sample names different between transcrx and protx'= .)

# unique(prot_df$sample_name)[!(unique(prot_df$sample_name) %in% unique(mrna_unip_df$sample_name))]

noprotsamples = unique(mrna_unip_df$sample_name)[!(unique(mrna_unip_df$sample_name) %in% unique(prot_df$sample_name))]
if (length(noprotsamples) > 0) {
  warning(paste(paste0(noprotsamples, collapse = ', '), 'did not have proteomics samples, only transcriptomics'))
  mrna_unip_df = mrna_unip_df %>% 
    dplyr::filter(!grepl(pattern = paste0(noprotsamples, collapse = '|'), sample_name))
}

notranscrsamples = unique(prot_df$sample_name)[!(unique(prot_df$sample_name) %in% unique(mrna_unip_df$sample_name))]
if (length(notranscrsamples) > 0) {
  warning(paste(paste0(notranscrsamples, collapse = ', '), 'did not have transcriptomics samples, only proteomics'))
  prot_df = prot_df %>% 
    dplyr::filter(!grepl(pattern = paste0(notranscrsamples, collapse = '|'), sample_name))
  
}

# mrna_unip_df_sum = mrna_unip_df[, 'TPM_value'] %>% 
#   aggregate.data.frame(by = list(uniprot_sample = mrna_unip_df$uniprot_sample), FUN = sum, na.rm = T)

mrna_unip_df_sum = mrna_unip_df[, c('uniprot_sample', 'TPM_value')] %>% 
  addVarsProt(fnc_list = c('mean', 'median', 'min', 'max', 'sum', 'sd'), by_str = 'uniprot_sample')

mrna_unip_df_sum_log2 = mrna_unip_df_sum[, -grep('uniprot_sample', colnames(mrna_unip_df_sum))] %>% 
  log2()
colnames(mrna_unip_df_sum_log2) = colnames(mrna_unip_df_sum_log2) %>% 
  paste0('_log2')

mrna_unip_df_sum = mrna_unip_df_sum %>% 
  cbind.data.frame(mrna_unip_df_sum_log2)

mrna_unip_df_sum[mrna_unip_df_sum == -Inf] = log2(2e-06)


mrna_unip_df_ids = mrna_unip_df[, c('uniprot_sample', 'uniprotswissprot')] %>% unique()

mrna_unip_df$noncds_length = mrna_unip_df$transcript_length - mrna_unip_df$cds_length
mrna_unip_df$proportion_noncds_length = mrna_unip_df$noncds_length / mrna_unip_df$transcript_length

mrna_unip_df_median = mrna_unip_df[, c('uniprotswissprot', 'strand', 
                                       'transcript_length', 
                                       'percentage_gene_gc_content', 
                                       'cds_length', 'noncds_length', 
                                       'proportion_noncds_length')] %>% 
  addVarsProt(fnc_list = c('mean', 'median', 'min', 'max', 'sum', 'sd'), 
              by_str = 'uniprotswissprot')

mrna_unip_df = merge.data.frame(mrna_unip_df_ids, mrna_unip_df_sum, 'uniprot_sample')
mrna_unip_df = merge.data.frame(mrna_unip_df, mrna_unip_df_median, 'uniprotswissprot')

mrna_unip_df[, 'TPM_value_sd'] = mrna_unip_df[, 'TPM_value_sd'] %>% 
  naToZero()
mrna_unip_df[, 'TPM_value_sd_log2'][is.na(mrna_unip_df[, 'TPM_value_sd_log2'])] = log2(2e-06)

