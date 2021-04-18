# Get mRNA data -----------------------------------------------------------

mrna_dir = '/share/analysis/hecatos/juantxo/mRNA/quant_salmon/Homo_sapiens.GRCh38.cdna.ncrna.circbase/'

mrna_dir_spe = paste0(mrna_dir, tissue, '/', compound)

mrna_df_file = paste0('../bigdata_ml_proteomics/salmon_', tissue, '.rds')

if (file.exists(mrna_df_file)) {
  mrna_df = readRDS(file = mrna_df_file)
} else {
  mrna_df = mergeFiles(files_patt = 'quant.sf', by_col = 'Name',
                       path = mrna_dir_spe, all_true = T, recursive = T,
                       header = T)
  
  saveRDS(object = mrna_df, file = mrna_df_file)
  
}


# Clean colnames and rownames ---------------------------------------------


mrna_df_counts = mrna_df %>% 
  remove_rownames() %>% 
  column_to_rownames('Name') %>% 
  dplyr::select(contains('NumReads'))

colnames(mrna_df_counts) = colnames(mrna_df_counts) %>% 
  gsub(pattern = paste0(mrna_dir_spe, '|', '/quant.sf|_quant|NumReads_', '|', '/'), replacement = '') 

mrna_df_counts = mrna_df_counts %>% 
  filterSamplesBySeqDepth()

mrna_df = mrna_df %>% 
  remove_rownames() %>% 
  column_to_rownames('Name') %>% 
  dplyr::select(contains('TPM'))

colnames(mrna_df) = colnames(mrna_df) %>% 
  gsub(pattern = mrna_dir_spe, replacement = '') %>% 
  gsub(pattern = '/quant.sf|_quant|TPM_', replacement = '') %>% 
  gsub(pattern = '/', replacement = '')


rownames(mrna_df) = rownames(mrna_df) %>% 
  gsub('\\..*', '', .) 

stopifnot(all.equal(colnames(mrna_df), colnames(mrna_df_counts)))


# Get protein IDs associated to transcripts -------------------------------



mrna_prot_ids_file = paste0('data/mrna_prot_ids_', tissue, '.rds')

if (file.exists(mrna_prot_ids_file)) {
  mrna_prot_ids = readRDS(file = mrna_prot_ids_file)
} else {
  mart = openMart2018()
  
  mrna_prot_ids = getBM(attributes = c('ensembl_transcript_id', 'uniprotswissprot', 'strand', 'transcript_length', 'percentage_gene_gc_content'),
                        filters = 'ensembl_transcript_id',
                        values = rownames(mrna_df),
                        mart = mart)
  
  mrna_prot_ids_2 = getBM(attributes = c('ensembl_transcript_id', 'cds_length'),
                          filters = 'ensembl_transcript_id',
                          values = rownames(mrna_df),
                          mart = mart)
  
  
  
  mrna_prot_ids_4 = merge.data.frame(x = mrna_prot_ids, y = mrna_prot_ids_2, by = 'ensembl_transcript_id', all = T)
  
  
  saveRDS(object = mrna_prot_ids_4, file = mrna_prot_ids_file)
  
  mrna_prot_ids = mrna_prot_ids_4
  
}

mrna_unip_df = mrna_df %>% 
  rownames_to_column('ensembl_transcript_id') %>% 
  merge.data.frame(y = mrna_prot_ids, by = 'ensembl_transcript_id', all.y = T) 

