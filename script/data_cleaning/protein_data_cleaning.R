# Get proteomics data -----------------------------------------------------

Tissue = tools::toTitleCase(tissue)

prot_dir = paste0('/ngs-data/data/hecatos/', Tissue, '/t0_controls/Protein/')

prot_df = mergeFiles(by_col = 'Row.Names', path = prot_dir, header = T, 
                     fill = F, sep = '\t', all_true = T)

# Clean rownames and colnames of proteomics data --------------------------


colnames(prot_df) = colnames(prot_df) %>% 
  gsub(pattern = paste0(prot_dir), replacement = '')

prot_df = prot_df %>% 
  column_to_rownames('Row.Names')

prot_df[, !grepl('log2', colnames(prot_df))] = 
  prot_df[, !grepl('log2', colnames(prot_df))] %>% 
  log2()

colnames(prot_df)[!grepl('log2', colnames(prot_df))] = 
  colnames(prot_df)[!grepl('log2', colnames(prot_df))] %>% 
  paste0('_log2')

if (!file.exists(paste0('data/biostudies/', tissue, '/biostudies_', tissue, '.rds'))) {
  file.remove(paste0('data/biostudies/', tissue, '/biostudies_', tissue, '.rds'))
  source(file = 'script/data_cleaning/renaming_proteomics_samples.R', local = T, echo = T)
}

files_rocheid = readRDS(paste0('data/biostudies/', tissue, '/biostudies_', tissue, '.rds')) %>% 
  as.data.frame() 

colnames(prot_df)[unlist(files_rocheid$colnum)] = files_rocheid$Files

colnames(prot_df) = 
  colnames(prot_df) %>% 
  gsub(pattern = '_/.*', 
       replacement = '')

# Subselect samples -------------------------------------------------------


t0_cols = 
  colnames(prot_df) %>% 
  subset(., grepl(pattern = 'UNTR|DF2|^DMSO|000', .)) 


prot_df = prot_df[, !duplicated(colnames(prot_df))] %>% 
  rownames_to_column() %>% 
  dplyr::select(matches('UNTR|DF2|^DMSO|000'), rowname) %>% 
  column_to_rownames() %>% 
  filterSamplesBySeqDepth() %>% 
  rownames_to_column() %>% 
  cleanProtIds() %>% 
  dplyr::select(!rowname)

prot_df = prot_df %>% 
  remove_rownames() %>% 
  column_to_rownames('uniprot_gn') 

prot_df = prot_df %>% 
  apply(2, unlist) %>% 
  apply(2, as.numeric) %>% 
  data.frame(row.names = rownames(prot_df)) %>% 
  normalizeProteomics() %>% 
  rownames_to_column('uniprotswissprot')

prot_df = prot_df %>% 
  melt.data.frame()

colnames(prot_df)[2:3] = c('sample_name', 'proteomics_value')

# mrna_prot_ids_3_file = paste0('data/mrna_prot_ids_3_', tissue, '.rds')
# 
# if (file.exists()) {
#   mrna_prot_ids_3 = readRDS(file = mrna_prot_ids_3_file)
# } else {
#   mrna_prot_ids_3 = getBM(attributes = c('uniprotswissprot', 'peptide'),
#                           filters = 'uniprotswissprot',
#                           values = unique(prot_df$uniprotswissprot),
#                           mart = mart)
#   
#   mrna_prot_ids_3$peptide_length = mrna_prot_ids_3['peptide'] %>% 
#     sapply(strsplit, split='') %>% 
#     sapply(length)
#   
#   mrna_prot_ids_3 = mrna_prot_ids_3[, -2]
#   
#   mrna_prot_ids_3 = mrna_prot_ids_3[, 'peptide_length'] %>% 
#     aggregate.data.frame(by = list(uniprotswissprot = mrna_prot_ids_3$uniprotswissprot), FUN = median, na.rm = T)
#   colnames(mrna_prot_ids_3)[2] = 'peptide_length'
#   
#   saveRDS(object = mrna_prot_ids_3, file = mrna_prot_ids_3_file)
# }
# 
# prot_df = merge.data.frame(x = prot_df, y = mrna_prot_ids_3, by = 'uniprotswissprot', all.x = T) 

if (tissue == 'cardiac') {
  prot_df$sample_name = prot_df$sample_name %>% 
    gsub(pattern = 'DF2_The', replacement = 'DF2')
  
  # Not have transcrx quantification yet
  prot_df = prot_df %>% 
    dplyr::filter(!grepl(pattern = 'DMSO|UNTR', x = sample_name))
  
}

prot_df$sample_name = prot_df$sample_name %>% 
  gsub(pattern = 'UNTR_The', replacement = 'UNTR') %>% 
  gsub(pattern = 'DMSO_The', replacement = 'DMSO') %>% 
  gsub(pattern = 'APAP', replacement = 'APA') %>% 
  

prot_df$uniprot_sample = paste(prot_df$uniprotswissprot, prot_df$sample_name, sep = '--')

if (tissue == 'hepatic') {
  prot_df$sample_name = prot_df$sample_name %>% 
    gsub(pattern = 'APAP', replacement = 'APA') %>% 
    gsub(pattern = 'DMSO_The_', replacement = 'DMSO_') %>% 
    gsub(pattern = 'UNTR_The_', replacement = 'UNTR_') 
}

# Add protein characteristics as features ---------------------------------

if (file.exists('data/uniprot_yourlist.rds')) {
  uniprot_yourlist = readRDS('data/uniprot_yourlist.rds')
} else {
  source('script/data_cleaning/uniprot_data.R')
  uniprot_yourlist = readRDS('data/uniprot_yourlist.rds')
}

prot_df = prot_df %>% 
  merge.data.frame(y = uniprot_yourlist, by.x = 'uniprotswissprot', by.y = 'Entry', all.x = T)


