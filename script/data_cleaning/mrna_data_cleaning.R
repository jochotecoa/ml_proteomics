source('script/functions/functions_JOA.R')
forceLibrary(c('dplyr', 'tibble', 'biomaRt', 'reshape'))


# Get mRNA data -----------------------------------------------------------

mrna_dir = '/share/analysis/hecatos/juantxo/mRNA/quant_salmon/Homo_sapiens.GRCh38.cdna.ncrna.circbase/'
tissue = 'cardiac'
compound = 't0_controls'

mrna_dir_spe = paste0(mrna_dir, tissue, '/', compound)

mrna_df = mergeFiles(files_patt = 'quant.sf', by_col = 'Name', 
                     path = mrna_dir_spe, all_true = T, recursive = T, 
                     header = T)

colnames(mrna_df) = colnames(mrna_df) %>% 
  gsub(pattern = '/share/analysis/hecatos/juantxo/mRNA/quant_salmon/Homo_sapiens.GRCh38.cdna.ncrna.circbase/cardiac/t0_controls/', replacement = '') %>% 
  gsub(pattern = '/quant.sf|_quant|TPM', replacement = '')

mrna_df = mrna_df %>% 
  column_to_rownames('Name') %>% 
  select(contains('TPM'))

rownames(mrna_df) = rownames(mrna_df) %>% 
  gsub('\\..*', '', .) 

mart = openMart2018()

mrna_prot_ids = getBM(attributes = c('ensembl_transcript_id', 'uniprot_gn'), 
                      filters = 'ensembl_transcript_id', 
                      values = rownames(mrna_df), 
                      mart = mart)

mrna_unip_df = mrna_df %>% 
  rownames_to_column('ensembl_transcript_id') %>% 
  merge.data.frame(y = mrna_prot_ids, by = 'ensembl_transcript_id', all.y = T)


# Get proteomics data -----------------------------------------------------

prot_dir = '/ngs-data/data/hecatos/Cardiac/t0_controls/Protein/'

# source(file = 'script/data_cleaning/fix_ANT_proteomics_file.R')

prot_df = mergeFiles(by_col = 'Row.Names', path = prot_dir, header = T, 
                     fill = F, sep = '\t', all_true = T)

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

# source(file = 'script/data_cleaning/renaming_proteomics_samples.R')

files_rocheid = readRDS('data/biostudies/cardiac/biostudies_cardiac.rds') %>% 
  as.data.frame()

files_rocheid$`Roche ID`[grep('AMI_.*_000', files_rocheid$Files)] = 
  files_rocheid$`Roche ID`[grep('AMI_.*_000', files_rocheid$Files)] %>% 
  paste0('_e1')

files_rocheid$`Roche ID`[grep('DOC_.*_000', files_rocheid$Files)] = 
  files_rocheid$`Roche ID`[grep('DOC_.*_000', files_rocheid$Files)] %>% 
  paste0('_e2')

files_rocheid$colnum = files_rocheid$`Roche ID` %>% 
  paste0('_', ., '_') %>% 
  as.data.frame() %>% 
  apply(MARGIN = 1, FUN = grep, colnames(prot_df)) 

files_rocheid$len = files_rocheid$colnum %>% 
  sapply(length)

files_rocheid$colnum[files_rocheid$len == 0] = 
  files_rocheid$`Roche ID`[files_rocheid$len == 0] %>% 
  paste0('_0', ., '_') %>% 
  as.data.frame() %>% 
  apply(MARGIN = 1, FUN = grep, colnames(prot_df)) 

files_rocheid$len = files_rocheid$colnum %>% 
  sapply(length)


files_rocheid = files_rocheid[files_rocheid$len != 0, ]

colnames(prot_df)[unlist(files_rocheid$colnum)] = files_rocheid$Files

colnames(prot_df) = 
  colnames(prot_df) %>% 
  gsub(pattern = '_/.*', 
       replacement = '')

t0_cols = 
  colnames(prot_df) %>% 
  subset(., grepl(pattern = 'UNTR|DF2|^DMSO|000', .)) 


prot_df = prot_df[, !duplicated(colnames(prot_df))] %>% 
  dplyr::select(matches('UNTR|DF2|^DMSO|000')) %>% 
  cleanProtIds() %>% 
  dplyr::select(!rowname) 

prot_df = prot_df %>% 
  remove_rownames() %>% 
  column_to_rownames('uniprot_gn') 

prot_df = 
  prot_df %>% 
  apply(2, unlist) %>% 
  apply(2, as.numeric) %>% 
  data.frame(row.names = rownames(prot_df)) %>% 
  normalizeProteomics() %>% 
  rownames_to_column('uniprot_gn')



prot_df = 
  prot_df %>% 
  melt.data.frame() %>% 
  head()


# Combine mRNA with proteins ----------------------------------------------

mrna_unip_df = mrna_unip_df %>% 
  merge.data.frame(y = prot_df)

# Divide transcripts per protein (after connecting to protein expr --------

dupl_mrna = mrna_unip_df$ensembl_transcript_id %>% 
  subset(., duplicated(.))

dupl_mrna = mrna_unip_df$ensembl_transcript_id[mrna_unip_df$ensembl_transcript_id %in% dupl_mrna] %>% 
  table() %>% 
  as.data.frame()

colnames(dupl_mrna) = c('ensembl_transcript_id', 'n_protein_per_transcript')

mrna_unip_df = mrna_unip_df %>% 
  merge.data.frame(y = dupl_mrna, all = T)

mrna_unip_df$n_protein_per_transcript[is.na(mrna_unip_df$n_protein_per_transcript)] = 1
