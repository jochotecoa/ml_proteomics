source('script/functions/functions_JOA.R')
forceLibrary(c('dplyr', 'tibble', 'biomaRt', 'reshape', 'caret'))


# Get mRNA data -----------------------------------------------------------

mrna_dir = '/share/analysis/hecatos/juantxo/mRNA/quant_salmon/Homo_sapiens.GRCh38.cdna.ncrna.circbase/'
tissue = 'hepatic'
compound = 't0_controls_ML'

mrna_dir_spe = paste0(mrna_dir, tissue, '/', compound)

mrna_df_file = paste0('data/salmon_', tissue, '.rds')

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
  
  mrna_prot_ids_3 = getBM(attributes = c('uniprotswissprot', 'peptide'),
                          filters = 'ensembl_transcript_id',
                          values = rownames(mrna_df),
                          mart = mart)
  
  mrna_prot_ids_4 = merge.data.frame(x = mrna_prot_ids, y = mrna_prot_ids_2, by = 'ensembl_transcript_id', all = T)
  
  # mrna_prot_ids_5 = merge.data.frame(x = mrna_prot_ids_4, y = mrna_prot_ids_3, by = 'uniprotswissprot')
  
  saveRDS(object = mrna_prot_ids_4, file = mrna_prot_ids_file)
  
}

mrna_unip_df = mrna_df %>% 
  rownames_to_column('ensembl_transcript_id') %>% 
  merge.data.frame(y = mrna_prot_ids, by = 'ensembl_transcript_id', all.y = T) 

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

if (tissue == 'cardiac') {
  prot_df$sample_name = prot_df$sample_name %>% 
    gsub(pattern = 'DF2_The', replacement = 'DF2')
  
  # Not have transcrx quantification yet
  prot_df = prot_df %>% 
    dplyr::filter(!grepl(pattern = 'DMSO|UNTR', x = sample_name))
  
}


prot_df$uniprot_sample = paste(prot_df$uniprotswissprot, prot_df$sample_name, sep = '--')

if (tissue == 'hepatic') {
  prot_df$sample_name = prot_df$sample_name %>% 
    gsub(pattern = 'APAP', replacement = 'APA') %>% 
    gsub(pattern = 'DMSO_The_', replacement = 'DMSO_') %>% 
    gsub(pattern = 'UNTR_The_', replacement = 'UNTR_') 
}


# Divide transcripts per protein (after connecting to protein expr --------

mrna_unip_df = mrna_unip_df[mrna_unip_df$uniprotswissprot %in% prot_df$uniprotswissprot, ]

dupl_mrna = mrna_unip_df$ensembl_transcript_id %>% 
  subset(., duplicated(.))

mrna_unip_df = mrna_unip_df[!(mrna_unip_df$ensembl_transcript_id %in% dupl_mrna), ]

# dupl_mrna = mrna_unip_df$ensembl_transcript_id[mrna_unip_df$ensembl_transcript_id %in% dupl_mrna] %>% 
#   table() %>% 
#   as.data.frame()
# 
# colnames(dupl_mrna) = c('ensembl_transcript_id', 'n_protein_per_transcript')
# 
# mrna_unip_df = mrna_unip_df %>% 
#   merge.data.frame(y = dupl_mrna, all = T)
# 
# mrna_unip_df$n_protein_per_transcript[is.na(mrna_unip_df$n_protein_per_transcript)] = 1
# 
# mrna_unip_df_ids = mrna_unip_df[, c(1, ncol(mrna_unip_df)-1, ncol(mrna_unip_df))]
# 
# mrna_unip_df = mrna_unip_df[, -c(1, ncol(mrna_unip_df)-1, ncol(mrna_unip_df))] / mrna_unip_df_ids$n_protein_per_transcript
# mrna_unip_df = cbind.data.frame(mrna_unip_df, mrna_unip_df_ids)
# mrna_unip_df = mrna_unip_df[, -ncol(mrna_unip_df)]


mrna_unip_df = mrna_unip_df %>% 
  melt.data.frame(id = c('ensembl_transcript_id', 'uniprotswissprot', 'strand', 'transcript_length', 'percentage_gene_gc_content', 'cds_length'))

colnames(mrna_unip_df)[ncol(mrna_unip_df)-1:ncol(mrna_unip_df)] = c('sample_name', 'TPM_value')

if (tissue == 'cardiac') {
  mrna_unip_df$sample_name = mrna_unip_df$sample_name %>% 
    gsub(pattern = '5FU', replacement = 'X5FU') %>% 
    gsub(pattern = 'con_DF2', replacement = 'DF2')
  
  mrna_unip_df = mrna_unip_df %>% 
    dplyr::filter(!grepl(pattern = 'DAU', x = sample_name))
  
}

if (tissue == 'hepatic') {
  mrna_unip_df$sample_name = mrna_unip_df$sample_name %>% 
    gsub(pattern = 'Con', replacement = '')
}


mrna_unip_df$uniprot_sample = paste(mrna_unip_df$uniprotswissprot, mrna_unip_df$sample_name, sep = '--')


# Check that we have the same samples in transcrx and protx ---------------



all(unique(prot_df$sample_name) %in% unique(mrna_unip_df$sample_name)) %>% 
  stopifnot('sample names different between transcrx and protx'= .)

# unique(prot_df$sample_name)[!(unique(prot_df$sample_name) %in% unique(mrna_unip_df$sample_name))]

noprotsamples = unique(mrna_unip_df$sample_name)[!(unique(mrna_unip_df$sample_name) %in% unique(prot_df$sample_name))]
warning(paste(paste0(noprotsamples, collapse = ', '), 'did not have proteomics samples, only transcriptomics'))
mrna_unip_df = mrna_unip_df %>% 
  dplyr::filter(!grepl(pattern = paste0(noprotsamples, collapse = '|'), sample_name))


mrna_unip_df = mrna_unip_df[, 'TPM_value'] %>% 
  aggregate.data.frame(by = list(uniprot_sample = mrna_unip_df$uniprot_sample), FUN = sum)

colnames(mrna_unip_df)[2] = c('TPM_value')

# Combine mRNA with proteins ----------------------------------------------

mrna_prot_df = mrna_unip_df %>% 
  merge.data.frame(y = prot_df, by = 'uniprot_sample', all.y = T)

mrna_prot_df = mrna_prot_df %>% 
  dplyr::select(!contains('.y'))

colnames(mrna_prot_df) = colnames(mrna_prot_df) %>% 
  gsub(pattern = '.x', replacement = '')

# mrna_prot_df = mrna_prot_df %>% 
#   dplyr::filter(!is.na(proteomics_value))

mrna_prot_df = mrna_prot_df %>%
  zeroToNa() %>% 
  na.omit()
# Model building ----------------------------------------------------------

# possible values: boot", "boot632", "cv", "repeatedcv", "LOOCV", "LGOCV"
fitControl <- trainControl(method = "repeatedcv",   
                           number = 10,     # number of folds
                           repeats = 10)    # repeated ten times

model.cv <- train(mpg ~ .,
                  data = mtcars,
                  method = "lasso",  # now we're using the lasso method
                  trControl = fitControl)  



mrna_prot_df$log2_TPM_value = log2(mrna_prot_df$TPM_value)

library(caret)# Simple linear regression model (lm means linear model)
model <- train(proteomics_value ~ log2_TPM_value,
               data = mrna_prot_df,
               method = "lm")

# 1. Add predictions 
pred.int <- predict(model, interval = "prediction")
mydata <- cbind(mrna_prot_df, pred.int)

# 2. Regression line + confidence intervals
library("ggplot2")
ggplot(mydata, aes(log2(TPM_value), proteomics_value)) +
  geom_point() +
  geom_line(aes(y = pred.int), color = "red")
  stat_smooth(method = lm)


mse <- mean((mydata$proteomics_value - mydata$pred.int)^2)
rmse <- sqrt(mse) # Residual standard error in summary(model)





# Multiple linear regression model
model <- train(mpg ~ .,
               data = mtcars,
               method = "lm")

# Ridge regression model
model <- train(mpg ~ .,
               data = mtcars,
               method = "ridge") # Try using "lasso"
