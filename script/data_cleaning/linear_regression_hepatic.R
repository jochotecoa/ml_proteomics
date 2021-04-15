source('script/functions/functions_JOA.R')
forceLibrary(c('dplyr', 'tibble', 'biomaRt', 'reshape', 'caret', 'caTools'))


# Get mRNA data -----------------------------------------------------------

mrna_dir = '/share/analysis/hecatos/juantxo/mRNA/quant_salmon/Homo_sapiens.GRCh38.cdna.ncrna.circbase/'
tissue = 'hepatic'
compound = 't0_controls_ML'

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


prot_df$uniprot_sample = paste(prot_df$uniprotswissprot, prot_df$sample_name, sep = '--')

if (tissue == 'hepatic') {
  prot_df$sample_name = prot_df$sample_name %>% 
    gsub(pattern = 'APAP', replacement = 'APA') %>% 
    gsub(pattern = 'DMSO_The_', replacement = 'DMSO_') %>% 
    gsub(pattern = 'UNTR_The_', replacement = 'UNTR_') 
}

if (file.exists('data/uniprot_yourlist.rds')) {
  uniprot_yourlist = readRDS('data/uniprot_yourlist.rds')
} else {
  source('script/data_cleaning/uniprot_data.R')
  uniprot_yourlist = readRDS('data/uniprot_yourlist.rds')
}

prot_df = prot_df %>% 
  merge.data.frame(y = uniprot_yourlist, by.x = 'uniprotswissprot', by.y = 'Entry', all.x = T)


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


mrna_unip_df_sum = mrna_unip_df[, 'TPM_value'] %>% 
  aggregate.data.frame(by = list(uniprot_sample = mrna_unip_df$uniprot_sample), FUN = sum, na.rm = T)

colnames(mrna_unip_df_sum)[2] = c('TPM_value')

mrna_unip_df_ids = mrna_unip_df[, c('uniprot_sample', 'uniprotswissprot')] %>% unique()

mrna_unip_df_median = mrna_unip_df[, c('strand', 'transcript_length', 'percentage_gene_gc_content', 'cds_length')] %>% 
  aggregate.data.frame(by = list(uniprotswissprot = mrna_unip_df$uniprotswissprot), FUN = median, na.rm = T)

mrna_unip_df = merge.data.frame(mrna_unip_df_ids, mrna_unip_df_sum, 'uniprot_sample')
mrna_unip_df = merge.data.frame(mrna_unip_df, mrna_unip_df_median, 'uniprotswissprot')


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

# Pre-processing ----------------------------------------------------------
# Zero- and Near Zero-Variance Predictors

nzv_m <- nearZeroVar(mrna_prot_df, saveMetrics = T)
nzv <- nearZeroVar(mrna_prot_df)

if (length(nzv) > 0) {
  zv_colnames = colnames(mrna_prot_df)[nzv_m$zeroVar]
  nzv_colnames = colnames(mrna_prot_df)[as.logical(nzv_m$nzv - nzv_m$zeroVar)]
  nzv_m[nzv, ] %>% print()
  # warning(paste(paste0(zv_colnames, collapse = ', '), 'presented zero variance'))
  # warning(paste(paste0(nzv_colnames, collapse = ', '), 'presented near zero variance'))
  mrna_prot_df <- mrna_prot_df[, -nzv]
}

# Creating Dummy Variables

dummies <- dummyVars(" ~ .", data = mrna_prot_df)
mrna_prot_df = predict(dummies, newdata = mrna_prot_df)
mrna_prot_df = mrna_prot_df %>% as.data.frame()


# Handle Missing Values ---------------------------------------------------

feat_missval = mrna_prot_df %>% apply(2, anyNA) %>% subset(., .) %>% names()
miss_per_feature = mrna_prot_df[, feat_missval] %>% summary %>% .[7, ] %>% 
  gsub(pattern = "NA's   :| ", replacement = "") %>% unlist() %>% as.numeric()
if (any(miss_per_feature >= (0.05*nrow(mrna_prot_df)))) {
  stop()
} else {
  prev_rows = nrow(mrna_prot_df)
  mrna_prot_df = mrna_prot_df %>% na.omit()
  warning(paste(prev_rows - (nrow(mrna_prot_df)), 'rows have been removed due to missing values'))
}


# Identifying Correlated Predictors

descrCor <- cor(mrna_prot_df)
highlyCor <- findCorrelation(descrCor, cutoff = .75)
if (length(highlyCor) > 0) {
  mrna_prot_df <- mrna_prot_df[,-highlyCor]
}


# Linear Dependencies

comboInfo <- findLinearCombos(mrna_prot_df)
if (length(comboInfo$remove) > 0) {
  mrna_prot_df = mrna_prot_df[, -comboInfo$remove]
}

mrna_prot_df = mrna_prot_df %>% as.data.frame()

saveRDS(object = mrna_prot_df, file = 'data/training_data.rds')

# 

trainIndex = createDataPartition(mrna_prot_df$proteomics_value, p=.8, list=F, 
                                 times = 1)

protTrain <- mrna_prot_df[ trainIndex,]
protTest  <- mrna_prot_df[-trainIndex,]

# Model building ----------------------------------------------------------
# [1] "Linear Regression (lm)"                                             "eXtreme Gradient Boosting (xgbDART)"                               
# [3] "Quantile Random Forest (qrf)"                                       "Relevance Vector Machines with Linear Kernel (rvmLinear)"          
# [5] "Partial Least Squares (pls)"                                        "Bayesian Regularized Neural Networks (brnn)"                       
# [7] "glmnet (glmnet)"                                                    "Adaptive-Network-Based Fuzzy Inference System (ANFIS)"             
# [9] "Gaussian Process with Radial Basis Function Kernel (gaussprRadial)" "Linear Regression with Stepwise Selection (lmStepAIC)"             
# [11] "Robust Linear Model (rlm)"                                         


library(caret)# Simple linear regression model (lm means linear model)


fitControl <- trainControl(method = "repeatedcv",   
                           number = 10,     # number of folds
                           repeats = 3, 
                           search = "random", 
                           allowParallel= T)    # repeated ten times

methods_caret = c('lm', 'xgbDART', 'qrf', 'rvmLinear', 'pls', 'brnn', 'glmnet', 'ANFIS', 'gaussprRadial', 'lmStepAIC', 'rlm')

models_list = list()

for (method_caret in methods_caret) {
  model <- train(proteomics_value ~ .,
                 data = mrna_prot_df,
                 method = method_caret, 
                 trControl = fitControl, 
                 preProcess = c('scale', 'center'), 
                 tuneLength = 10)
  
  models_list$new = model
  names(models_list)[length(models_list)] = method_caret
  
}

model_lm <- train(proteomics_value ~ log2_TPM_value,
               data = mrna_prot_df,
               method = "lm", 
               trControl = fitControl, 
               preProcess = c('scale', 'center'))

model_mlm <- train(proteomics_value ~ .,
                   data = mrna_prot_df,
                   method = "lm",
                   trControl = fitControl, preProcess = c('scale', 'center'))


# possible values: boot", "boot632", "cv", "repeatedcv", "LOOCV", "LGOCV"

model_lasso <- train(proteomics_value ~ .,
                  data = mrna_prot_df,
                  method = "lasso",  # now we're using the lasso method
                  trControl = fitControl, preProcess = c('scale', 'center'))


# Multiple linear regression model


# Ridge regression model
model_ridge <- train(proteomics_value ~ .,
               data = mrna_prot_df,
               method = "ridge", # Try using "lasso"
               trControl = fitControl, preProcess = c('scale', 'center'))


resamps <- resamples(list(LM = model_lm,
                          MLM = model_mlm,
                          LSS = model_lasso,
                          RDG = model_ridge))

theme1 <- trellis.par.get()
theme1$plot.symbol$col = rgb(.2, .2, .2, .4)
theme1$plot.symbol$pch = 16
theme1$plot.line$col = rgb(1, 0, 0, .7)
theme1$plot.line$lwd <- 2
trellis.par.set(theme1)
bwplot(resamps, layout = c(4, 1))

summary(model_lm)$sigm
min(model_lasso$results$RMSE)
summary(model_mlm)$sigma
min(model_ridge$results$RMSE)

rmse_matrix = cbind(
  rmse_matrix, 
  c(
    summary(model_lm)$sigma,
    min(model_lasso$results$RMSE),
    summary(model_mlm)$sigma,
    min(model_ridge$results$RMSE)
  )
)

# rmse_matrix = matrix(
#   c(
#     summary(model_lm)$sigma,
#     min(model_lasso$results$RMSE),
#     summary(model_mlm)$sigma,
#     min(model_ridge$results$RMSE)
#   )
# )

# rownames(rmse_matrix) = c('linear', 'lasso', 'multiple linear', 'ridge')
# barplot(t(rmse_matrix)) 

ggplot(rmse_matrix %>% t %>% melt, aes(x = X1, y = value, col = X2)) +           # Draw line plot with ggplot2
  geom_line()
