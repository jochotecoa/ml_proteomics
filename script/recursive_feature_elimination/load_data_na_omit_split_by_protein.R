# load the library
source('script/functions/functions_JOA.R')
forceLibrary(c('mlbench', 'caret', 'doParallel', 'dplyr', 'RANN'))

# cl <- makePSOCKcluster(15)
# registerDoParallel(cl)
# load the data


# file_X_train = 'data/training_data_preds_na_omit_by_prot.rds'
# file_Y = 'data/training_data_target_na_omit_by_prot.rds'

X = readRDS(file = 'data/whole_data/preds_na_omit_hepatic.rds')
Y = readRDS(file = 'data/whole_data/target_na_omit_hepatic.rds') %>% 
  unlist()

 
  # dplyr::select(!contains('strand'))


# dplyr::select(!contains('strand')) # including the strand column gives always an error, independently of the algorithm used

# sizes = seq(1, ncol(X), (ncol(X)-1)/5) %>% as.integer()
sizes = as.integer(seq(2, ncol(X), floor(ncol(X)/4)))[-5]

path_output = '../output_rfe/split_by_prot/na_omit/hepatic/72_features'
if (!dir.exists(path_output)) {
  dir.create(path_output, recursive = T)
}

train_prots = X %>% rownames() %>% strsplit('--') %>% sapply('[[', 1)

folds = groupKFold(train_prots, k = 10)

