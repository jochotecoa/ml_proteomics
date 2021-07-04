# load the library
source('script/functions/functions_JOA.R')
forceLibrary(c('mlbench', 'caret', 'doParallel', 'dplyr', 'RANN'))

cl <- makePSOCKcluster(15)
registerDoParallel(cl)
# load the data


# file_X_train = 'data/training_data_preds_na_omit_by_prot.rds'
# file_Y = 'data/training_data_target_na_omit_by_prot.rds'

X = readRDS(file = 'data/whole_data_preds_na_omit.rds')
Y = readRDS(file = 'data/whole_data_target_na_omit.rds') %>% 
  unlist()


# dplyr::select(!contains('strand'))


# dplyr::select(!contains('strand')) # including the strand column gives always an error, independently of the algorithm used

# sizes = seq(1, ncol(X), (ncol(X)-1)/5) %>% as.integer()

path_output = '../output_rfe/split_by_sample/na_omit/cardiac'
if (!dir.exists(path_output)) {
  dir.create(path_output, recursive = T)
}

train_samples = X_clean %>% rownames() %>% strsplit('--') %>% sapply('[[', 2)

sizes = c(1:10, seq(15, ncol(X), by = 5))

folds = groupKFold(train_samples, k = 10)

repeats = 5
