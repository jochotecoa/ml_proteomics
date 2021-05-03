# load the library
source('script/functions/functions_JOA.R')
forceLibrary(c('mlbench', 'caret', 'doParallel', 'dplyr', 'RANN'))

cl <- makePSOCKcluster(5)
registerDoParallel(cl)
# load the data


file_X_train = 'data/training_data_preds.rds'
file_Y = 'data/training_data_target.rds'

X_train = readRDS(object = , file = 'data/training_data_preds.rds')
Y_train = readRDS(object = , file = 'data/training_data_target.rds')

X_test = readRDS(object = , file = 'data/test_data_preds.rds')
Y_test = readRDS(object = , file = 'data/test_data_target.rds')


X = X_train %>% 
  dplyr::select(!contains('strand'))
Y = Y_train %>% 
  unlist()


# dplyr::select(!contains('strand')) # including the strand column gives always an error, independently of the algorithm used

# sizes = seq(1, ncol(X), (ncol(X)-1)/5) %>% as.integer()
sizes = seq(1, ncol(X), 1) %>% as.integer()

path_output = '../output_rfe/data_split'
if (!dir.exists(path_output)) {
  dir.create(path_output, recursive = T)
}


