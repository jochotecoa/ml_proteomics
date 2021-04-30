# load the library
source('script/functions/functions_JOA.R')
forceLibrary(c('mlbench', 'caret', 'doParallel', 'dplyr', 'RANN'))

cl <- makePSOCKcluster(5)
registerDoParallel(cl)
# load the data

file_X = 'data/training_data_preds_na_omit.rds'
file_Y = 'data/training_data_target_na_omit.rds'

X = readRDS(file = file_X) %>% 
  dplyr::select(!contains('strand'))
Y = readRDS(file = file_Y) %>% 
  unlist()

# dplyr::select(!contains('strand')) # including the strand column gives always an error, independently of the algorithm used

# sizes = seq(1, ncol(X), (ncol(X)-1)/5) %>% as.integer()
sizes = seq(2, ncol(X), 1) %>% as.integer()

path_output = '../output_rfe/na_omit'
if (!dir.exists(path_output)) {
  dir.create(path_output, recursive = T)
}



source('script/recursive_feature_elimination/linear_model_rfe.R',)
source('script/recursive_feature_elimination/glmnet_rfe.R')

source('script/recursive_feature_elimination/bagged_cart_rfe.R')
source('script/recursive_feature_elimination/blackboost_rfe.R')
source('script/recursive_feature_elimination/bstTree_rfe.R')
source('script/recursive_feature_elimination/cubist_rfe.R')
source('script/recursive_feature_elimination/random_forest_rfe.R')

source('script/recursive_feature_elimination/kknn_rfe.R')
source('script/recursive_feature_elimination/knn_rfe.R')
source('script/recursive_feature_elimination/nnet_rfe.R')
source('script/recursive_feature_elimination/neuralnet_rfe.R')
source('script/recursive_feature_elimination/mars_bagEarth_rfe.R')
source('script/recursive_feature_elimination/support_vector_machine_rfe.R')
