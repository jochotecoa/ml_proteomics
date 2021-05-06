source('script/recursive_feature_elimination/load_data_split_rfe_na_omit.R.R')

source('script/recursive_feature_elimination/linear_model_rfe.R')
source('script/recursive_feature_elimination/bagged_cart_rfe.R')
source('script/recursive_feature_elimination/glmnet_rfe.R')
source('script/recursive_feature_elimination/kknn_rfe.R')
source('script/recursive_feature_elimination/random_forest_rfe.R')
source('script/recursive_feature_elimination/blackboost_rfe.R')
source('script/recursive_feature_elimination/cubist_rfe.R')
source('script/recursive_feature_elimination/bstTree_rfe.R')

# source('script/recursive_feature_elimination/knn_rfe.R')
# source('script/recursive_feature_elimination/mars_bagEarth_rfe.R')
source('script/recursive_feature_elimination/nnet_rfe.R')
source('script/recursive_feature_elimination/support_vector_machine_rfe.R')
source('script/recursive_feature_elimination/neuralnet_rfe.R')

stopCluster(cl)

methods = c('bag', 
            # 'blackboost', 'bstTree', 'cubist', 'glmnet', 'kknn', 
            'lm', 'rf')
            # 'knn', 'bagEarth', 'neuralnet', 'nnet')
results = list()

for (method_i in methods) {
  file_rds = paste0(path_output, '/', method_i, 'Profile.rds')
  prof = readRDS(file_rds)
  if (method_i == 'rf') {
    prot_pred = predict(prof$fit, X_test, type="response",
                        norm.votes=TRUE, predict.all=FALSE, proximity=FALSE, nodes=FALSE)
    
  } else {
    prot_pred = predict(prof$fit, X_test)
  }
  
  result = postResample(pred = prot_pred, obs = Y_test)
  print(method_i)
  print(result)
  results = c(results, result)
}
