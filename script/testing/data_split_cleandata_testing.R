source('script/recursive_feature_elimination/load_data_split_rfe_na_omit.R')

source('script/recursive_feature_elimination/linear_model_rfe.R')
source('script/recursive_feature_elimination/glmnet_rfe.R')
source('script/recursive_feature_elimination/bagged_cart_rfe.R')
source('script/recursive_feature_elimination/kknn_rfe.R')
source('script/recursive_feature_elimination/random_forest_rfe.R')
source('script/recursive_feature_elimination/blackboost_rfe.R')
source('script/recursive_feature_elimination/cubist_rfe.R')
source('script/recursive_feature_elimination/support_vector_machine_rfe.R')
source('script/recursive_feature_elimination/bstTree_rfe.R')
source('script/recursive_feature_elimination/nnet_rfe.R')
# source('script/recursive_feature_elimination/knn_rfe.R')
# source('script/recursive_feature_elimination/mars_bagEarth_rfe.R')
# source('script/recursive_feature_elimination/neuralnet_rfe.R')

stopCluster(cl)

methods = c(
            'glmnet', 'lm', 'bag', 'kknn'
            , 'rf'
            # , 'blackboost', 'cubist' 
            # 'bstTree', 
            # 'nnet', 'svm'
            )
            # 'knn', 'bagEarth', 'neuralnet')
all_results = data.frame()

for (method_i in methods) {
  file_rds = paste0(path_output, '/', method_i, 'Profile.rds')
  prof = readRDS(file_rds)
  if (method_i == 'rf') {
    forceLibrary('randomForest')
    prot_pred = predict(prof$fit, X_test, type="response",
                        norm.votes=TRUE, predict.all=FALSE, proximity=FALSE, nodes=FALSE)
    
  } else {
    prot_pred = predict(prof$fit, X_test)
  }
  
  test_result = postResample(pred = prot_pred, obs = Y_test) 
  names(test_result) = paste0('test_', names(test_result))
  train_result = prof$results[prof$results$Variables == prof$bestSubset, 2:4] %>% unlist()
  names(train_result) = paste0('train_', names(train_result))
  
  all_result = c(train_result, test_result) %>% as.data.frame() %>% t()
  rownames(all_result) = method_i
  
    print(method_i)
  print(test_result)
  all_results = rbind.data.frame(all_results, all_result)
}
