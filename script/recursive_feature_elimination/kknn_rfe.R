# kknn
source('script/recursive_feature_elimination/load_data_rfe.R')
# KNN
forceLibrary('kknn')
sizes = as.integer(seq(2, ncol(X), (ncol(X)-2)/4))[-5]

knnProfile <- rfe(X, Y,
                  sizes = sizes,
                  rfeControl = rfeControl(functions = caretFuncs, 
                                          verbose = T, 
                                          method = 'repeatedcv'),
                  ## pass options to train()
                  method = "knn")
saveRDS(knnProfile, '../output_rfe/knnProfile.rds')
stopCluster(cl)
