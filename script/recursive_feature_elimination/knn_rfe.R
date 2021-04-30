# source('script/recursive_feature_elimination/load_data_rfe.R')
# KNN

sizes = as.integer(seq(2, ncol(X), (ncol(X)-2)/4))[-5]

knnProfile <- rfe(X, Y,
                     sizes = sizes,
                     rfeControl = rfeControl(functions = caretFuncs, 
                                             verbose = T, 
                                             method = 'repeatedcv'),
                     ## pass options to train()
                       method = "knn")
file_rds = paste0(path_output, '/knnProfile.rds')
saveRDS(knnProfile, file_rds)
stopCluster(cl)
