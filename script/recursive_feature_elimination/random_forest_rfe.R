# source('script/recursive_feature_elimination/load_data_rfe.R')
forceLibrary('randomForest')

# sizes = as.integer(seq(2, ncol(X), (ncol(X)-2)/4))[-5]
sizes = as.integer(seq(2, 35, 1))

rfProfile <- rfe(X, Y,
                 sizes = sizes,
                 rfeControl = rfeControl(functions = rfFuncs,
                                         verbose = T,
                                         allowParallel = T,
                                         method = 'repeatedcv', 
                                         index = folds))

file_rds = paste0(path_output, '/rfProfile.rds')
saveRDS(rfProfile, file_rds)
# stopCluster(cl)
