# source('script/recursive_feature_elimination/load_data_rfe.R')
forceLibrary('randomForest')

rfProfile <- rfe(X, Y,
                 sizes = sizes,
                 rfeControl = rfeControl(functions = rfFuncs,
                                         method = 'cv', 
                                         verbose = T,
                                         index = folds,
                                         allowParallel = T))

file_rds = paste0(path_output, '/rfProfile.rds')
saveRDS(rfProfile, file_rds)
# stopCluster(cl)
