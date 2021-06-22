# source('script/recursive_feature_elimination/load_data_rfe.R')
forceLibrary('randomForest')

Xi = X[,-grep('seq_depth_mirna', colnames(X)), F]

stopifnot(ncol(Xi)<ncol(X))

rfProfile <- rfe(Xi, Y,
                 sizes = sizes,
                 rfeControl = rfeControl(functions = rfFuncs,
                                         method = 'repeatedcv', 
                                         repeats = repeats,
                                         verbose = T,
                                         index = folds,
                                         allowParallel = T))

file_rds = paste0(path_output, '/rfProfile_repeatedcv_wout_seq_depth_mirna.rds')
saveRDS(rfProfile, file_rds)
# stopCluster(cl)
