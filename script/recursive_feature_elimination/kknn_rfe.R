# kknn
# source('script/recursive_feature_elimination/load_data_rfe.R')
# KNN
forceLibrary('kknn')

kknnProfile <- rfe(X, Y,
                  sizes = sizes,
                  rfeControl = rfeControl(functions = caretFuncs, 
                                          verbose = T, 
                                          method = 'cv', 
                                          index = folds),
                  ## pass options to train()
                  method = "kknn")
file_rds = paste0(path_output, '/kknnProfile.rds')
saveRDS(kknnProfile, file_rds)
# stopCluster(cl)
