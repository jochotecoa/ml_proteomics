# kknn
# source('script/recursive_feature_elimination/load_data_rfe.R')
# KNN
forceLibrary('kknn')
sizes = as.integer(seq(2, ncol(X), (ncol(X)-2)/4))[-5]

kknnProfile <- rfe(X, Y,
                  sizes = sizes,
                  rfeControl = rfeControl(functions = caretFuncs, 
                                          verbose = T, 
                                          method = 'repeatedcv'),
                  ## pass options to train()
                  method = "kknn")
file_rds = paste0(path_output, '/kknnProfile.rds')
saveRDS(kknnProfile, file_rds)
stopCluster(cl)
