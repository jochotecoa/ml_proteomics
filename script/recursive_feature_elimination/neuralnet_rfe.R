source('script/recursive_feature_elimination/load_data_rfe.R')
# neural network
forceLibrary('neuralnet')
sizes = as.integer(seq(2, ncol(X), (ncol(X)-2)/4))[-5]

neuralnetProfile <- rfe(X, Y,
                  sizes = sizes,
                  rfeControl = rfeControl(functions = caretFuncs, 
                                          verbose = T, 
                                          method = 'repeatedcv'),
                  ## pass options to train()
                  method = "neuralnet")
file_rds = paste0(path_output, '/neuralnetProfile.rds')
saveRDS(neuralnetProfile, file_rds)
stopCluster(cl)
