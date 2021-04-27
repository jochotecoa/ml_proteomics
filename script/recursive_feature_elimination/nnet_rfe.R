source('script/recursive_feature_elimination/load_data_rfe.R')
# neural network
forceLibrary('nnet')
sizes = as.integer(seq(2, ncol(X), (ncol(X)-2)/4))[-5]

nnetProfile <- rfe(X, Y,
                        sizes = sizes,
                        rfeControl = rfeControl(functions = caretFuncs, 
                                                verbose = T, 
                                                method = 'repeatedcv'),
                        ## pass options to train()
                        method = "nnet")
saveRDS(nnetProfile, '../output_rfe/nnetProfile.rds')
stopCluster(cl)
