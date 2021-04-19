source('script/recursive_feature_elimination/load_data_rfe.R')

forceLibrary('kernlab')

sizes = as.integer(seq(2, ncol(X), (ncol(X)-2)/4))[-5]

  svmProfile <- rfe(X, Y,
                    sizes = sizes,
                    rfeControl = rfeControl(functions = caretFuncs, 
                                            verbose = T, 
                                            method = 'repeatedcv'),
                    ## pass options to train()
                    method = "svmRadial")
  saveRDS(svmProfile, '../output/rfe/svmProfile.rds')
  stopCluster(cl)
  