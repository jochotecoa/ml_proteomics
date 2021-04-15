source('script/recursive_feature_elimination/load_data_rfe.R')

forceLibrary('kernlab')

  svmProfile <- rfe(X, Y,
                    sizes = sizes,
                    rfeControl = rfeControl(functions = caretFuncs, verbose = T),
                    ## pass options to train()
                    method = "svmRadial")
  saveRDS(svmProfile, '../output/rfe/svmProfile.rds')
  stopCluster(cl)