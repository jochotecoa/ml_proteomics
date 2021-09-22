# source('script/recursive_feature_elimination/load_data_rfe.R')

forceLibrary(c('ipred', 'plyr', 'e1071')) # Needed for bagged trees

  bagProfile <- rfe(X, Y,
                    sizes = sizes,
                    rfeControl = rfeControl(functions = treebagFuncs,
                                            verbose = T,
                                            allowParallel = T, 
                                            method = 'cv', 
                                            index = folds))
  file_rds = paste0(path_output, '/bagProfile.rds')
  saveRDS(bagProfile, file_rds)
  # stopCluster(cl)