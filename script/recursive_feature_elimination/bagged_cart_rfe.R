source('script/recursive_feature_elimination/load_data_rfe.R')

forceLibrary(c('ipred', 'plyr', 'e1071')) # Needed for bagged trees

  
  bagProfile <- rfe(X, Y,
                    sizes = sizes,
                    rfeControl = rfeControl(functions = treebagFuncs,
                                            verbose = T,
                                            allowParallel = T, number = 1))
  saveRDS(bagProfile, 'output/rfe/bagProfile.rds')
