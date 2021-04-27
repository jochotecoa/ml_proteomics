source('script/recursive_feature_elimination/load_data_rfe.R')
# MARS

forceLibrary('earth')

sizes = as.integer(seq(2, ncol(X), (ncol(X)-2)/4))[-5]

bagEarthProfile <- rfe(X, Y,
                     sizes = sizes,
                     rfeControl = rfeControl(functions = caretFuncs, 
                                             verbose = T, 
                                             method = 'repeatedcv'),
                     ## pass options to train()
                     method = "bagEarth")
saveRDS(bagEarthProfile, '../output_rfe/bagEarthProfile.rds')
stopCluster(cl)
