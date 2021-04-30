source('script/recursive_feature_elimination/load_data_rfe.R')

forceLibrary(c('party', 'mboost', 'plyr', 'partykit'))
# blackboost

sizes = as.integer(seq(2, ncol(X), (ncol(X)-2)/4))[-5]

blackboostProfile <- rfe(X, Y,
                     sizes = sizes,
                     rfeControl = rfeControl(functions = caretFuncs, 
                                             verbose = T, 
                                             method = 'repeatedcv'),
                     ## pass options to train()
                     method = "blackboost")
file_rds = paste0(path_output, '/blackboostProfile.rds')
saveRDS(blackboostProfile, file_rds)

stopCluster(cl)
