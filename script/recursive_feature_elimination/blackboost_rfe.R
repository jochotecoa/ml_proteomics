# source('script/recursive_feature_elimination/load_data_rfe.R')

forceLibrary(c('party', 'mboost', 'plyr', 'partykit'))
# blackboost

sizes = 1:ncol(X)
# sizes = as.integer(seq(2, ncol(X), (ncol(X)-2)/4))[-5]

blackboostProfile <- rfe(X, Y,
                     sizes = sizes,
                     rfeControl = rfeControl(functions = caretFuncs, 
                                             verbose = T, 
                                             method = 'repeatedcv', 
                                             index = folds),
                     ## pass options to train()
                     method = "blackboost")
file_rds = paste0(path_output, '/blackboostProfile_wout_TPM_value_mean_log2.rds')
saveRDS(blackboostProfile, file_rds)

# stopCluster(cl)
