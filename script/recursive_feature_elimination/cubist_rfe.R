# source('script/recursive_feature_elimination/load_data_rfe.R')
# cubist
forceLibrary('Cubist')

sizes = as.integer(seq(2, ncol(X), (ncol(X)-2)/4))[-5]

cubistProfile <- rfe(X, Y,
                   sizes = sizes,
                   rfeControl = rfeControl(functions = caretFuncs, 
                                           verbose = T, 
                                           method = 'repeatedcv'),
                   ## pass options to train()
                   method = "cubist")
file_rds = paste0(path_output, '/cubistProfile.rds')
saveRDS(cubistProfile, file_rds)
# stopCluster(cl)
