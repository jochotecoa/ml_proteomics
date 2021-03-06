# source('script/recursive_feature_elimination/load_data_rfe.R')

forceLibrary(c('bst', 'plyr'))
# bstTree

sizes = as.integer(seq(2, ncol(X), (ncol(X)-2)/4))[-5]

bstTreeProfile <- rfe(X, Y,
                         sizes = sizes,
                         rfeControl = rfeControl(functions = caretFuncs, 
                                                 verbose = T, 
                                                 method = 'repeatedcv', 
                                                 index = folds),
                         ## pass options to train()
                         method = "bstTree")
file_rds = paste0(path_output, '/bstTreeProfile.rds')
saveRDS(bstTreeProfile, file_rds)

# stopCluster(cl)
