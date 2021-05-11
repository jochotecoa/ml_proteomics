# glmnet
# source('script/recursive_feature_elimination/load_data_rfe.R')

forceLibrary('glmnet', 'Matrix')

sizes = as.integer(seq(2, ncol(X), (ncol(X)-2)/4))[-5]

glmnetProfile <- rfe(X, Y,
                  sizes = sizes,
                  rfeControl = rfeControl(functions = caretFuncs, 
                                          verbose = T, 
                                          method = 'repeatedcv', 
                                          index = folds),
                  ## pass options to train()
                  method = "glmnet")
file_rds = paste0(path_output, '/glmnetProfile.rds')
saveRDS(glmnetProfile, file_rds)
# stopCluster(cl)
