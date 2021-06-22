# source('script/recursive_feature_elimination/load_data_rfe.R')
forceLibrary('randomForest')

sizes = 1:ncol(X)
# sizes = as.integer(seq(2, ncol(X), (ncol(X)-2)/4))[-5]

RMSE_df = data.frame()

while (ncol(X) >= 1) {
  
  rfFit <- train(X, Y,
                 method = 'rf',
                 trControl  = trainControl(method = 'cv', 
                                           index = folds))
  
  RMSE_num = c(ncol(X), min(rfFit$results$RMSE))
  
  RMSE_df = rbind.data.frame(RMSE_df, RMSE_num)
  colnames(RMSE_df) = c('n_vars', 'RMSE')
  
  saveRDS(object = RMSE_df, file = 'output/rfFit_naomit_samplesplit_RMSEs_df.rds')
  saveRDS(object = rfFit, file = 'output/rfFit_latest_naomit_samplesplit.rds')
  
  X = X[, -ncol(X), F]
}


file_rds = paste0(path_output, '/rfProfile_wout_seq_depth_mirna.rds')
saveRDS(rfProfile, file_rds)
stopCluster(cl)

# version 2 ---------------------------------------------------------------


forceLibrary('randomForest')

RMSE_df = data.frame()

for (i in ncol(X):1) {
  Xi = X[, -i, F]
  
  rfFit <- train(Xi, Y,
                 method = 'rf',
                 trControl  = trainControl(method = 'cv', 
                                           index = folds))
  
  RMSE_num = c(i, min(rfFit$results$RMSE))
  
  RMSE_df = rbind.data.frame(RMSE_df, RMSE_num)
  colnames(RMSE_df) = c('n_vars', 'RMSE')
  
  saveRDS(object = RMSE_df, file = 'output/RMSEs_rfFit_naomit_samplesplit_v2.1_df.rds')
  saveRDS(object = rfFit, file = 'output/rfFit_latest_naomit_samplesplit_v2.1.rds')
  
}


# version 3 ---------------------------------------------------------------

forceLibrary('randomForest')

RMSE_df = data.frame()

for (i in 1:ncol(X)) {
  Xi = X [, c(i, grep('Version', colnames(X))), F]
  
  
  rfFit <- train(Xi, Y,
                 method = 'rf',
                 trControl  = trainControl(method = 'cv', 
                                           index = folds))
  
  RMSE_num = c(i, min(rfFit$results$RMSE))
  
  RMSE_df = rbind.data.frame(RMSE_df, RMSE_num)
  colnames(RMSE_df) = c('n_vars', 'RMSE')
  
  saveRDS(object = RMSE_df, file = 'output/RMSEs_rfFit_naomit_samplesplit_v3_df.rds')
  saveRDS(object = rfFit, file = 'output/rfFit_latest_naomit_samplesplit_v3.rds')
  
}
stopCluster(cl)

