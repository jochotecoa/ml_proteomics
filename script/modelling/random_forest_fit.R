# source('script/recursive_feature_elimination/load_data_rfe.R')
forceLibrary('randomForest')

sizes = 1:ncol(X)
# sizes = as.integer(seq(2, ncol(X), (ncol(X)-2)/4))[-5]

RMSE_df = data.frame()

X_clean_2 = X %>% 
  dplyr::select(circ_sum, circ_mean, circ_mean_log2, linear_density, TPM_value_mean_log2, mirna_mean_stringent, mirna_mean_stringent_log2, mirna_sum_stringent, Aa_M_prop, Aa_D_prop)

  rfFit <- train(X_clean_2, Y,
                 method = 'rf',
                 trControl  = trainControl(method = 'cv', 
                                           index = folds))

dir.create('output/cardiac/na_omit/split_by_sample/72_features/')
saveRDS(object = rfFit, file = 'output/cardiac/na_omit/split_by_sample/72_features/rfFit.rds')

