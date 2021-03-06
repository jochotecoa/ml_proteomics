
sizes = seq(1, ncol(X)) %>% as.integer()
# while (optimum == F) {
lmProfile <- rfe(X, Y,
                 sizes = sizes,
                 rfeControl = rfeControl(functions = lmFuncs,
                                         verbose = T,
                                         method = 'repeatedcv', 
                                         allowParallel = T, 
                                         index = folds))

train_samples = X_train %>% rownames() %>% strsplit('--') %>% sapply('[[', 2)

folds = groupKFold(train_samples, k = 10)


lmProfile_sample <- rfe(X, Y,
                        sizes = sizes,
                        rfeControl = rfeControl(functions = lmFuncs,
                                                verbose = T,
                                                method = 'repeatedcv', 
                                                allowParallel = T, 
                                                index = folds))




old_RMSE = 99
new_RMSE = 98



while (old_RMSE > new_RMSE) {
  old_RMSE = new_RMSE
  
  if (old_RMSE != 98) {
    bad_feats = lst$var_rm[lst$diff > 0] %>% unlist %>% .[-1]
    
    X_train_2 = X_train_2 %>% 
      dplyr::select(!(!!bad_feats))
    
    X_test_2 = X_test_2 %>% 
      dplyr::select(!(!!bad_feats))
    
  } else {
    X_train_2 = X_train
    X_test_2 = X_test
  }
  
  
  sizes = seq(2, ncol(X_train_2)) %>% as.integer()
  # while (optimum == F) {
  lmProfile_2 <- rfeIter(x = X_train_2, y = Y, testX = X_test_2, testY = Y_test,
                         sizes = sizes,
                         rfeControl = rfeControl(functions = lmFuncs,
                                                 verbose = T,
                                                 method = 'repeatedcv', 
                                                 allowParallel = T))
  
  lst = data.frame()
  for (n_vars in 2:max(lmProfile_2$pred$Variables)) {
    
    mumr = postResample(pred = lmProfile_2$pred$pred[lmProfile_2$pred$Variables == n_vars], obs = lmProfile_2$pred$obs[lmProfile_2$pred$Variables == n_vars]) %>% as.data.frame() %>% t()
    rownames(mumr) = n_vars
    
    lst = rbind.data.frame(lst, mumr)
    
  }
  
  new_RMSE = min(lst$RMSE)
  
  if (old_RMSE > new_RMSE) {
    oldProf = lmProfile_2
  } else {
    newProf = lmProfile_2
  }
  
  lst$diff = NA
  lst$diff[2:(max(lmProfile_2$pred$Variables) - 1)] = lst$RMSE %>% diff
  
  lst$var_rm = NA
  
  for (nvars in 2:(max(lmProfile_2$pred$Variables) - 1)) {
    lmp_i = (max(lmProfile_2$pred$Variables) + 1) - nvars 
    lst$var_rm[nvars] = setdiff(rownames(lmProfile_2$finalVariables[[lmp_i-1]]), rownames(lmProfile_2$finalVariables[[lmp_i]]))
  }
  
  
}


sizes = seq(2, ncol(X)) %>% as.integer()
# while (optimum == F) {
lmProfile <- rfeIter(x = X, y = Y, testX = X_test, testY = Y_test,
                     sizes = sizes,
                     rfeControl = rfeControl(functions = lmFuncs,
                                             verbose = T,
                                             method = 'repeatedcv', 
                                             allowParallel = T))

lst = data.frame()
for (n_vars in 2:42) {
  
  mumr = postResample(pred = lmProfile$pred$pred[lmProfile$pred$Variables == n_vars], obs = lmProfile$pred$obs[lmProfile$pred$Variables == n_vars]) %>% as.data.frame() %>% t()
  rownames(mumr) = n_vars
  
  lst = rbind.data.frame(lst, mumr)
  
}

lst$diff = NA
lst$diff[2:41] = lst$RMSE %>% diff

lst$var_rm = NA

for (nvars in 2:41) {
  lmp_i = 43 - nvars 
  lst$var_rm[nvars] = setdiff(rownames(lmProfile$finalVariables[[lmp_i-1]]), rownames(lmProfile$finalVariables[[lmp_i]]))
}

bad_feats = lst$var_rm[lst$diff > 0] %>% unlist %>% .[-1]

X_train_2 = X_train %>% 
  dplyr::select(!(!!bad_feats))

X_test_2 = X_test %>% 
  dplyr::select(!(!!bad_feats))

sizes = seq(2, ncol(X_train_2)) %>% as.integer()
# while (optimum == F) {
lmProfile_2 <- rfeIter(x = X_train_2, y = Y, testX = X_test_2, testY = Y_test,
                       sizes = sizes,
                       rfeControl = rfeControl(functions = lmFuncs,
                                               verbose = T,
                                               method = 'repeatedcv', 
                                               allowParallel = T))

lst_2 = data.frame()
for (n_vars in 2:ncol(X_train_2)) {
  
  mumr = postResample(pred = lmProfile_2$pred$pred[lmProfile_2$pred$Variables == n_vars], obs = lmProfile_2$pred$obs[lmProfile_2$pred$Variables == n_vars]) %>% as.data.frame() %>% t()
  rownames(mumr) = n_vars
  
  lst_2 = rbind.data.frame(lst_2, mumr)
  
}

