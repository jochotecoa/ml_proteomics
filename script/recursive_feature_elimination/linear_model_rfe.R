# source('script/recursive_feature_elimination/load_data_rfe.R')

generate_new_sizes <- function(bestSubset, ncol_X, sizes) {
  if(length(sizes) == 1){return(bestSubset)}
  besub_n = bestSubset
  pos_besub = grep(besub_n, sizes)
  if (besub_n == ncol_X) {
    if(length(pos_besub) == 0){
      sizes = c(sizes, ncol_X)
      pos_besub = length(sizes)
    }
  }
  mini = 1
  # mini = sizes[pos_besub - 1]
  # if (length(mini) == 0) {
  #   mini = min(sizes)
  # }
  maxi = besub_n
  # maxi = sizes[pos_besub + 1]
  # if (is.na(maxi)) {
  #   if (besub_n == ncol_X) {
  #     maxi = ncol_X
  #   } else {
  #     maxi = max(sizes)
  #   }
  # }
  sizes = seq(mini, maxi, (maxi - mini)/5) %>% round() %>% unique()
  if (maxi == ncol_X) {
    sizes = sizes[-6]
    }
  return(sizes)
  
}

  # b <- rfe(x = mtcars[,-1],
  #                  y = mtcars[,1],
  #                  sizes = 1:ncol(mtcars),
  #                  rfeControl = rfeControl(functions = lmFuncs, verbose = T))
  
# besub = 0
# besubs = NULL
# optimum = F
# sizes = as.integer(c(seq(2, .36*ncol(X), 1), seq((.36*ncol(X))+5, ncol(X)-5, 5)))
sizes = seq(1, ncol(X)) %>% as.integer()
# while (optimum == F) {
  lmProfile <- rfe(X, Y,
                   sizes = sizes,
                   rfeControl = rfeControl(functions = lmFuncs,
                                           verbose = T,
                                           method = 'repeatedcv', 
                                           allowParallel = T, 
                                           index = folds))
#   besub_n = lmProfile$bestSubset
#   if (besub == besub_n) {
#     optimum = T
#   } else {
#     besub = besub_n
#     besubs = c(besubs, besub)
#   }
# 
#   sizes = generate_new_sizes(bestSubset = besub_n,
#                              ncol_X = ncol(X),
#                              sizes = sizes)
# 
#   print(besubs)
#   print(sizes)
# 
# }

  file_rds = paste0(path_output, '/lmProfile.rds')
  saveRDS(lmProfile, file_rds)

# With rerank -------------------------------------------------------------


  # source('script/recursive_feature_elimination/load_data_rfe.R')
  
  # besub = 0
  # besubs = NULL
  # optimum = F
  # 
  # while (optimum == F) {
    
    # lmProfile2 <- rfe(X, Y,
    #                   sizes = sizes,
    #                   rfeControl = rfeControl(functions = lmFuncs,
    #                                           rerank = TRUE, 
    #                                           method = 'repeatedcv',
    #                                           verbose = F))
  #   besub_n = lmProfile2$bestSubset
  #   if (besub == besub_n) {
  #     optimum = T
  #   } else {
  #     besub = besub_n
  #     besubs = c(besubs, besub)
  #   }
  #   
  #   sizes = generate_new_sizes(bestSubset = besub_n, 
  #                              ncol_X = ncol(X), 
  #                              sizes = sizes)
  #   
  #   print(besubs)
  #   print(sizes)
  # }
  
  # saveRDS(lmProfile2, '../output_rfe/lmProfile2.rds')
  # stopCluster(cl)
  