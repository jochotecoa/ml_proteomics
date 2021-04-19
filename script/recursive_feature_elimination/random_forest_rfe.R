source('script/recursive_feature_elimination/load_data_rfe.R')
forceLibrary('randomForest')

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

sizes = as.integer(seq(2, ncol(X), (ncol(X)-2)/4))[-5]

# besub = 0
# besubs = NULL
# optimum = F
# 
# while (optimum == F) {
#   
# 
rfProfile <- rfe(X, Y,
                 sizes = sizes,
                 rfeControl = rfeControl(functions = rfFuncs,
                                         verbose = F,
                                         allowParallel = T,
                                         method = 'repeatedcv'))
# besub_n = rfProfile$bestSubset
# if (besub == besub_n) {
#   optimum = T
# } else {
#   besub = besub_n
#   besubs = c(besubs, besub)
# }
# 
# sizes = generate_new_sizes(bestSubset = besub_n,
#                            ncol_X = ncol(X),
#                            sizes = sizes)
# 
# print(besubs)
# print(sizes)
# 
# }

saveRDS(rfProfile, '../output_rfe/rfProfile.rds')
stopCluster(cl)
