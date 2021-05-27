library(lattice)
library(dplyr)

adaptRmseResults <- function(profile, reference, metric = 'RMSE') {
  rmse_list = unlist(reference$results[metric])
  rmse_list[reference$results$Variables %in% profile$results$Variables] = unlist(profile$results[metric])[profile$results$Variables %in% reference$results$Variables]
  rmse_list[!reference$results$Variables %in% profile$results$Variables] = NA
  names(rmse_list) = NULL
  rmse_list %>% return()
}

  
metric = 'Rsquared' # Rsquared


glmnetProfile = readRDS('../output_rfe/split_by_sample/na_omit/glmnetProfile.rds') 

Generalized_linear_model = glmnetProfile %>% 
  adaptRmseResults(glmnetProfile, metric = metric)

Linear_model = readRDS('../output_rfe/split_by_sample/na_omit/lmProfile.rds') %>% 
  adaptRmseResults(glmnetProfile, metric = metric)

Bagged_trees = readRDS('../output_rfe/split_by_sample/na_omit/bagProfile.rds') %>% 
  adaptRmseResults(glmnetProfile, metric = metric)

K_nearest_neighbours = readRDS('../output_rfe/split_by_sample/na_omit/kknnProfile.rds') %>% 
  adaptRmseResults(glmnetProfile, metric = metric)
    
Random_forest = readRDS('../output_rfe/split_by_sample/na_omit/rfProfile.rds') %>% 
  adaptRmseResults(glmnetProfile, metric = metric)
  
Blackboost = readRDS('../output_rfe/split_by_prot/na_omit/blackboostProfile.rds') %>% 
  adaptRmseResults(glmnetProfile, metric = metric)

Cubist = readRDS('../output_rfe/split_by_sample/na_omit/cubistProfile.rds') %>% 
  adaptRmseResults(glmnetProfile, metric = metric)

Support_vector_machines = readRDS('../output_rfe/split_by_sample/na_omit/svmProfile.rds') %>% 
  adaptRmseResults(glmnetProfile, metric = metric)

Boosted_trees = readRDS('../output_rfe/split_by_sample/na_omit/bstTreeProfile.rds') %>% 
  adaptRmseResults(glmnetProfile, metric = metric)

Neural_networks = readRDS('../output_rfe/split_by_sample/na_omit/nnetProfile.rds') %>% 
  adaptRmseResults(glmnetProfile, metric = metric)





png('output/xyplot_rfe_na_omit.png', width = 1024, height = 1024)
xyplot(Generalized_linear_model + 
         Linear_model +
         Bagged_trees +
         K_nearest_neighbours +
         Random_forest +
         Blackboost +
         Cubist +
         # nnetProfile +
         Support_vector_machines +
         Boosted_trees ~
         glmnetProfile$results$Variables, xlab = 'Number of Variables', ylab = metric,
       type = c("g", "p", "l"),
       auto.key = TRUE)
dev.off()


varis = rfp %>% varImp() %>% head(7) %>% rownames
caret::featurePlot(x = X[, varis], y = Y)
