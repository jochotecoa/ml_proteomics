library(lattice)
library(dplyr)

adaptRmseResults <- function(profile, reference) {
  rmse_list = reference$results$RMSE
  rmse_list[reference$results$Variables %in% profile$results$Variables] = profile$results$RMSE[profile$results$Variables %in% reference$results$Variables]
  rmse_list[!reference$results$Variables %in% profile$results$Variables] = NA
  return(rmse_list)
}

if (file.exists('../output_rfe/rfProfile.rds')) {
  rfProfile = readRDS('../output_rfe/rfProfile.rds')
} else {
  source('script/recursive_feature_elimination/random_forest_rfe.R')
}

if (file.exists('../output_rfe/lmProfile.rds')) {
  lmProfile = readRDS('../output_rfe/lmProfile.rds') %>% 
    adaptRmseResults(rfProfile)
} else {
  source('script/recursive_feature_elimination/linear_model_rfe.R')
}

# if (file.exists('../output_rfe/lmProfile2.rds')) {
#   lmProfile2 = readRDS('../output_rfe/lmProfile2.rds')
# } else {
#   source('script/recursive_feature_elimination/linear_model_rfe.R')
# }



if (file.exists('../output_rfe/bagProfile.rds')) {
  bagProfile = readRDS('../output_rfe/bagProfile.rds') %>% 
    adaptRmseResults(rfProfile)
  
} else {
  source('script/recursive_feature_elimination/bagged_cart_rfe.R')
}

if (file.exists('../output_rfe/svmProfile.rds')) {
  svmProfile = readRDS('../output_rfe/svmProfile.rds') %>% 
    adaptRmseResults(rfProfile)
  
} else {
  source('script/recursive_feature_elimination/support_vector_machine_rfe.R')
}


blackboostProfile = readRDS('../output_rfe/blackboostProfile.rds') %>% 
  adaptRmseResults(rfProfile)

bstTreeProfile = readRDS('../output_rfe/bstTreeProfile.rds') %>% 
  adaptRmseResults(rfProfile)

glmnetProfile = readRDS('../output_rfe/glmnetProfile.rds') %>% 
  adaptRmseResults(rfProfile)

nnetProfile = readRDS('../output_rfe/nnetProfile.rds') %>% 
  adaptRmseResults(rfProfile)

cubistProfile = readRDS('../output_rfe/cubistProfile.rds') %>% 
  adaptRmseResults(rfProfile)



png('output/xyplot_rfe.png', width = 1024, height = 1024)
xyplot(lmProfile + 
         rfProfile$results$RMSE + 
         bagProfile +
         svmProfile +
         bstTreeProfile +
         cubistProfile +
         glmnetProfile +
         # nnetProfile +
         blackboostProfile ~
         rfProfile$results$Variables, xlab = 'Number of Variables', ylab = 'Root-mean-square error (RMSE)',
       type = c("g", "p", "l"),
       auto.key = TRUE)
dev.off()
