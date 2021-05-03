library(lattice)
library(dplyr)

adaptRmseResults <- function(profile, reference) {
  rmse_list = reference$results$RMSE
  rmse_list[reference$results$Variables %in% profile$results$Variables] = profile$results$RMSE[profile$results$Variables %in% reference$results$Variables]
  rmse_list[!reference$results$Variables %in% profile$results$Variables] = NA
  return(rmse_list)
}

if (file.exists('../output_rfe/na_omit/rfProfile.rds')) {
  rfProfile = readRDS('../output_rfe/na_omit/rfProfile.rds')
} else {
  source('script/recursive_feature_elimination/random_forest_rfe.R')
}

if (file.exists('../output_rfe/na_omit/lmProfile.rds')) {
  lmProfile = readRDS('../output_rfe/na_omit/lmProfile.rds') %>% 
    adaptRmseResults(rfProfile)
} else {
  source('script/recursive_feature_elimination/linear_model_rfe.R')
}

# if (file.exists('../output_rfe/na_omit/lmProfile2.rds')) {
#   lmProfile2 = readRDS('../output_rfe/na_omit/lmProfile2.rds')
# } else {
#   source('script/recursive_feature_elimination/linear_model_rfe.R')
# }



if (file.exists('../output_rfe/na_omit/bagProfile.rds')) {
  bagProfile = readRDS('../output_rfe/na_omit/bagProfile.rds') %>% 
    adaptRmseResults(rfProfile)
  
} else {
  source('script/recursive_feature_elimination/bagged_cart_rfe.R')
}

# if (file.exists('../output_rfe/na_omit/svmProfile.rds')) {
#   svmProfile = readRDS('../output_rfe/na_omit/svmProfile.rds') %>% 
#     adaptRmseResults(rfProfile)
#   
# } else {
#   source('script/recursive_feature_elimination/support_vector_machine_rfe.R')
# }


blackboostProfile = readRDS('../output_rfe/na_omit/blackboostProfile.rds') %>% 
  adaptRmseResults(rfProfile)

bstTreeProfile = readRDS('../output_rfe/na_omit/bstTreeProfile.rds') %>% 
  adaptRmseResults(rfProfile)

glmnetProfile = readRDS('../output_rfe/na_omit/glmnetProfile.rds') %>% 
  adaptRmseResults(rfProfile)

nnetProfile = readRDS('../output_rfe/na_omit/nnetProfile.rds') %>% 
  adaptRmseResults(rfProfile)

cubistProfile = readRDS('../output_rfe/na_omit/cubistProfile.rds') %>% 
  adaptRmseResults(rfProfile)

kknnProfile = readRDS('../output_rfe/na_omit/kknnProfile.rds') %>% 
  adaptRmseResults(rfProfile)



png('output/xyplot_rfe_na_omit.png', width = 1024, height = 1024)
xyplot(lmProfile + 
         rfProfile$results$RMSE + 
         bagProfile +
         # svmProfile +
         bstTreeProfile +
         cubistProfile +
         glmnetProfile +
         kknnProfile +
         # nnetProfile +
         blackboostProfile ~
         rfProfile$results$Variables, xlab = 'Number of Variables', ylab = 'Root-mean-square error (RMSE)',
       type = c("g", "p", "l"),
       auto.key = TRUE)
dev.off()


varis = rfp %>% varImp() %>% head(7) %>% rownames
caret::featurePlot(x = X[, varis], y = Y)
