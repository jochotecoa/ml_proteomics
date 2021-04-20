library(lattice)

if (file.exists('../output_rfe/lmProfile.rds')) {
  lmProfile = readRDS('../output_rfe/lmProfile.rds')
} else {
  source('script/recursive_feature_elimination/linear_model_rfe.R')
}

if (file.exists('../output_rfe/lmProfile2.rds')) {
  lmProfile2 = readRDS('../output_rfe/lmProfile2.rds')
} else {
  source('script/recursive_feature_elimination/linear_model_rfe.R')
}

if (file.exists('../output_rfe/rfProfile.rds')) {
  rfProfile = readRDS('../output_rfe/rfProfile.rds')
} else {
  source('script/recursive_feature_elimination/random_forest_rfe.R')
}


if (file.exists('../output_rfe/bagProfile.rds')) {
  bagProfile = readRDS('../output_rfe/bagProfile.rds')
} else {
  source('script/recursive_feature_elimination/bagged_cart_rfe.R')
}

if (file.exists('../output_rfe/svmProfile.rds')) {
  svmProfile = readRDS('../output_rfe/svmProfile.rds')
} else {
  source('script/recursive_feature_elimination/support_vector_machine_rfe.R')
}


png('output/xyplot_rfe.png')
xyplot(lmProfile$results$RMSE[rfProfile$results$Variables] + 
         # lmProfile2$results$RMSE + 
         rfProfile$results$RMSE + 
         bagProfile$results$RMSE ~
         # svmProfile$results$RMSE ~
         rfProfile$results$Variables, xlab = 'Number of Variables', ylab = 'Root-mean-square deviation (RMSE)',
       type = c("g", "p", "l"),
       auto.key = TRUE)
dev.off()
