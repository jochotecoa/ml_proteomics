source('script/functions/functions_JOA.R')
forceLibrary(c('dplyr', 'tibble', 'biomaRt', 'reshape', 'caret', 'caTools'))

tissue = 'hepatic'
compound = 't0_controls_ML'

source('script/data_cleaning/mrna_data_cleaning.R')

source('script/data_cleaning/protein_data_cleaning.R')

source('script/data_cleaning/transcript_expression_per_protein.R')

source('script/data_cleaning/testing_sample_compability_between_transcrx_and_protx.R')

# if (file.exists('data/whole_raw_dataset.rds')) {
#   mrna_prot_df = readRDS('data/whole_raw_dataset.rds')
# } else {
  source('script/data_cleaning/combine_transcrx_with_protx.R')
# }

source('script/data_cleaning/combine_mirna_with_protx.R')

source('script/data_cleaning/combine_circrna_with_protx.R')

source('script/data_cleaning/combine_prot_stability_with_protx.R')

# Pre-processing ----------------------------------------------------------
# Zero- and Near Zero-Variance Predictors

nzv_m <- nearZeroVar(mrna_prot_df, saveMetrics = T)
nzv <- nearZeroVar(mrna_prot_df)
nzv_m[nzv, ] %>% View()

if (length(nzv) > 0) {
  zv_colnames = colnames(mrna_prot_df)[nzv_m$zeroVar]
  nzv_colnames = colnames(mrna_prot_df)[as.logical(nzv_m$nzv - nzv_m$zeroVar)]
  nzv_m[nzv, ] %>% print()
  # warning(paste(paste0(zv_colnames, collapse = ', '), 'presented zero variance'))
  # warning(paste(paste0(nzv_colnames, collapse = ', '), 'presented near zero variance'))
  mrna_prot_df <- mrna_prot_df[, -nzv]
}

# Creating Dummy Variables

dummies <- dummyVars(" ~ .", data = mrna_prot_df)
mrna_prot_df = predict(dummies, newdata = mrna_prot_df)
mrna_prot_df = mrna_prot_df %>% as.data.frame()


# Handle Missing Values ---------------------------------------------------

# feat_missval = mrna_prot_df %>% apply(2, anyNA) %>% subset(., .) %>% names()
# miss_per_feature = mrna_prot_df[, feat_missval] %>% summary %>% .[7, ] %>% 
#   gsub(pattern = "NA's   :| ", replacement = "") %>% unlist() %>% as.numeric()
# if (any(miss_per_feature >= (0.05*nrow(mrna_prot_df)))) {
#   stop()
# } else {
#   prev_rows = nrow(mrna_prot_df)
#   mrna_prot_df = mrna_prot_df %>% na.omit()
#   warning(paste(prev_rows - (nrow(mrna_prot_df)), 'rows have been removed due to missing values'))
# }

# mrna_prot_df = readRDS('data/training_data.rds')
# 
X = mrna_prot_df[, -grep('proteomics', colnames(mrna_prot_df))] %>% as.data.frame()
Y = mrna_prot_df[, grep('proteomics', colnames(mrna_prot_df))]


# Identifying Correlated Predictors

X = X %>% 
  dplyr::select(-circ_score_min) # it only has 2 values: 0 (rows with NAs) and 8


descrCor <- cor(na.omit(X))
highlyCor <- findCorrelation(descrCor, cutoff = .75)
corrplot::corrplot(descrCor[-highlyCor, highlyCor])
if (length(highlyCor) > 0) {
  X <- X[,-highlyCor]
}


# Linear Dependencies

comboInfo <- findLinearCombos(na.omit(X))
if (length(comboInfo$remove) > 0) {
  X = X[, -comboInfo$remove]
}

X = X %>% as.data.frame()

normalization <- preProcess(X, verbose = T, method = c("center", "scale", "knnImpute"))
X <- predict(normalization, X) %>% 
  as.data.frame() 


saveRDS(object = X, file = 'data/training_data_preds.rds')
saveRDS(object = Y, file = 'data/training_data_target.rds')

# 
# 
# trainIndex = createDataPartition(mrna_prot_df$proteomics_value, p=.8, list=F, 
#                                  times = 1)
# 
# protTrain <- mrna_prot_df[ trainIndex,]
# protTest  <- mrna_prot_df[-trainIndex,]
# 
# # Model building ----------------------------------------------------------
# # [1] "Linear Regression (lm)"                                             "eXtreme Gradient Boosting (xgbDART)"                               
# # [3] "Quantile Random Forest (qrf)"                                       "Relevance Vector Machines with Linear Kernel (rvmLinear)"          
# # [5] "Partial Least Squares (pls)"                                        "Bayesian Regularized Neural Networks (brnn)"                       
# # [7] "glmnet (glmnet)"                                                    "Adaptive-Network-Based Fuzzy Inference System (ANFIS)"             
# # [9] "Gaussian Process with Radial Basis Function Kernel (gaussprRadial)" "Linear Regression with Stepwise Selection (lmStepAIC)"             
# # [11] "Robust Linear Model (rlm)"                                         
# 
# 
# library(caret)# Simple linear regression model (lm means linear model)
# 
# 
# fitControl <- trainControl(method = "repeatedcv",   
#                            number = 10,     # number of folds
#                            repeats = 3, 
#                            search = "random", 
#                            allowParallel= T)    # repeated ten times
# 
# methods_caret = c('lm', 'xgbDART', 'qrf', 'rvmLinear', 'pls', 'brnn', 'glmnet', 'ANFIS', 'gaussprRadial', 'lmStepAIC', 'rlm')
# 
# models_list = list()
# 
# for (method_caret in methods_caret) {
#   model <- train(proteomics_value ~ .,
#                  data = mrna_prot_df,
#                  method = method_caret, 
#                  trControl = fitControl, 
#                  preProcess = c('scale', 'center'), 
#                  tuneLength = 10)
#   
#   models_list$new = model
#   names(models_list)[length(models_list)] = method_caret
#   
# }
# 
# model_lm <- train(proteomics_value ~ log2_TPM_value,
#                data = mrna_prot_df,
#                method = "lm", 
#                trControl = fitControl, 
#                preProcess = c('scale', 'center'))
# 
# model_mlm <- train(proteomics_value ~ .,
#                    data = mrna_prot_df,
#                    method = "lm",
#                    trControl = fitControl, preProcess = c('scale', 'center'))
# 
# 
# # possible values: boot", "boot632", "cv", "repeatedcv", "LOOCV", "LGOCV"
# 
# model_lasso <- train(proteomics_value ~ .,
#                   data = mrna_prot_df,
#                   method = "lasso",  # now we're using the lasso method
#                   trControl = fitControl, preProcess = c('scale', 'center'))
# 
# 
# # Multiple linear regression model
# 
# 
# # Ridge regression model
# model_ridge <- train(proteomics_value ~ .,
#                data = mrna_prot_df,
#                method = "ridge", # Try using "lasso"
#                trControl = fitControl, preProcess = c('scale', 'center'))
# 
# 
# resamps <- resamples(list(LM = model_lm,
#                           MLM = model_mlm,
#                           LSS = model_lasso,
#                           RDG = model_ridge))
# 
# theme1 <- trellis.par.get()
# theme1$plot.symbol$col = rgb(.2, .2, .2, .4)
# theme1$plot.symbol$pch = 16
# theme1$plot.line$col = rgb(1, 0, 0, .7)
# theme1$plot.line$lwd <- 2
# trellis.par.set(theme1)
# bwplot(resamps, layout = c(4, 1))
# 
# summary(model_lm)$sigm
# min(model_lasso$results$RMSE)
# summary(model_mlm)$sigma
# min(model_ridge$results$RMSE)
# 
# rmse_matrix = cbind(
#   rmse_matrix, 
#   c(
#     summary(model_lm)$sigma,
#     min(model_lasso$results$RMSE),
#     summary(model_mlm)$sigma,
#     min(model_ridge$results$RMSE)
#   )
# )
# 
# # rmse_matrix = matrix(
# #   c(
# #     summary(model_lm)$sigma,
# #     min(model_lasso$results$RMSE),
# #     summary(model_mlm)$sigma,
# #     min(model_ridge$results$RMSE)
# #   )
# # )
# 
# # rownames(rmse_matrix) = c('linear', 'lasso', 'multiple linear', 'ridge')
# # barplot(t(rmse_matrix)) 
# 
# ggplot(rmse_matrix %>% t %>% melt, aes(x = X1, y = value, col = X2)) +           # Draw line plot with ggplot2
#   geom_line()
