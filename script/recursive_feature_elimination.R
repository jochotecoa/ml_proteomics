# load the library
forceLibrary(c('mlbench', 'caret', 'doParallel'))
forceLibrary(c('ipred', 'plyr', 'e1071')) # Needed for bagged trees
forceLibrary('kernlab')

cl <- makePSOCKcluster(24)
registerDoParallel(cl)
# load the data
X = mrna_prot_df[, -grep('proteomics', colnames(mrna_prot_df))] %>% as.data.frame()
Y = mrna_prot_df[, grep('proteomics', colnames(mrna_prot_df))]

normalization <- preProcess(X, verbose = T)
X <- predict(normalization, X)
X <- as.data.frame(X)


# # define the control using a random forest selection function
# control <- rfeControl(functions=rfFuncs, method="cv", number=10, verbose = T)
# # run the RFE algorithm
# results <- rfe(X, Y, sizes=c(1:ncol(X)), rfeControl=control)
# # summarize the results
# print(results)
# # list the chosen features
# predictors(results)
# # plot the results
# plot(results, type=c("g", "o"))


sizes = 1:ncol(X)

# if (file.exists('lmProfile.rds')) {
#   lmProfile = readRDS('lmProfile.rds')
# } else {
  lmProfile <- rfe(X, Y,
                   sizes = sizes,
                   rfeControl = rfeControl(functions = lmFuncs, verbose = T))
  saveRDS(lmProfile, 'lmProfile.rds')
  
# }

# normalization <- preProcess(X)
# X_pp <- predict(normalization, X)
# X_pp <- as.data.frame(X_pp)

# pp_lmProfile <- rfe(X_pp, Y,
#                  sizes = sizes,
#                  rfeControl = rfeControl(functions = lmFuncs, verbose = T))
# saveRDS(pp_lmProfile, 'pp_lmProfile.rds')



if (file.exists('lmProfile2.rds')) {
  lmProfile2 = readRDS('lmProfile2.rds')
} else {
  
  lmProfile2 <- rfe(X, Y,
                    sizes = sizes,
                    rfeControl = rfeControl(functions = lmFuncs,
                                            rerank = TRUE, verbose = T))
  saveRDS(lmProfile2, 'lmProfile2.rds')
}

if (file.exists('rfProfile.rds')) {
  rfProfile = readRDS('rfProfile.rds')
} else {
  rfProfile <- rfe(X, Y,
                   sizes = sizes,
                   rfeControl = rfeControl(functions = rfFuncs, verbose = T))
  saveRDS(rfProfile, 'rfProfile.rds')
  
}

if (file.exists('bagProfile.rds')) {
  bagProfile = readRDS('bagProfile.rds')
} else {
  
  bagProfile <- rfe(X, Y,
                    sizes = sizes,
                    rfeControl = rfeControl(functions = treebagFuncs, verbose = T))
  saveRDS(bagProfile, 'bagProfile.rds')
}

svmProfile <- rfe(X, Y,
                  sizes = sizes,
                  rfeControl = rfeControl(functions = caretFuncs, verbose = T, number = 2),
                  ## pass options to train()
                  method = "svmRadial")

xyplot(lmProfile$results$RMSE + lmProfile2$results$RMSE + rfProfile$results$RMSE + bagProfile$results$RMSE  ~
         lmProfile$results$Variables, xlab = 'Number of Variables', ylab = 'Root-mean-square deviation (RMSE)',
       type = c("g", "p", "l"),
       auto.key = TRUE)

