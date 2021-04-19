bagPv1 = readRDS('../rfe/bagProfile.rds')
bagPredsv1 = predictors(bagP) %>% unlist %>% as.character()
lmPv1 = readRDS('../rfe/lmProfile.rds')
lmPredsv1 = predictors(lmP) %>% unlist %>% as.character()
lmP2v1 = readRDS('../rfe/lmProfile2.rds')
lmP2redsv1 = predictors(lmP2) %>% unlist %>% as.character()
rfPv1 = readRDS('../rfe/rfProfile.rds')
rfPredsv1 = predictors(rfP) %>% unlist %>% as.character()
svmPv1 = readRDS('../rfe/svmProfile.rds')
svmPredsv1 = predictors(svmP) %>% unlist %>% as.character()
allPredsv1 = rfP$variables$var %>% unique

bagPv2 = readRDS('../output_rfe/bagProfile.rds')
bagPredsv2 = predictors(bagP) %>% unlist %>% as.character()
lmPv2 = readRDS('../output_rfe/lmProfile.rds')
lmPredsv2 = predictors(lmPv2) %>% unlist %>% as.character()
lmP2v2 = readRDS('../output_rfe/lmProfile2.rds')
lmP2redsv2 = predictors(lmP2) %>% unlist %>% as.character()
rfPv2 = readRDS('../output_rfe/rfProfile.rds')
rfPredsv2 = predictors(rfP) %>% unlist %>% as.character()
svmPv2 = readRDS('../output_rfe/svmProfile.rds')
svmPredsv2 = predictors(svmP) %>% unlist %>% as.character()
allPredsv2 = rfP$variables$var %>% unique



lt = list(bagPreds = bagPreds,
         lmPreds = lmPreds,
         lmP2reds = lmP2reds,
         rfPreds = rfPreds,
         svmPreds = svmPreds,
         allPreds = allPreds)
lt = lt %>% sapply(gsub, pattern = '\\.', replacement = '')
preds_uniq = rfP$variables$var %>% unique %>% sapply(gsub, pattern = '\\.', replacement = '')
lt_df = data.frame(matrix(nrow = length(lt), ncol = length(preds_uniq)))
colnames(lt_df) = preds_uniq
rownames(lt_df) = names(lt)

for(coln in 1:ncol(lt_df)){lt_df[,coln] = grepl(paste('\\b', colnames(lt_df)[coln], '\\b', sep = ''), lt)}
lt_df_t = lt_df %>% t
lt_df[, colSums(lt_df) == 5]

forceLibrary(c('UpSetR', 'tibble', 'reshape2', 'turner'), dependencies = T)
lt_df = fromList(lt)
rownames(lt_df) = rfP$variables$var %>% unique

upset(fromList(lt), nsets = 6, order.by = 'freq')
