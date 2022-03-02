# load the library
source('script/functions/functions_JOA.R')
forceLibrary(c('mlbench', 'caret', 'doMC', 'dplyr', 'RANN'))

doMC::registerDoMC(10)
  
X = readRDS(file = 'data/whole_data/preds_na_omit_cardiac.rds')
Y = readRDS(file = 'data/whole_data/target_na_omit_cardiac.rds') %>% 
  unlist()
dataset = data.frame(X, proteomics_value = Y)

all_samples = X %>% rownames() %>% strsplit('--') %>% sapply('[[', 2)
sample_test = all_samples[1]

dataset_test = dataset[grepl(sample_test, rownames(dataset)), ]
dataset_train = dataset[!grepl(sample_test, rownames(dataset)), ]
  
dataset_train = dataset_train %>% 
  dplyr::select(circ_sum, circ_mean, circ_mean_log2, linear_density, 
                TPM_value_mean_log2, mirna_mean_stringent, 
                mirna_mean_stringent_log2, mirna_sum_stringent, Aa_M_prop, 
                Aa_D_prop, proteomics_value)
dataset_test = dataset_test %>% 
  dplyr::select(circ_sum, circ_mean, circ_mean_log2, linear_density, 
                TPM_value_mean_log2, mirna_mean_stringent, 
                mirna_mean_stringent_log2, mirna_sum_stringent, Aa_M_prop, 
                Aa_D_prop, proteomics_value)

samples_train = dataset_train %>% rownames() %>% strsplit('--') %>% sapply('[[', 2)
folds = groupKFold(samples_train, k = 10)


rfFit <- train(proteomics_value ~ ., 
               data = dataset_train,
               method = 'rf',
               trControl  = trainControl(method = 'cv', 
                                         index = folds))

# rfFit = readRDS(file = 'output/cardiac/na_omit/split_by_sample/72_features/rfFit.rds')

pred = predict(rfFit, dataset_test)
obsVSpred = dataset_test$proteomics_value %>% 
  cbind.data.frame(obs = ., pred)


rsquared_all_prot = (cor(obsVSpred[1], obsVSpred[2]))^2

obsVSpred$uniprotswissprot = obsVSpred %>% rownames() %>% strsplit('--') %>% 
  sapply('[[', 1)


library(biomaRt)

mart=useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")


prot_goterms = getBM(attributes = c('name_1006', 'uniprotswissprot'), 
                     filters = 'uniprotswissprot', 
                     values = obsVSpred$uniprotswissprot, 
                     mart = mart)

populated_goterms = table(prot_goterms$name_1006) >= 10
populated_goterms = populated_goterms[populated_goterms] %>% 
  names()


obsVSpred_goterms = obsVSpred %>% 
  merge.data.frame(prot_goterms, by = 'uniprotswissprot')

rsquared_go_results = data.frame()

# goterm = NULL, rsquared_goterm = NULL
for (goterm in populated_goterms) {
  obsVSpred_goterm = obsVSpred_goterms[obsVSpred_goterms$name_1006 == goterm, ]
  
  rsquared_goterm = (cor(obsVSpred_goterm[2], obsVSpred_goterm[3]))^2
  
  rsquared_go_result = data.frame(goterm, rsquared_goterm)
  
  rsquared_go_results = rsquared_go_results %>% 
    rbind(rsquared_go_result)
  
}

rsquared_go_results$rsquared = rsquared_go_results$pred

worst_goterms = rsquared_go_results[order(rsquared_go_results$rsquared), ] %>% head()
best_goterms = rsquared_go_results[order(rsquared_go_results$rsquared), ] %>% tail()
barplot_goterms = rbind(worst_goterms, best_goterms)


pdf()

par(mar=c(4,14,4,4))
barplot(barplot_goterms$rsquared, 
        names.arg = d$goterm, 
        horiz = T, 
        las =2, 
        xlim = c(0, 1), 
        main = 'Best and worse GO terms by Rsquared')# rsquared_go_results %<>% na.omit()


dev.off()

tiff("output/plots/best_worse_goterms.tiff", units="in", width=10, height=5, res=300)
par(mar=c(4,14,4,8))
barplot(barplot_goterms$rsquared, 
        names.arg = d$goterm, 
        horiz = T, 
        las =2, 
        xlim = c(0, 1), 
        main = 'Best and worse predicted proteins grouped by GO terms (Rsquared)')# rsquared_go_results %<>% na.omit()

dev.off()
