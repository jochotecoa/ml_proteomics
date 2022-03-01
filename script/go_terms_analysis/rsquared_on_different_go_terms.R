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

pred = predict(rfFit, dataset_sample)
obsVSpred = dataset_sample$proteomics_value %>% 
  cbind.data.frame(obs = ., pred)


a = cor(obsVSpred[1], obsVSpred[2])
