source('script/functions/functions_JOA.R')
forceLibrary(c('dplyr', 'tibble', 'biomaRt', 'reshape', 'caret', 'caTools'))

tissue = 'cardiac'
compound = 't0_controls_ML'

source('script/data_cleaning/mrna_data_cleaning.R')

source('script/data_cleaning/protein_data_cleaning.R')

source('script/data_cleaning/transcript_expression_per_protein.R')

source('script/data_cleaning/testing_sample_compability_between_transcrx_and_protx.R')

# Combine mRNA with proteins ----------------------------------------------

mrna_prot_df = mrna_unip_df %>% 
  merge.data.frame(y = prot_df, by = 'uniprot_sample') # all.y = T if there's something else aside TPM

mrna_prot_df %>% dplyr::select(where(is.numeric)) %>% na.omit %>% cor

mrna_prot_df = mrna_prot_df %>% 
  dplyr::select(!matches('\\.y'))

colnames(mrna_prot_df) = colnames(mrna_prot_df) %>% 
  gsub(pattern = '\\.x', replacement = '')

mrna_prot_df = mrna_prot_df %>% 
  dplyr::select(!c(uniprotswissprot, sample_name)) %>% 
  remove_rownames() %>% 
  column_to_rownames('uniprot_sample')

mrna_prot_df$strand_median = mrna_prot_df$strand_median %>% as.factor()
mrna_prot_df$strand_min = mrna_prot_df$strand_min %>% as.factor()
mrna_prot_df$strand_max = mrna_prot_df$strand_max %>% as.factor()


source('script/data_cleaning/combine_mirna_with_protx.R')

source('script/data_cleaning/combine_circrna_with_protx.R')

source('script/data_cleaning/combine_prot_stability_with_protx.R')

clean_rows = mrna_prot_df[, -grep('proteomics_value', colnames(mrna_prot_df))] %>% 
  na.omit() %>% 
  rownames()
  
dataset = mrna_prot_df[clean_rows, ]

dataset_train = dataset[!is.na(dataset$proteomics_value), ]
dataset_test = dataset[is.na(dataset$proteomics_value), ]

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


dataset_train = dataset_train[grep('UNTR', rownames(dataset_train)), , F]
dataset_test = dataset_test[grep('UNTR', rownames(dataset_test)), , F]


doMC::registerDoMC(10)

rfFit <- train(proteomics_value ~ ., 
               data = dataset_train,
               method = 'rf',
               trControl  = trainControl(method = 'cv', verbose = T))



proteins_train = dataset_train %>% rownames() %>% strsplit('--') %>% 
  sapply('[[', 1) %>% unique()
proteins_test = dataset_test %>% rownames() %>% strsplit('--') %>% 
  sapply('[[', 1)



pred = predict(rfFit, dataset_test)
obs = dataset_train['proteomics_value']

# pred_untr = pred[grep('UNTR', names(pred))]
# obs_untr = obs[grep('UNTR', rownames(obs)), , F]

protein_pred_untr = pred_untr %>% names() %>% strsplit('--') %>% 
  sapply('[[', 1) %>% table
protein_obs_untr = obs_untr %>% rownames() %>% strsplit('--') %>% 
  sapply('[[', 1) %>% table


protein_untr_Nimp = protein_obs_untr[protein_obs_untr == 17] %>% 
  sample(1) %>% 
  names()

protein_obs_untr_Nimp = obs_untr[grep(protein_untr_Nimp, rownames(obs_untr)), , F]
protein_obs_untr_Nimp[order(rownames(protein_obs_untr_Nimp)), , F] %>% print()


for (Nimp in 1:17) {
  protein_untr_Nimp = protein_pred_untr[protein_pred_untr == Nimp] 
  if (length(protein_untr_Nimp) == 0) {
  }
  protein_untr_Nimp = protein_untr_Nimp %>% 
    sample(1) %>% 
    names()
  
  protein_pred_untr_Nimp = pred_untr[grep(protein_untr_Nimp, names(pred_untr))] 
  names(protein_pred_untr_Nimp) = names(protein_pred_untr_Nimp) %>% 
    paste0('_imputed')
  protein_obs_untr_Nimp = obs_untr[grep(protein_untr_Nimp, rownames(obs_untr)), , F]
  
  protein_obs_untr_Nimp_num = protein_obs_untr_Nimp %>% unlist()
  
  names(protein_obs_untr_Nimp_num) = rownames(protein_obs_untr_Nimp)
  
  protein_untr_Nimp = c(protein_pred_untr_Nimp, protein_obs_untr_Nimp_num)
  
  protein_untr_Nimp[order(names(protein_untr_Nimp))] %>% as.data.frame %>% print()
  
  readline(prompt="Press [enter] to continue")
  
}



