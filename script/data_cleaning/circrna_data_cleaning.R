source('script/functions/functions_JOA.R')
forceLibrary(c('dplyr', 'tibble', 'biomaRt', 'reshape', 'caret', 'caTools'))

tissue = 'hepatic'
compound = 't0_controls_ML'

source('script/data_cleaning/mrna_data_cleaning.R')


circ_df = mrna_df %>% 
  rownames_to_column() %>% 
  dplyr::filter(!grepl('ENST', rowname)) %>% 
  dplyr::rename(circBase_ID = rowname)

cir_mir = readRDS('../ciRmiR_unique_highsponging.rds')
cir_mir2 = readRDS('../ciRmiR_unique_nseedfrequency.rds')

cir_mir[, 'miRBase_ID'] = cir_mir[, 'Var1'] %>% 
  as.character() %>% 
  strsplit(' ') %>% 
  sapply('[[', 1) %>% 
  as.character()
cir_mir2['miRBase_ID'] = cir_mir2[, 'Var1'] %>% 
  as.character() %>% 
  strsplit(' ') %>% 
  sapply('[[', 1) %>% 
  as.character()

cir_mir['circBase_ID'] = cir_mir[, 'Var1'] %>% 
  as.character() %>% 
  strsplit(' ') %>% 
  sapply('[[', 2) %>% 
  as.character()
cir_mir2['circBase_ID'] = cir_mir2[, 'Var1'] %>% 
  as.character() %>% 
  strsplit(' ') %>% 
  sapply('[[', 2) %>% 
  as.character()


colnames(cir_mir)[2] = 'circ_score'
colnames(cir_mir2)[2] = 'circ_score'

dir.create('../ml_bigdata/')

saveRDS(cir_mir2, '../ml_bigdata/cir_mir2.rds')

circ_df['circBase_ID'] = circ_df[, 'circBase_ID'] %>% 
  cleanCircNames()

circ_df = circ_df[rowSums(circ_df[-1]) > 0, , F]

# miRNA filtering ---------------------------------------------------------


mirna_counts = read.csv('data/miRNA/miR.Counts.csv')
mirna_counts = mirna_counts %>% 
  column_to_rownames('miRNA') %>% 
  filterSamplesBySeqDepth()
filt_cols = colnames(mirna_counts)

mirna_rpm = read.csv('data/miRNA/miR.RPM.csv') %>% 
  column_to_rownames('miRNA')
mirna_rpm = mirna_rpm[, filt_cols]
mirna_rpm = mirna_rpm[rowSums(mirna_rpm) > 0, , F]
mirna_rpm = mirna_rpm %>% 
  rownames_to_column('miRBase_ID')

cir_mir = cir_mir[cir_mir$miRBase_ID %in% mirna_rpm$miRBase_ID, , F]
cir_mir2 = cir_mir2[cir_mir2$miRBase_ID %in% mirna_rpm$miRBase_ID, , F]


# circRNA filtering -------------------------------------------------------



cir_mir = cir_mir[cir_mir$circBase_ID %in% circ_df$circBase_ID, , F]
cir_mir2 = cir_mir2[cir_mir2$circBase_ID %in% circ_df$circBase_ID, , F]




# saveRDS(circ_df_stringent, '../circ_df_stringent.rds')
# saveRDS(circ_df_all, '../circ_df_all.rds')
# circ_df_stringent = readRDS('../circ_df_stringent.rds')
# circ_df_all = readRDS('../circ_df_all.rds')


# Get protein-miRNA table -------------------------------------------------


enst_mir = readRDS('data/mti/mrna_mir/mirna_transcript_strict_gene.rds') 

unzip(zipfile = 'data/mti/mrna_mir/miRDB_hsa_v6.0_prediction_result.txt.zip', exdir = 'data/mti/mrna_mir/')
enst_mir_all = read.table('data/mti/mrna_mir/miRDB_hsa_v6.0_prediction_result.txt')
file.remove('data/mti/mrna_mir/miRDB_hsa_v6.0_prediction_result.txt')
colnames(enst_mir_all) = c('miRBase_ID', 'refseq_mrna', 'score')


mart = openMart2018()
refseq_swiss = getBM(filters = 'refseq_mrna', 
                     values = unique(enst_mir_all$refseq_mrna), 
                     attributes = c('refseq_mrna', 'ensembl_transcript_id', 'uniprotswissprot'), 
                     mart = mart)

swiss_mir = enst_mir %>% 
  merge.data.frame(refseq_swiss, 'refseq_mrna')
swiss_mir_all = enst_mir_all %>% 
  merge.data.frame(refseq_swiss, 'refseq_mrna')

enst_df = mrna_df %>% 
  rownames_to_column() %>% 
  dplyr::filter(grepl('ENST', rowname)) %>% 
  column_to_rownames()

enst_df = enst_df[rowSums(enst_df) > 0, , F]

swiss_mir = swiss_mir[swiss_mir$ensembl_transcript_id %in% rownames(enst_df), , F]
swiss_mir_all = swiss_mir_all[swiss_mir_all$ensembl_transcript_id %in% rownames(enst_df), , F]

cir_mir = cir_mir[-1]
cir_mir2 = cir_mir2[-1]

swiss_mir = swiss_mir %>% 
  dplyr::select(-c(refseq_mrna, ensembl_gene_id, ensembl_transcript_id, score)) %>% 
  unique.data.frame()
swiss_mir_all = swiss_mir_all %>% 
  dplyr::select(-c(refseq_mrna, score, ensembl_transcript_id)) %>% 
  unique.data.frame()


# Filter by proteomics ----------------------------------------------------

source('script/data_cleaning/protein_data_cleaning.R')


swiss_mir = swiss_mir[swiss_mir$uniprotswissprot %in% prot_df$uniprotswissprot, ]
swiss_mir_all = swiss_mir_all[swiss_mir_all$uniprotswissprot %in% prot_df$uniprotswissprot, ]

cir_mir = cir_mir %>% 
  merge.data.frame(swiss_mir, 'miRBase_ID')
cir_mir2 = cir_mir2 %>% 
  merge.data.frame(swiss_mir_all, 'miRBase_ID')

cir_mir = cir_mir %>% 
  dplyr::select(-c(miRBase_ID)) %>% 
  unique.data.frame()

cir_mir = cir_mir[cir_mir$uniprotswissprot %in% prot_df$uniprotswissprot, ]

saveRDS(cir_mir, '../ml_bigdata/cir_mir.rds')

# Combine protein IDs with circ IDs ---------------------------------------

circ_df_stringent = circ_df %>% 
  merge.data.frame(cir_mir, 'circBase_ID')
# circ_df_all = circ_df %>% 
#   merge.data.frame(cir_mir2, 'circBase_ID')

saveRDS(cir_mir, '../ml_bigdata/circ_prot_scores_stringent.rds')  
# circ_scores_all = cir_mir2 

saveRDS(circ_df_stringent, '../ml_bigdata/circ_df_stringent.rds')

source('script/data_cleaning/circrna_measures.R')






