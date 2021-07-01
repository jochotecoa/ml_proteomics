source('script/functions/functions_JOA.R')
forceLibrary(c('dplyr', 'tibble', 'reshape2'))
addVarsProt <- function(x, fnc_list, by_str) {
  fnc_str = fnc_list[1]
  fnc = get(fnc_str)
  transf_cols = colnames(x)[!grepl(by_str, colnames(x))]
  by_lst = x[, by_str] %>% list()
  names(by_lst) = by_str
  df = x %>% 
    .[, !grepl(by_str, colnames(x)), F] %>% 
    aggregate.data.frame(by = by_lst, FUN = fnc, na.rm = TRUE) 
  colnames(df)[!grepl(by_str, colnames(df))] = paste0(transf_cols, '_', fnc_str)
  if (length(fnc_list) > 1) {
    for (fnc_str in fnc_list[-1]) {
      fnc = get(fnc_str)
      transf_cols = colnames(x)[!grepl(by_str, colnames(x))]
      by_lst = x[, by_str] %>% list()
      names(by_lst) = by_str
      x_2 = x %>% 
        .[, !grepl(by_str, colnames(x)), F] %>% 
        aggregate.data.frame(by = by_lst, FUN = fnc, na.rm = TRUE) 
      colnames(x_2)[!grepl(by_str, colnames(x_2))] = paste0(transf_cols, '_', fnc_str)
      df = x_2 %>% 
        merge.data.frame(x = df, by = by_str)
    }
  }
  return(df)
}


# Identify samples with low sequencing depth ------------------------------

if (tissue == 'cardiac') {
  mir_dir = '/ngs-data/analysis/hecatos/juantxo/miRNA/miRge2-2021/Cardiac/t0_controls_ML/miRge_altoguether/'
  
  mirna_counts = mergeFilesCsv(files_patt = 'miR.Counts', by_col = 'miRNA', path = mir_dir, all_true = T, recursive = T)
  
  colnames(mirna_counts) = colnames(mirna_counts) %>% 
    gsub(pattern='_R.*|.fastq.*', replacement='')  
} else {
  mirna_counts = read.csv('data/miRNA/miR.Counts.csv')
  mirna_counts_tmp = read.csv('data/miRNA/miR.Counts_tmp.csv')
  mirna_counts = mirna_counts %>% 
    merge.data.frame(mirna_counts_tmp)
  
}

mirna_counts = mirna_counts[!grepl('miRNAtotal', rownames(mirna_counts)), ]
mirna_counts = mirna_counts %>% 
  column_to_rownames('miRNA') %>% 
  filterSamplesBySeqDepth()
filt_cols = colnames(mirna_counts)


# Remove samples with low sequencing depth --------------------------------

if (tissue == 'cardiac') {
  mir_dir = '/ngs-data/analysis/hecatos/juantxo/miRNA/miRge2-2021/Cardiac/t0_controls_ML/miRge_altoguether/'
  
  mirna_rpm = mergeFilesCsv(files_patt = 'miR.RPM', by_col = 'miRNA', path = mir_dir, all_true = T, recursive = T)
  
  colnames(mirna_rpm) = colnames(mirna_rpm) %>% 
    gsub(pattern='_R.*|.fastq.*', replacement='')  
} else {
  mirna_rpm = read.csv('data/miRNA/miR.RPM.csv') 
  mirna_rpm_tmp = read.csv('data/miRNA/miR.RPM_tmp.csv')
  mirna_rpm = mirna_rpm %>% 
    merge.data.frame(mirna_rpm_tmp)
  
}

mirna_rpm = mirna_rpm %>% 
  column_to_rownames('miRNA')
mirna_rpm = mirna_rpm[, filt_cols]

# Rename samples ----------------------------------------------------------

if (tissue == 'hepatic') {
  biostudies_dataset = readRDS(paste0('data/biostudies/', tissue, '/biostudies_', tissue, '.rds'))
  colnames(mirna_rpm) = colnames(mirna_rpm) %>% 
    gsub(pattern = '98', replacement = '098') # VPA samples, not ISO
  colnum_mir = sapply(biostudies_dataset$`Roche ID`, grep, colnames(mirna_rpm)) %>% unlist()
  stopifnot(!any(duplicated(colnum_mir)))
  colnum_biost = sapply(paste("^",names(colnum_mir),"$", sep=""), grep, biostudies_dataset$`Roche ID`) %>% unlist()
  stopifnot(!any(duplicated(colnum_biost)))
  colnames(mirna_rpm)[colnum_mir] = biostudies_dataset$Files[colnum_biost]
  colnames(mirna_counts)[colnum_mir] = biostudies_dataset$Files[colnum_biost]
  
}


colnames(mirna_rpm) = colnames(mirna_rpm) %>% 
  gsub(pattern = '.fastq|Con_', replacement = '')
colnames(mirna_counts) = colnames(mirna_counts) %>% 
  gsub(pattern = '.fastq|Con_', replacement = '')

# Creating seq_depth_mir feature ------------------------------------------

seq_depth_mir = mirna_counts %>% 
  colSums(na.rm = T) %>% 
  data.frame(seq_depth = .) %>% 
  rownames_to_column('sample_name')

saveRDS(seq_depth_mir, paste0('data/miRNA/seq_depth_mir', tissue, '.rds'))

# Connect the miR IDs with the ENSTs --------------------------------------

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

# mrna_df = readRDS('data/mrna/mrna_df.rds')

enst_df = mrna_df %>% 
  rownames_to_column() %>% 
  dplyr::filter(grepl('ENST', rowname)) %>% 
  column_to_rownames()

enst_df = enst_df[rowSums(enst_df) > 0, , F]

swiss_mir = swiss_mir[swiss_mir$ensembl_transcript_id %in% rownames(enst_df), , F]
swiss_mir_all = swiss_mir_all[swiss_mir_all$ensembl_transcript_id %in% rownames(enst_df), , F]

mirna_swiss = mirna_rpm %>% 
  rownames_to_column('miRBase_ID') %>% 
  merge.data.frame(swiss_mir, 'miRBase_ID')
mirna_swiss_all = mirna_rpm %>% 
  rownames_to_column('miRBase_ID') %>% 
  merge.data.frame(swiss_mir_all, 'miRBase_ID')


# Create features to represent the miRNA effect ---------------------------
score_per_prot = mirna_swiss %>% 
  dplyr::select(uniprotswissprot, score)
score_per_prot_all = mirna_swiss_all %>% 
  dplyr::select(uniprotswissprot, score)

score_feats = score_per_prot %>% 
  addVarsProt(fnc_list = c('mean', 'median', 'min', 'max', 'sum', 'sd'), by_str = 'uniprotswissprot')
score_feats_all = score_per_prot_all %>% 
  addVarsProt(fnc_list = c('mean', 'median', 'min', 'max', 'sum', 'sd'), by_str = 'uniprotswissprot')

colnames(score_feats)[-1] = paste0(colnames(score_feats)[-1], '_stringent')
colnames(score_feats_all)[-1] = paste0(colnames(score_feats_all)[-1], '_all')


score_feats = score_feats[score_feats$uniprotswissprot != '', , F]
score_feats_all = score_feats_all[score_feats_all$uniprotswissprot != '', , F]

saveRDS(score_feats, paste0('data/miRNA/score_feats', tissue, '.rds'))
saveRDS(score_feats_all, paste0('data/miRNA/score_feats_all', tissue, '.rds'))


mirna_feats = mirna_swiss %>% 
  dplyr::select(!c(ensembl_gene_id,miRBase_ID, refseq_mrna, score, ensembl_transcript_id)) %>% 
  addVarsProt(fnc_list = c('mean', 'median', 'min', 'max', 'sum', 'sd'), by_str = 'uniprotswissprot')
mirna_feats_all = mirna_swiss_all %>% 
  dplyr::select(!c(miRBase_ID, refseq_mrna, score, ensembl_transcript_id)) %>% 
  addVarsProt(fnc_list = c('mean', 'median', 'min', 'max', 'sum', 'sd'), by_str = 'uniprotswissprot')

saveRDS(mirna_feats, paste0('data/miRNA/mirna_feats_', tissue, '.rds'))
saveRDS(mirna_feats_all, paste0('data/miRNA/mirna_feats_all_', tissue, '.rds'))

mirna_feats_melt = mirna_feats %>% 
  melt()
mirna_feats_melt_all = mirna_feats_all %>% 
  melt()


mirna_feats_melt = mirna_feats_melt[mirna_feats_melt$uniprotswissprot != '', , F]
mirna_feats_melt$measure = NA
mirna_feats_melt$measure[grepl('mean', mirna_feats_melt$variable)] = 'mean'
mirna_feats_melt$measure[grepl('median', mirna_feats_melt$variable)] = 'median'
mirna_feats_melt$measure[grepl('min', mirna_feats_melt$variable)] = 'min'
mirna_feats_melt$measure[grepl('max', mirna_feats_melt$variable)] = 'max'
mirna_feats_melt$measure[grepl('sum', mirna_feats_melt$variable)] = 'sum'
mirna_feats_melt$measure[grepl('sd', mirna_feats_melt$variable)] = 'sd'
mirna_feats_melt$measure = mirna_feats_melt$measure %>% 
  paste0('mirna_', .)

mirna_feats_melt_all = mirna_feats_melt_all[mirna_feats_melt_all$uniprotswissprot != '', , F]
mirna_feats_melt_all$measure = NA
mirna_feats_melt_all$measure[grepl('mean', mirna_feats_melt_all$variable)] = 'mean'
mirna_feats_melt_all$measure[grepl('median', mirna_feats_melt_all$variable)] = 'median'
mirna_feats_melt_all$measure[grepl('min', mirna_feats_melt_all$variable)] = 'min'
mirna_feats_melt_all$measure[grepl('max', mirna_feats_melt_all$variable)] = 'max'
mirna_feats_melt_all$measure[grepl('sum', mirna_feats_melt_all$variable)] = 'sum'
mirna_feats_melt_all$measure[grepl('sd', mirna_feats_melt_all$variable)] = 'sd'
mirna_feats_melt_all$measure = mirna_feats_melt_all$measure %>% 
  paste0('mirna_', .)


mirna_feats_melt$sample = mirna_feats_melt$variable %>% 
  gsub(pattern = '_mean|_median|_min|_max|_sum|_sd', replacement = '')
mirna_feats_melt_all$sample = mirna_feats_melt_all$variable %>% 
  gsub(pattern = '_mean|_median|_min|_max|_sum|_sd', replacement = '')

mirna_feats_melt$uniprot_sample = paste0(mirna_feats_melt$uniprotswissprot, '--', mirna_feats_melt$sample)
mirna_feats_melt_all$uniprot_sample = paste0(mirna_feats_melt_all$uniprotswissprot, '--', mirna_feats_melt_all$sample)

mirna_feats_dcast = mirna_feats_melt %>% 
  dplyr::select(-c(uniprotswissprot, variable, sample)) %>% 
  dcast(uniprot_sample ~ measure) 
mirna_feats_dcast_all = mirna_feats_melt_all %>% 
  dplyr::select(-c(uniprotswissprot, variable, sample)) %>% 
  dcast(uniprot_sample ~ measure) 

colnames(mirna_feats_dcast)[-1] = 
  colnames(mirna_feats_dcast)[-1] %>% 
  paste0('_stringent')
colnames(mirna_feats_dcast_all)[-grep('prot', colnames(mirna_feats_dcast_all))] = 
  colnames(mirna_feats_dcast_all)[-grep('prot', colnames(mirna_feats_dcast_all))] %>% 
  paste0('_all')

mirna_feats_dcast_log2 = mirna_feats_dcast[, -grep('uniprot_sample', colnames(mirna_feats_dcast))] %>% 
  log2()
colnames(mirna_feats_dcast_log2) = colnames(mirna_feats_dcast_log2) %>% 
  paste0('_log2')
mirna_feats_dcast = mirna_feats_dcast %>% 
  cbind.data.frame(mirna_feats_dcast_log2)
mirna_feats_dcast[mirna_feats_dcast == -Inf] = log2(2e-06)

mirna_feats_dcast_all_log2 = mirna_feats_dcast_all[, -grep('uniprot_sample', colnames(mirna_feats_dcast_all))] %>% 
  log2()
colnames(mirna_feats_dcast_all_log2) = colnames(mirna_feats_dcast_all_log2) %>% 
  paste0('_log2')
mirna_feats_dcast_all = mirna_feats_dcast_all %>% 
  cbind.data.frame(mirna_feats_dcast_all_log2)
mirna_feats_dcast_all[mirna_feats_dcast_all == -Inf] = log2(2e-06)



saveRDS(mirna_feats_dcast, paste0('data/miRNA/mirna_feats_dcast_', tissue, '.rds'))
saveRDS(mirna_feats_dcast_all, paste0('data/miRNA/mirna_feats_dcast_all_', tissue, '.rds'))
