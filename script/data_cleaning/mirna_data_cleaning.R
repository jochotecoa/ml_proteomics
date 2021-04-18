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

mirna_counts = read.csv('data/miRNA/miR.Counts.csv')
mirna_counts = mirna_counts %>% 
  column_to_rownames('miRNA') %>% 
  filterSamplesBySeqDepth()
filt_cols = colnames(mirna_counts)


# Remove samples with low sequencing depth --------------------------------

mirna_rpm = read.csv('data/miRNA/miR.RPM.csv') %>% 
  column_to_rownames('miRNA')
mirna_rpm = mirna_rpm[, filt_cols]


# Rename samples ----------------------------------------------------------

biostudies_dataset = readRDS('data/biostudies/hepatic/biostudies_hepatic.rds')
colnum_mir = sapply(biostudies_dataset$`Roche ID`, grep, colnames(mirna_rpm)) %>% unlist()
stopifnot(!any(duplicated(colnum_mir)))
colnum_biost = sapply(paste("^",names(colnum_mir),"$", sep=""), grep, biostudies_dataset$`Roche ID`) %>% unlist()
stopifnot(!any(duplicated(colnum_biost)))
colnames(mirna_rpm)[colnum_mir] = biostudies_dataset$Files[colnum_biost]
colnames(mirna_counts)[colnum_mir] = biostudies_dataset$Files[colnum_biost]

colnames(mirna_rpm) = colnames(mirna_rpm) %>% 
  gsub(pattern = '.fastq|Con_', replacement = '')
colnames(mirna_counts) = colnames(mirna_counts) %>% 
  gsub(pattern = '.fastq|Con_', replacement = '')

# Creating seq_depth_mir feature ------------------------------------------

seq_depth_mir = mirna_counts %>% 
  colSums() %>% 
  data.frame(seq_depth = .) %>% 
  rownames_to_column('sample_name')

saveRDS(seq_depth_mir, 'data/miRNA/seq_depth_mir.rds')

# Connect the miR IDs with the ENSTs --------------------------------------

enst_mir = readRDS('data/mti/mrna_mir/mirna_transcript_strict_gene.rds') 

unzip(zipfile = 'data/mti/mrna_mir/miRDB_hsa_v6.0_prediction_result.txt.zip', exdir = 'data/mti/mrna_mir/')
enst_mir_all = read.table('data/mti/mrna_mir/miRDB_hsa_v6.0_prediction_result.txt')
colnames(enst_mir_all) = c('miRBase_ID', 'refseq_mrna', 'score')


mart = openMart2018()
refseq_swiss = getBM(filters = 'refseq_mrna', 
                    values = unique(enst_mir_all$refseq_mrna), 
                    attributes = c('refseq_mrna', 'uniprotswissprot'), 
                    mart = mart)

swiss_mir = enst_mir %>% 
  merge.data.frame(refseq_swiss, 'refseq_mrna')
swiss_mir_all = enst_mir_all %>% 
  merge.data.frame(refseq_swiss, 'refseq_mrna')

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


saveRDS(score_feats, 'data/miRNA/score_feats.rds')
saveRDS(score_feats_all, 'data/miRNA/score_feats_all.rds')


mirna_feats = mirna_swiss %>% 
  dplyr::select(!c(ensembl_gene_id,miRBase_ID, refseq_mrna)) %>% 
  addVarsProt(fnc_list = c('mean', 'median', 'min', 'max', 'sum', 'sd'), by_str = 'uniprotswissprot')
mirna_feats_all = mirna_swiss_all %>% 
  dplyr::select(!c(miRBase_ID, refseq_mrna)) %>% 
  addVarsProt(fnc_list = c('mean', 'median', 'min', 'max', 'sum', 'sd'), by_str = 'uniprotswissprot')

saveRDS(mirna_feats, 'data/miRNA/mirna_feats.rds')
saveRDS(mirna_feats_all, 'data/miRNA/mirna_feats_all.rds')

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

saveRDS(mirna_feats_dcast, 'data/miRNA/mirna_feats_dcast.rds')
saveRDS(mirna_feats_dcast_all, 'data/miRNA/mirna_feats_dcast_all.rds')
