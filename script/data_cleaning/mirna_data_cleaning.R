source('script/functions/functions_JOA.R')
forceLibrary(c('dplyr', 'tibble', 'reshape2'))


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

mart = openMart2018()
refseq_swiss = getBM(filters = 'refseq_mrna', 
                    values = unique(enst_mir$refseq_mrna), 
                    attributes = c('refseq_mrna', 'uniprotswissprot'), 
                    mart = mart)
swiss_mir = enst_mir %>% 
  merge.data.frame(refseq_swiss, 'refseq_mrna')

mirna_swiss = mirna_rpm %>% 
  rownames_to_column('miRBase_ID') %>% 
  merge.data.frame(swiss_mir, 'miRBase_ID')

# Create features to represent the miRNA effect ---------------------------
score_per_prot = mirna_swiss %>% 
  dplyr::select(uniprotswissprot, score)


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

score_feats = score_per_prot %>% 
  addVarsProt(fnc_list = c('mean', 'median', 'min', 'max', 'sum', 'sd'), by_str = 'uniprotswissprot')


mirna_feats = mirna_swiss %>% 
  dplyr::select(!c(ensembl_gene_id,miRBase_ID, refseq_mrna)) %>% 
  addVarsProt(fnc_list = c('mean', 'median', 'min', 'max', 'sum', 'sd'), by_str = 'uniprotswissprot')
 
saveRDS(mirna_feats, 'data/miRNA/mirna_feats.rds')

mirna_feats_melt = mirna_feats %>% 
  melt()
