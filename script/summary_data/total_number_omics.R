# load the library
source('script/functions/functions_JOA.R')
forceLibrary(c('mlbench', 'caret', 'doParallel', 'dplyr', 'RANN', 'tibble'))
tissue = 'cardiac'
compound = 't0_controls_ML'

proteomics_cleaning = function(tissue) {
  # Get proteomics data -----------------------------------------------------
  
  Tissue = tools::toTitleCase(tissue)
  
  prot_dir = paste0('/ngs-data/data/hecatos/', Tissue, '/t0_controls/Protein/')
  if (Tissue == 'Cardiac') {
    prot_dir = paste0('/ngs-data/data/hecatos/', Tissue, '/t0_controls_ML/Protein/')
  }
  
  prot_df = mergeFiles(by_col = 'Row.Names', path = prot_dir, header = T, 
                       fill = F, sep = '\t', all_true = T)
  
  # Clean rownames and colnames of proteomics data --------------------------
  
  
  colnames(prot_df) = colnames(prot_df) %>% 
    gsub(pattern = paste0(prot_dir), replacement = '')
  
  prot_df = prot_df %>% 
    column_to_rownames('Row.Names')
  
  prot_df[, !grepl('log2', colnames(prot_df))] = 
    prot_df[, !grepl('log2', colnames(prot_df))] %>% 
    log2()
  
  colnames(prot_df)[!grepl('log2', colnames(prot_df))] = 
    colnames(prot_df)[!grepl('log2', colnames(prot_df))] %>% 
    paste0('_log2')
  
  if (!file.exists(paste0('data/biostudies/', tissue, '/biostudies_', tissue, '.rds'))) {
    file.remove(paste0('data/biostudies/', tissue, '/biostudies_', tissue, '.rds'))
    source(file = 'script/data_cleaning/renaming_proteomics_samples.R', local = T, echo = T)
  }
  
  files_rocheid = readRDS(paste0('data/biostudies/', tissue, '/biostudies_', tissue, '.rds')) %>% 
    as.data.frame() 
  
  # if (Tissue == 'Cardiac') {
  #   a = files_rocheid[grep('DF2_The_002', files_rocheid$Files), ]
  #   b = a[a$`Roche ID` %in% 958:960, ]
  #   files_rocheid = files_rocheid[-as.numeric(rownames(b)), ]
  #   
  # }
  
  
  colnames(prot_df)[unlist(files_rocheid$colnum)] = files_rocheid$Files
  
  colnames(prot_df) = 
    colnames(prot_df) %>% 
    gsub(pattern = '_/.*', 
         replacement = '')
  
  # Subselect samples -------------------------------------------------------
  
  
  t0_cols = 
    colnames(prot_df) %>% 
    subset(., grepl(pattern = 'UNTR|DF2|^DMSO|000', .)) 
  
  
  prot_df = prot_df[, !duplicated(colnames(prot_df))] %>% 
    rownames_to_column() %>% 
    dplyr::select(matches('UNTR|DF2|^DMSO|000'), rowname) %>% 
    column_to_rownames() %>% 
    filterSamplesBySeqDepth() %>% 
    rownames_to_column() %>% 
    cleanProtIds() %>% 
    dplyr::select(!rowname)
  
  prot_df = prot_df %>% 
    remove_rownames() %>% 
    column_to_rownames('uniprot_gn') 
  
  prot_df = prot_df %>% 
    apply(2, unlist) %>% 
    apply(2, as.numeric) %>% 
    data.frame(row.names = rownames(prot_df)) %>% 
    normalizeProteomics() %>% 
    rownames_to_column('uniprotswissprot')
  
  return(prot_df)
}

transcriptomics_cleaning <- function(tissue, compound) {
  # Get mRNA data -----------------------------------------------------------
  
  mrna_dir = '/share/analysis/hecatos/juantxo/mRNA/quant_salmon/Homo_sapiens.GRCh38.cdna.ncrna.circbase/'
  
  mrna_dir_spe = paste0(mrna_dir, tissue, '/', compound)
  
  mrna_df_file = paste0('../bigdata_ml_proteomics/salmon_', tissue, '.rds')
  
  if (file.exists(mrna_df_file)) {
    mrna_df = readRDS(file = mrna_df_file)
  } else {
    mrna_df = mergeFiles(files_patt = 'quant.sf', by_col = 'Name',
                         path = mrna_dir_spe, all_true = T, recursive = T,
                         header = T)
    
    saveRDS(object = mrna_df, file = mrna_df_file)
    
  }
  
  
  # Clean colnames and rownames ---------------------------------------------
  
  
  mrna_df_counts = mrna_df %>% 
    remove_rownames() %>% 
    column_to_rownames('Name') %>% 
    dplyr::select(contains('NumReads'))
  
  colnames(mrna_df_counts) = colnames(mrna_df_counts) %>% 
    gsub(pattern = paste0(mrna_dir_spe, '|', '/quant.sf|_quant|NumReads_', '|', '/'), replacement = '') 
  
  mrna_df_counts = mrna_df_counts %>% 
    filterSamplesBySeqDepth()
  
  mrna_df = mrna_df %>% 
    remove_rownames() %>% 
    column_to_rownames('Name') %>% 
    dplyr::select(contains('TPM'))
  
  colnames(mrna_df) = colnames(mrna_df) %>% 
    gsub(pattern = mrna_dir_spe, replacement = '') %>% 
    gsub(pattern = '/quant.sf|_quant|TPM_', replacement = '') %>% 
    gsub(pattern = '/', replacement = '')
  
  mrna_df = mrna_df[colnames(mrna_df_counts)]
  
  rownames(mrna_df) = rownames(mrna_df) %>% 
    gsub('\\..*', '', .) 
  
  stopifnot(all.equal(colnames(mrna_df), colnames(mrna_df_counts)))
  
  
  return(mrna_df_counts)
}

mirnaomics_cleaning <- function(tissue) {
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
  
  return(mirna_counts)
}

# all_data_hepatic = readRDS('data/whole_data_preds_na_omit.rds')


prot_df = proteomics_cleaning(tissue)
all_transcripts = transcriptomics_cleaning(tissue, compound)
all_mirna = mirnaomics_cleaning(tissue)

prot_df = prot_df %>% 
  column_to_rownames("uniprotswissprot")

all_circ = all_transcripts[grepl('circ', rownames(all_transcripts)), ]
all_linear = all_transcripts[!grepl('circ', rownames(all_transcripts)), ]

nas_sample_linear = all_linear %>% zeroToNa %>% is.na
expressed_linear_sample = !nas_sample_linear

prot_df_na_count = prot_df %>% apply(1, is.na) %>% t %>% rowSums()
prot_df_expr = prot_df[prot_df_na_count < ncol(prot_df), ]
all_circ_expr = all_circ[rowSums(all_circ, na.rm = T) > 0, ]
all_linear_expr = all_linear[rowSums(all_linear, na.rm = T) > 0, ]
all_mirna_expr = all_mirna[rowSums(all_mirna, na.rm = T) > 0, ]


prot_df_const = prot_df %>% na.omit()
all_circ_const = all_circ_expr %>% zeroToNa() %>% na.omit()
all_linear_const = all_linear_expr %>% zeroToNa() %>% na.omit()
all_mirna_const = all_mirna_expr %>% zeroToNa() %>% na.omit()


asd = character()



for (variable in 1:nrow(a)) {
  asd_i = colnames(a)[a[variable, ]] %>% paste0(collapse = '__')
  asd = c(asd, asd_i)
  print((variable/nrow(a))*100)
}

fdas = table(asd)

fdas_2 = fdas %>% names %>% sapply(strsplit, '__') %>% sapply(length)

fdas_3 = fdas %>% names %>% sapply(strsplit, '__')

ddd = fdas[order(fdas_2, decreasing = T)]
ddd = ddd[!which.max(ddd)]


fdas[order(fdas_2, decreasing = T)]

fdas[order(fdas_2, decreasing = T)][2:100] %>% .[which.max(.)]
