source('script/functions/functions_JOA.R')
forceLibrary(c('dplyr', 'tibble', 'biomaRt', 'readxl'))

biost_dir = paste0('data/biostudies/', tissue, '/')
xl_files = list.files(path = biost_dir, full.names = T)
xl_files = xl_files[!grepl('biostudies_cardiac.rds', xl_files)]

biost_df = data.frame()



for (xl_file in xl_files) {
  biostudies_xl <- read_excel(xl_file, 
                              skip = 20) %>% 
    .[, 1:2, F]
  
  if (any(grepl('Files', unlist(biostudies_xl[, 1])))) {
    row_colnames = grep('Files', unlist(biostudies_xl[, 1]))
    colnames(biostudies_xl) = biostudies_xl[row_colnames, ]
    biostudies_xl = biostudies_xl[-(1:row_colnames), ]
  }
  
  biostudies_xl$Files = biostudies_xl$Files %>% 
    gsub(pattern = paste0('Total_RNA/', Tissue, '/.*/.*/|_R.*'), replacement = '')
  
  stopifnot(!any(grepl('Total', biostudies_xl$Files)))
  
  biostudies_xl = biostudies_xl %>% 
    unique()
  
  biost_df = rbind.data.frame(biost_df, biostudies_xl)
  
}

biost_df = biost_df %>% 
  na.omit() %>% 
  filter(grepl('_', Files))

files_rocheid = biost_df

if (tissue == 'hepatic') {
  files_rocheid = files_rocheid %>% 
    dplyr::filter(!grepl('CYC_Tox_000.*', Files))
}

if (tissue == 'cardiac') {
  files_rocheid$`Roche ID`[grep('AMI_.*_000', files_rocheid$Files)] = 
    files_rocheid$`Roche ID`[grep('AMI_.*_000', files_rocheid$Files)] %>% 
    paste0('_e1')
  
  files_rocheid$`Roche ID`[grep('DOC_.*_000', files_rocheid$Files)] = 
    files_rocheid$`Roche ID`[grep('DOC_.*_000', files_rocheid$Files)] %>% 
    paste0('_e2')
}

files_rocheid = files_rocheid %>% 
  unique.data.frame()

files_rocheid$colnum = files_rocheid$`Roche ID` %>% 
  paste0('_', .) %>% 
  as.data.frame() %>% 
  apply(MARGIN = 1, FUN = grep, colnames(prot_df)) 

files_rocheid$len = files_rocheid$colnum %>% 
  sapply(length)

files_rocheid$colnum[files_rocheid$len == 0] = 
  files_rocheid$`Roche ID`[files_rocheid$len == 0] %>% 
  paste0('_0', ., '_') %>% 
  as.data.frame() %>% 
  apply(MARGIN = 1, FUN = grep, colnames(prot_df)) 

files_rocheid$len = files_rocheid$colnum %>% 
  sapply(length)

colnum_dupl = files_rocheid$colnum[files_rocheid$len != 0] %>% duplicated %>% files_rocheid$colnum[files_rocheid$len != 0][.] %>% unlist

if (length(colnum_dupl) != 0) {
  files_rocheid$colnum[files_rocheid$len == 0] = 
    files_rocheid$`Roche ID`[files_rocheid$len == 0] %>% 
    paste0('_0', .) %>% 
    as.data.frame() %>% 
    apply(MARGIN = 1, FUN = grep, colnames(prot_df)) 
  
  files_rocheid$len = files_rocheid$colnum %>% 
    sapply(length)
  
}

if (tissue == 'hepatic') {
  files_rocheid$colnum[files_rocheid$len == 0] = 
    files_rocheid$`Roche ID`[files_rocheid$len == 0] %>% 
    paste0('_', ., '_') %>% 
    as.data.frame() %>% 
    apply(MARGIN = 1, FUN = grep, colnames(prot_df)) 
  
  files_rocheid$len = files_rocheid$colnum %>% 
    sapply(length)
}

colnum_dupl = files_rocheid$colnum[files_rocheid$len != 0] %>% duplicated %>% files_rocheid$colnum[files_rocheid$len != 0][.] %>% unlist

stopifnot(length(colnum_dupl) == 0)

# files_rocheid[files_rocheid$colnum %in% colnum_dupl, ]

files_rocheid = files_rocheid[files_rocheid$len != 0, ]


saveRDS(object = files_rocheid, file = paste0('data/biostudies/', tissue, '/biostudies_', tissue, '.rds'))
