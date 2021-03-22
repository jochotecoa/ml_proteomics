source('script/functions/functions_JOA.R')
forceLibrary(c('dplyr', 'tibble', 'biomaRt', 'readxl'))

biost_dir = paste0('data/biostudies/', tissue, '/')
xl_files = list.files(path = biost_dir, full.names = T)

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

saveRDS(object = biost_df, file = paste0('data/biostudies/', tissue, '/biostudies_', tissue, '.rds'))
