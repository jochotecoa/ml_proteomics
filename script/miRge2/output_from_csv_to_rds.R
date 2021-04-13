library(dplyr)

count_files = list.files(
  path = "/ngs-data/analysis/hecatos/juantxo/miRNA/miRge2-2020",
  pattern='Counts', 
  recursive=T, 
  full.names=T)

for (count_file in count_files) {
  mirna_counts = read.csv(file = count_file)
  comp = count_file %>% 
    strsplit('/') %>% 
    unlist() %>% 
    .[8] %>% 
    gsub(pattern = '-', replacement = '')
  filenaam = paste0('data/mirge2_mirna/', comp, '.miR.Counts.rds')
  saveRDS(object = mirna_counts, 
          file = filenaam)
}