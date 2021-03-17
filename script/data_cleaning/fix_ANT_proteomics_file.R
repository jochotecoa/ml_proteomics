source('script/functions/functions_JOA.R')
forceLibrary(c('dplyr', 'tibble', 'biomaRt'))

prot_dir = '/ngs-data/data/hecatos/Cardiac/t0_controls/Protein/'

files = list.files(path = prot_dir, full.names = T)

# Hecatos_Cardio_Px_ANTvsFluctDMSO_log2.txt

prot1 = read.table(file = files[9], header = T, fill = F, sep = ' ')

colnames(prot1)[1] = 'Row.Names'

write.table(x = prot1, file = files[9], sep = '\t')