source('script/functions/functions_JOA.R')

mrna_dir = '/share/analysis/hecatos/juantxo/mRNA/quant_salmon/Homo_sapiens.GRCh38.cdna.ncrna.circbase/'
tissue = 'cardiac'
compound = 't0_controls'

mrna_dir_spe = paste0(mrna_dir, tissue, '/', compound)

mrna_df = mergeFiles(files_patt = 'quant.sf', by_col = 'Name', path = mrna_dir_spe, all_true = T, recursive = T)
