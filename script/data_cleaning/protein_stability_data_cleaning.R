source('script/functions/functions_JOA.R')
forceLibrary(c('dplyr', 'tibble', 'biomaRt', 'reshape', 'caTools'))

protein_stability_data = readxl::read_excel('data/protein_stability_data.xlsx')
protein_stability_data = protein_stability_data %>% 
  unique.data.frame()

mart = openMart2018()

# prot_stab_data = protein_stability_data$Symbol %>% 
#   getBM(attributes = c('external_gene_name', 'uniprotswissprot'), filters = 'external_gene_name', values = ., mart = mart) %>% 
#   unique()

prot_stab_data = protein_stability_data$`Gene ID` %>% 
  getBM(attributes = c('entrezgene', 'uniprotswissprot'), filters = 'entrezgene', values = ., mart = mart) %>% 
  unique()

# symbols = protein_stability_data$Symbol %>% unique()
# (sum(symbols %in% unique(prot_stab_data$external_gene_name)))/length(symbols)
# symbols2 = protein_stability_data$`Gene ID` %>% unique()
# (sum(symbols2 %in% unique(prot_stab_data2$entrezgene)))/length(symbols2)

prot_stab_data = prot_stab_data %>% 
  merge.data.frame(y = protein_stability_data, by.x = 'entrezgene', by.y = 'Gene ID')

prot_stab_data = prot_stab_data[prot_stab_data$uniprotswissprot != '', , F]

prot_stab_feats = prot_stab_data %>% 
  dplyr::select(-c(entrezgene, Symbol)) %>% 
  addVarsProt(fnc_list = c('mean', 'median', 'min', 'max', 'sd'), by_str = 'uniprotswissprot') %>% 
  unique.data.frame()

dir.create('data/stability/')

saveRDS(object = prot_stab_feats, file = 'data/stability/prot_stab_feats.rds')
