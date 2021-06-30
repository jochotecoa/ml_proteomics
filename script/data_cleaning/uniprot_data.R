prot_df$uniprotswissprot %>% unique %>% 
  write.table('../prot_ids.txt', row.names=F, col.names=F, quote=F)

# Input your protein list in UniProt and select the desired features
if (tissue == 'hepatic') {
  uniprot.yourlist <- read.delim2("data/uniprot-yourlist M202104068471C63D39733769F8E060B506551E121ABF7EE(1).tab")
}

if (tissue == 'cardiac') {
  uniprot.yourlist <- read.delim2('data/uniprot_output/uniprot-yourlist_M20210625_cardiac.tab')
}
uniprot.yourlist_ori = uniprot.yourlist


# Count the number of aminoacids per protein ------------------------------

seq_unip = uniprot.yourlist$Sequence[1] %>% unlist() %>% as.character %>% 
  strsplit(split = '') %>% unlist() %>% table() %>% as.data.frame()
colnames(seq_unip) = c('aa', as.character(unlist(uniprot.yourlist$Entry[1])))

for (i in seq(2, nrow(uniprot.yourlist))) {
  seq_unip_i = uniprot.yourlist$Sequence[i] %>% unlist() %>% as.character %>% 
    strsplit(split = '') %>% unlist() %>% table() %>% as.data.frame()
  colnames(seq_unip_i) = c('aa', as.character(unlist(uniprot.yourlist$Entry[i])))
  
  seq_unip = merge.data.frame(x = seq_unip, y = seq_unip_i, by = 'aa', all = T)
  
}

seq_unip_ori = seq_unip
# seq_unip = seq_unip_ori

seq_unip = seq_unip %>% t %>% gsub(pattern = ' ', replacement = '') %>% 
  as.data.frame()

colnames(seq_unip) = seq_unip[1, ] %>% unlist() %>% as.character() %>% 
  paste0('Aa_', .)

seq_unip = seq_unip[-1, ]

seq_rnames = rownames(seq_unip)

seq_unip = seq_unip %>% apply(2, as.numeric) %>% naToZero() %>% as.data.frame()

rownames(seq_unip) = seq_rnames

uniprot.yourlist = uniprot.yourlist %>% 
  merge.data.frame(y = rownames_to_column(seq_unip, 'Entry'), by = 'Entry')

# Remove useless features -------------------------------------------------


uniprot.yourlist = uniprot.yourlist[, !(colnames(uniprot.yourlist) %in% c('Entry.name', 'Status', 'Fragment', 'Gene.names', 'RNA.editing', 'Sequence.uncertainty', 'Mass.spectrometry', 'Non.adjacent.residues', 'Sequence'))]


# Check the classes of the features ---------------------------------------

a = NULL

for(i in colnames(uniprot.yourlist)) {a = c(a, uniprot.yourlist[, i] %>% sapply(class) %>% table)}

a = names(a)
names(a) = colnames(uniprot.yourlist)



# Fix other features ------------------------------------------------------



uniprot.yourlist$Mass = uniprot.yourlist$Mass %>% gsub(',', '', .) %>% 
  unlist %>% as.numeric

uniprot.yourlist$Version..sequence. = uniprot.yourlist$Version..sequence. %>% 
  as.factor()

# Create proportion of aminoacids features --------------------------------


aa_cols = uniprot.yourlist %>% 
  dplyr::select(contains('Aa'))
aa_prop_cols = aa_cols / uniprot.yourlist$Length
colnames(aa_prop_cols) = colnames(aa_prop_cols) %>% 
  paste0('_prop')

uniprot.yourlist = uniprot.yourlist %>% 
  cbind.data.frame(aa_prop_cols)


# Generate linear density -------------------------------------------------

uniprot.yourlist$linear_density = uniprot.yourlist$Mass / uniprot.yourlist$Length

saveRDS(object = uniprot.yourlist, file = paste0('data/uniprot_yourlist', tissue, '.rds'))
