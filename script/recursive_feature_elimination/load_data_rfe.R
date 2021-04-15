# load the library
source('script/functions/functions_JOA.R')
forceLibrary(c('mlbench', 'caret', 'doParallel', 'dplyr'))

cl <- makePSOCKcluster(5)
registerDoParallel(cl)
# load the data

mrna_prot_df = readRDS('data/training_data.rds')

X = mrna_prot_df[, -grep('proteomics', colnames(mrna_prot_df))] %>% as.data.frame()
Y = mrna_prot_df[, grep('proteomics', colnames(mrna_prot_df))]

normalization <- preProcess(X, verbose = T)
X <- predict(normalization, X) %>% 
  as.data.frame() %>% 
  dplyr::select(!contains('strand')) # including the strand column gives always an error, independently of the algorithm used

# sizes = seq(1, ncol(X), (ncol(X)-1)/5) %>% as.integer()
sizes = seq(1, ncol(X), 1) %>% as.integer()
if (!dir.exists('output/rfe')) {
  dir.create('output/rfe', recursive = T)
}
