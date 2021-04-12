# load the library
source('script/functions/functions_JOA.R')
forceLibrary(c('mlbench', 'caret', 'doParallel', 'dplyr'))

# cl <- makePSOCKcluster(24)
registerDoParallel(5)
# load the data

mrna_prot_df = readRDS('data/training_data.rds')

X = mrna_prot_df[, -grep('proteomics', colnames(mrna_prot_df))] %>% as.data.frame()
Y = mrna_prot_df[, grep('proteomics', colnames(mrna_prot_df))]

normalization <- preProcess(X, verbose = T)
X <- predict(normalization, X) %>% 
  as.data.frame() %>% 
  .[, -2, F] # including the strand column gives always an error, independently of the algorithm used

sizes = seq(1, ncol(X), ncol(X)/5) %>% as.integer()

if (!dir.exists('output/rfe')) {
  dir.create('output/rfe', recursive = T)
}
