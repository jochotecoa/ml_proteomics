# Load library
# install.packages('h2o')
library('h2o')

# start h2o cluster
invisible(h2o.init())

# load the library
source('script/functions/functions_JOA.R')
forceLibrary(c('mlbench', 'caret', 'doParallel', 'dplyr', 'RANN'))

# cl <- parallel::makeCluster(15, setup_strategy = "sequential")
# registerDoParallel(cl)
# load the data


# file_X_train = 'data/training_data_preds_na_omit_by_prot.rds'
# file_Y = 'data/training_data_target_na_omit_by_prot.rds'

X = readRDS(file = 'data/whole_data/preds_na_omit_cardiac.rds')
Y = readRDS(file = 'data/whole_data/target_na_omit_cardiac.rds') %>% 
  unlist()

train_proteins = X %>% rownames() %>% strsplit('--') %>% sapply('[[', 1)

sizes = as.integer(seq(2, ncol(X), floor(ncol(X)/4)))[-5]

folds = groupKFold(train_proteins, k = 10)

train_rows = folds$Fold01

dataset = cbind.data.frame(X, proteomics_value = Y)



# train_rows = createDataPartition(dataset$proteomics_value, p = .8, 
#                                  list = FALSE, 
#                                  times = 1)
train = dataset[train_rows, ]
test = dataset[-train_rows, ]

# convert data as h2o type
train_h = as.h2o(train)
test_h = as.h2o(test)

# set label type
y = 'proteomics_value'
pred = setdiff(names(train), y)

#convert variables to factors
# train[,y] = as.factor(train[,y])
# test[,y] = as.factor(test[,y])

# Run AutoML for 20 base models
aml = h2o.automl(x = pred, y = y,
                 training_frame = train_h,
                 # max_models = 20,
                 seed = 1,
                 # max_runtime_secs = 3600
)

# AutoML Leaderboard
lb = aml@leaderboard
lb

# prediction result on test data
prediction_df = h2o.predict(aml@leader, test_h) %>%
  as.data.frame() %>% 
  unlist()

# create a confusion matrix
a = as.data.frame(test_h$proteomics_value) %>% unlist()
b = data.frame(obs = a, pred = prediction_df) 
d = cor(b)
d^2 %>% print()

caret::confusionMatrix(test$Species, prediction$predict)

# close h2o connection
h2o.shutdown(prompt = F)
