library(caret)

tag <- read.csv("C:/Users/j.ochoteco/Downloads/tag_data.csv", row.names = 1)
tag <- as.matrix(tag)

## Select only models for regression
regModels <- tag[tag[,"Regression"] == 1,]

all <- 1:nrow(regModels)
## Seed the analysis with the SVM model
start <- grep("(lm)", rownames(regModels), fixed = TRUE)
pool <- all[all != start]

## Select 4 model models by maximizing the Jaccard
## dissimilarity between sets of models
nextMods <- maxDissim(regModels[start,,drop = FALSE], 
                      regModels[pool, ], 
                      method = "Jaccard",
                      n = 10)

rownames(regModels)[c(start, nextMods)]
