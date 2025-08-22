#
# AI.R
#
# copyright (c) 2020-2030 - Danny Arends
#
#  Use random forest to predict longevity based on the 3 major covariates (sex, cohort, treatment)
#

setwd("/home/rqdt9/Github/UM-HET3")
source("adjustXprobs.R")
setwd("/home/rqdt9/OneDrive/Documents/HU-Berlin/UM-HET3/files")

# Read cross object
library(qtl)
library(caret)
mcross <- read.cross(format="csvr", file="um-het3-rqtl.csvr", genotypes=NULL, na.strings=c("-", "NA"))
mcross <- calc.genoprob(mcross, step = 0)
mcross <- adjustXprobs(mcross)
gtsp <- pull.genoprob(mcross)

full <- cbind(pull.pheno(mcross)[, c("sex", "cohort", "treatment", "longevity")], gtsp)

training <- sample(nrow(full), nrow(full) * .50)
testing <- which(!1:nrow(full) %in% training)

train <- full[training, ]
test  <- full[testing, ]

fitControl <- trainControl(method = "cv")
rfFit <- train(longevity ~ ., data = train, 
                 method = "rf", 
                 trControl = fitControl, na.action = na.pass)

plot(test[,"longevity"], predict(rfFit, newdata = test))
varImp(rfFit)

