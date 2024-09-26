#Model selection
library(MASS)

options <- commandArgs(trailingOnly = TRUE)

csvName = options[1]
outputName = options[2]

#Reading other covar
info <- read.table(csvName, sep = "\t", header = T)

#Merging covar with PCA and removing the ID to variable selection
allCovarNonID = subset(info,select =  -c(IID))

#Creating two models: Null and Full
modelFull <- glm(DISEASE~ ., data = allCovarNonID)
modelNull <- glm(DISEASE ~ 1, data = allCovarNonID)

#Selecting variables
X = stepAIC(modelFull, direction = "both")
matrix = model.matrix(X)
names = colnames(matrix)[-1]

#Output
write.table(names, paste(outputName,"_variables.tsv", sep = ""), row.names=FALSE, col.names=FALSE, quote=FALSE, sep = "\t")