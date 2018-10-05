################
# This code tests some of the checkDataObject functionality.
# Needs to be updated to check additional functionality.
################

library(RUnit)

#############
# Check the base cases that should on tinyMetaObject
#############
testTinyMeta <- list()
testTinyMeta$originalData <- tinyMetaObject$originalData
#Object is pre-analysis, but not pre- or post-filter
checkTrue(checkDataObject(testTinyMeta, "Meta", "Pre-Analysis"))
checkTrue(!(checkDataObject(testTinyMeta, "Meta", "Pre-Filter")))
checkTrue(!(checkDataObject(testTinyMeta, "Meta", "Post-Filter")))


testTinyResults <- runMetaAnalysis(testTinyMeta, maxCores=1)
checkTrue(checkDataObject(testTinyResults, "Meta", "Pre-Analysis"))
checkTrue(checkDataObject(testTinyResults, "Meta", "Pre-Filter"))
checkTrue(!(checkDataObject(testTinyResults, "Meta", "Post-Filter")))

testTinyFilter <- filterGenes(testTinyResults)
checkTrue(checkDataObject(testTinyFilter, "Meta", "Pre-Analysis"))
checkTrue(checkDataObject(testTinyFilter, "Meta", "Pre-Filter"))
checkTrue(checkDataObject(testTinyFilter, "Meta", "Post-Filter"))

checkTrue(!(checkDataObject(tinyMetaObject, "Dataset")))
checkTrue(!(checkDataObject(tinyMetaObject, "MetaAnalysis")))
checkTrue(!(checkDataObject(tinyMetaObject, "MetaFilter")))

checkTrue(checkDataObject(tinyMetaObject$originalData$PBMC.Study.1, "Dataset"))
checkTrue(checkDataObject(tinyMetaObject$metaAnalysis, "MetaAnalysis"))
checkTrue(checkDataObject(tinyMetaObject$filterResults[[1]], "MetaFilter"))

############
# Check when some of the pieces are null
############
tinyDataset <- tinyMetaObject$originalData[[1]]
checkTrue(checkDataObject(tinyDataset, "Dataset"))

tinyDataset$expr <- NULL
checkTrue(!(checkDataObject(tinyDataset, "Dataset")))

tinyDataset <- tinyMetaObject$originalData[[1]]
tinyDataset$keys <- NULL
checkTrue(!(checkDataObject(tinyDataset, "Dataset")))

tinyDataset <- tinyMetaObject$originalData[[1]]
tinyDataset$formattedName <- NULL
checkTrue(!(checkDataObject(tinyDataset, "Dataset")))

#Class is allowed to be null for generation of plots
tinyDataset <- tinyMetaObject$originalData[[1]]
tinyDataset$class <- NULL
checkTrue(checkDataObject(tinyDataset, "Dataset"))

tinyDataset <- tinyMetaObject$originalData[[1]]
tinyDataset$expr <- as.data.frame(tinyDataset$expr)
checkTrue(!(checkDataObject(tinyDataset, "Dataset")))

tinyDataset <- tinyMetaObject$originalData[[1]]
tinyDataset$keys <- as.matrix(tinyDataset$keys)
checkTrue(!(checkDataObject(tinyDataset, "Dataset")))

tinyDataset <- tinyMetaObject$originalData[[1]]
names(tinyDataset$class) <- NULL
checkTrue(!(checkDataObject(tinyDataset, "Dataset")))

###########
# Check some of the odd cases
##########

# Ensure checkDataObject returns false when all $keys are NULL
naKeys <- tinyMetaObject
naKeys$originalData$Whole.Blood.Study.1$keys <- rep(NA, length(tinyMetaObject$originalData$Whole.Blood.Study.1$keys))
names(naKeys$originalData$Whole.Blood.Study.1$keys) <- names(tinyMetaObject$originalData$Whole.Blood.Study.1$keys)
checkTrue(!(checkDataObject(naKeys, "Meta", "Pre-Analysis")))

# Ensure checkDataObject reutrns false when there is an infinite value in $expr
infExpr <- tinyMetaObject
infExpr$originalData$Whole.Blood.Study.1$expr[3,3] <- Inf
checkTrue(!(checkDataObject(infExpr, "Meta", "Pre-Analysis")))


