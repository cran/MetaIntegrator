## ----import package, message=FALSE, echo=FALSE, warning=FALSE------------
#setwd("../")
library(MetaIntegrator) # <- import library

# set to TRUE if R code should be executed 
eval=TRUE

## ----style, echo = FALSE, results = 'asis', warning=FALSE, message=FALSE----
BiocStyle::markdown()

## ----env, message=FALSE, echo=FALSE--------------------------------------
# Biocpkg("IRanges")

## ----example dsobj, eval=eval--------------------------------------------
dataObj1 <- tinyMetaObject$originalData$PBMC.Study.1
str(dataObj1, max.level = 1)

## ----checkExpression, eval=eval, fig.height=4, fig.width=10--------------
boxplot(dataObj1$expr[,1:15]) # -> shows samples 1-15, to see all run: boxplot(dataObj1$expr) 

## ----check1, eval=eval---------------------------------------------------
checkDataObject(dataObj1, "Dataset")

## ----example disc_data, eval=eval----------------------------------------
# use the additional 2 example datasets from tinyMetaObject
dataObj2 = tinyMetaObject$originalData$Whole.Blood.Study.1
dataObj3 = tinyMetaObject$originalData$Whole.Blood.Study.2
# and create the metaObject
discovery_datasets <- list(dataObj1, dataObj2, dataObj3)
names(discovery_datasets) = c(dataObj1$formattedName, dataObj2$formattedName, dataObj3$formattedName)
exampleMetaObj=list() 
exampleMetaObj$originalData <- discovery_datasets

## ----check2, eval=TRUE---------------------------------------------------
checkDataObject(exampleMetaObj, "Meta", "Pre-Analysis")

## ----runMetaAnalysis1, eval=eval, message=FALSE, warning=FALSE-----------
exampleMetaObj <- runMetaAnalysis(exampleMetaObj, maxCores=1)

## ----runMetaAnalysis2, eval=eval, message=FALSE, warning=FALSE-----------
str(exampleMetaObj, max.level = 2)

## ----filterGenes, eval=eval----------------------------------------------
exampleMetaObj <- filterGenes(exampleMetaObj, isLeaveOneOut = TRUE, FDRThresh = 0.001)

## ----summarizeFilterResults1, eval=FALSE---------------------------------
#  summarizeFilterResults(exampleMetaObj, "FDR0.001_es0_nStudies1_looaTRUE_hetero0")

## ----summarizeFilterResults2, eval=eval----------------------------------
summarizeFilterResults(exampleMetaObj, getMostRecentFilter(exampleMetaObj))

## ----violinPlot, message=FALSE, warning=FALSE, eval=eval, fig.height=8, fig.width=8----
violinPlot(exampleMetaObj$filterResults$FDR0.001_es0_nStudies1_looaTRUE_hetero0, dataObj2, labelColumn = 'group')

## ----rocPlot, eval=eval, fig.height=8, fig.width=8-----------------------
rocPlot(exampleMetaObj$filterResults$FDR0.001_es0_nStudies1_looaTRUE_hetero0, dataObj2, title = "ROC plot for discovery dataset2, FDR: 0.001")

## ----forestPlot, eval=eval, fig.height=5, fig.width=7--------------------
forestPlot(exampleMetaObj, "Gene27")

## ----calculateScore, eval=eval-------------------------------------------
calculateScore(exampleMetaObj$filterResults$FDR0.001_es0_nStudies1_looaTRUE_hetero0, dataObj2)

## ---- eval=FALSE---------------------------------------------------------
#  metaObject <- runMetaAnalysis(metaObject)

## ---- eval=FALSE---------------------------------------------------------
#  metaObject <- filterGenes(metaObject, filterParameter)

## ---- eval=FALSE---------------------------------------------------------
#  summarizeFilterResults(metaObject, metaFilterLabel)

## ---- eval=FALSE---------------------------------------------------------
#  calculateScore(datasetObject, filterObject)

## ---- eval=FALSE---------------------------------------------------------
#  forestPlot(metaObject, geneName)

## ---- eval=FALSE---------------------------------------------------------
#  violinPlot(filterObject, datasetObject, labelColumn)

## ---- eval=FALSE---------------------------------------------------------
#  rocPlot(filterObject, datasetObject, title = "ROC Plot")

## ---- eval=FALSE---------------------------------------------------------
#  forwardSearch(metaObject, geneList, yes.pos = NULL, yes.neg = NULL, forwardThresh = 0)

## ---- eval=FALSE---------------------------------------------------------
#  #Run a forward search
#  forwardRes <- forwardSearch(tinyMetaObject, tinyMetaObject$filterResults[[1]], forwardThresh = 0)

## ---- eval=FALSE---------------------------------------------------------
#  backwardSearch(metaObject, geneList, backThresh = 0)

## ---- eval=FALSE---------------------------------------------------------
#  #Run a backward search
#  backwardRes <- backwardSearch(tinyMetaObject, tinyMetaObject$filterResults[[1]], backThresh = -3)

## ---- eval=FALSE---------------------------------------------------------
#  checkDataObject(object, objectType, objectStage)

## ---- eval=FALSE---------------------------------------------------------
#  # check a datasetObject
#  checkDataObject(tinyMetaObject$originalData$Whole.Blood.Study.1, "Dataset")
#  
#  # check a metaObject before running the meta-analysis
#  checkDataObject(tinyMetaObject, "Meta", "Pre-Analysis")
#  
#  # check a metaObject after running the meta-analysis with runMetaAnalysis()
#  checkDataObject(tinyMetaObject, "Meta", "Pre-Filter")
#  
#  # check a metaObject after filtering the meta-analysis results with filterGenes()
#  checkDataObject(tinyMetaObject, "Meta", "Post-Filter")
#  
#  # check a metaAnalysisObject
#  checkDataObject(tinyMetaObject$metaAnalysis, "MetaAnalysis")
#  
#  # check a filterObject
#  checkDataObject(tinyMetaObject$filterResults[[1]], "MetaFilter")

## ---- eval=FALSE---------------------------------------------------------
#  getMostRecentFilter(metaObject)

## ---- eval=FALSE---------------------------------------------------------
#  calculateROC(labels, predictions, AUConly = F)

## ---- eval=FALSE---------------------------------------------------------
#  getSampleLevelGeneData(datasetObject, geneNames)

