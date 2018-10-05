#This file was written by Winn Haynes
#Last updated: April 8, 2015

library(RUnit)

# source("R/getSampleLevelGeneData.R")

#########
#These should pass
#########
sampleResults <- getSampleLevelGeneData(datasetObject=tinyMetaObject$originalData[[1]], 
																				geneNames=c(tinyMetaObject$filterResults[[1]]$posGeneNames, tinyMetaObject$filterResults[[1]]$negGeneNames))
checkEquals(sampleResults["Gene27", "Sample.1"], 7.016342, tolerance=1e-3)
checkEquals(sampleResults["Gene42", "Sample.49"], 9.196237, tolerance=1e-3)
