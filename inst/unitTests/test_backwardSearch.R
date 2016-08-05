#This file was written by Winn Haynes
#Last updated: April 8, 2015

library(RUnit)

#Load the RData object we are testing on
# load("data/tinyMetaObject.RData")
# 
# source("R/ROC_and_search_functions.R")


#########
#These should pass
#########

backwardRes <- backwardSearch(tinyMetaObject, tinyMetaObject$filterResults[[1]], backThresh = -3) 
checkEquals(length(backwardRes$posGeneNames), 12)
checkEquals(length(backwardRes$negGeneNames), 10)

#XAF1 is in the old results...
checkTrue("Gene56" %in% tinyMetaObject$filterResults[[1]]$negGeneNames)

#but has been removed from the new results
checkTrue(!("Gene56" %in% backwardRes$negGeneNames[[1]]))
