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

forwardRes <- forwardSearch(tinyMetaObject, tinyMetaObject$filterResults[[1]], forwardThresh = 0) 
checkEquals(length(forwardRes$posGeneNames), 1)
checkEquals(length(forwardRes$negGeneNames), 1)
checkEquals(forwardRes$posGeneNames[[1]], "Gene13")
checkEquals(forwardRes$negGeneNames[[1]], "Gene45")
