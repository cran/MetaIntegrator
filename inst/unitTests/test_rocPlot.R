#This file was written by Winn Haynes
#Last updated: April 8, 2015

library(RUnit)

#Load the RData object we are testing on
# load("data/tinyMetaObject.RData")
# 
# source("R/ROC_and_search_functions.R")

#Generate some test plots
#This output should be compared to the same file in "goldStandardOutput"
#pdf(file.path("..","MetaIntegrator","unitTests","testOutput","rocPlots.pdf"))

#################
#These plots should all work
#################
rocPlot(tinyMetaObject$filterResults[[1]], tinyMetaObject$originalData[[1]])
tt <- try(rocPlot(tinyMetaObject$filterResults[[1]], tinyMetaObject$originalData[[1]]))
checkTrue(!is(tt,"try-error"))

rocPlot(tinyMetaObject$filterResults[[1]], tinyMetaObject$originalData[[2]])
tt <- try(rocPlot(tinyMetaObject$filterResults[[1]], tinyMetaObject$originalData[[2]]))
checkTrue(!is(tt,"try-error"))

rocPlot(tinyMetaObject$filterResults[[1]], tinyMetaObject$originalData[[3]])
tt <- try(rocPlot(tinyMetaObject$filterResults[[1]], tinyMetaObject$originalData[[3]]))
checkTrue(!is(tt,"try-error"))

#dev.off()
