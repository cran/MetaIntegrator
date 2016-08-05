#This file was written by Winn Haynes
#Last updated: April 20, 2015

library(RUnit)
# source("R/calculateScore.R")
# load("data/tinyMetaObject.RData")

#########
#These should pass
#########
scoreResults <- calculateScore(tinyMetaObject$filterResults[[1]], tinyMetaObject$originalData[[1]]) 
checkEquals(scoreResults[[1]], -0.03687, tolerance=1e-3)
checkEquals(scoreResults[[45]], 2.12351, tolerance=1e-3)
