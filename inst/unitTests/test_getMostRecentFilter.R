#######################################################################################
#metaIntegrator Unit Testing: getMostRecentFilter
#2015/04/20 11:00am @ Stanford
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#the purpose of this code is to run the unit testing for the getMostRecentFilter function
#######################################################################################

#load RUnit testing library
library(RUnit)

# #Load the RData object we are testing on
# load("data/tinyMetaObject.RData")
# 
# #load runMetaAnalysis library
# source("R/getMostRecentFilter.R")

#######################################################################################
#Positive controls
#######################################################################################

#test
getMostRFOut <- try(getMostRecentFilter(tinyMetaObject))
checkTrue(!is(getMostRFOut,"try-error"))

#######################################################################################
#Negative controls
#######################################################################################

#Only one case in class
wrongTinyMetaObject               <- tinyMetaObject
wrongTinyMetaObject$filterResults <- NULL

#test
wrongGetMostRFOut <- try(getMostRecentFilter(wrongTinyMetaObject))
checkTrue(is(wrongGetMostRFOut,"try-error"))
