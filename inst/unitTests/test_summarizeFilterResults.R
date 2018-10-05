#######################################################################################
#metaIntegrator Unit Testing: summarizeFilterResult
#2015/04/20 11:03am @ Stanford
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#the purpose of this code is to run the unit testing for the summarizeFilterResult function
#######################################################################################

#load RUnit testing library
library(RUnit)
# 
# #Load the RData object we are testing on
# load("data/tinyMetaObject.RData")
# 
# #load runMetaAnalysis library
# source("R/summarizeFilterResults.R")

#######################################################################################
#Positive controls
#######################################################################################

#test
summarizePos <- try(summarizeFilterResults(tinyMetaObject,"pValueFDR0.05_es0_nStudies1_looaTRUE_hetero0"))
checkTrue(!is(summarizePos,"try-error"))

#######################################################################################
#Negative controls
#######################################################################################

#test
summarizeNeg <- try(summarizeFilterResults(tinyMetaObject,"BIGFAKELABEL"))
checkTrue(is(summarizeNeg,"try-error"))
