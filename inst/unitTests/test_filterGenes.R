#######################################################################################
#metaIntegrator Unit Testing: filterGenes
#2015/04/20 10:54am @ Stanford
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#the purpose of this code is to run the unit testing for the filterGenes function
#######################################################################################

#load RUnit testing library
library(RUnit)

# #Load the RData object we are testing on
# load("data/tinyMetaObject.RData")
# 
# #load runMetaAnalysis library
# source("R/filterGenes.R")

#######################################################################################
#Positive controls
#######################################################################################

#test without LeaveOneOut analysis
filterNoLOOA <- try(filterGenes(tinyMetaObject,isLeaveOneOut = FALSE))
checkTrue(!is(filterNoLOOA,"try-error"))

#test without LeaveOneOut analysis
filterLOOA   <- try(filterGenes(tinyMetaObject,isLeaveOneOut = TRUE))
checkTrue(!is(filterLOOA,"try-error"))

#######################################################################################
#Negative controls
#######################################################################################

#Only one case in class
wrongTinyMetaObject              <- tinyMetaObject
wrongTinyMetaObject$metaAnalysis <- NULL

#test without LeaveOneOut analysis
wrongFilterNoLOOA <- try(filterGenes(wrongTinyMetaObject,isLeaveOneOut = FALSE))
checkTrue(is(wrongFilterNoLOOA,"try-error"))

#test without LeaveOneOut analysis
wrongFilterLOOA   <- try(filterGenes(wrongTinyMetaObject,isLeaveOneOut = TRUE))
checkTrue(is(wrongFilterLOOA,"try-error"))


