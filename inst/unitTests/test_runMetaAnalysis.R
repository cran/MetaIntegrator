#######################################################################################
#metaIntegrator Unit Testing: runMetaAnalysis
#2015/04/20 10:38am @ Stanford
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#the purpose of this code is to run the unit testing for the runMetaAnalysis function
#######################################################################################

#load RUnit testing library
library(RUnit)

# #Load the RData object we are testing on
# load("data/tinyMetaObject.RData")
# 
# #load runMetaAnalysis library
# source("R/runMetaAnalysis.R")

#######################################################################################
#Positive controls
#######################################################################################

#test without LeaveOneOut analysis
runMetaNoLOOA <- try(runMetaAnalysis(tinyMetaObject,runLeaveOneOutAnalysis = FALSE, maxCores = 1))
checkTrue(!is(runMetaNoLOOA,"try-error"))

#test with LeaveOneOut analysis
runMetaLOOA <- try(runMetaAnalysis(tinyMetaObject,runLeaveOneOutAnalysis = TRUE, maxCores = 1))
checkTrue(!is(runMetaLOOA,"try-error"))

#######################################################################################
#Negative controls
#######################################################################################

#Only one case in class
wrongTinyMetaObject <- tinyMetaObject
wrongTinyMetaObject$originalData$Whole.Blood.Study.1$class[1:length(wrongTinyMetaObject$originalData$Whole.Blood.Study.1$class)] <- 0
wrongTinyMetaObject$originalData$Whole.Blood.Study.1$class[1]<- 1

#test without LeaveOneOut analysis
wrongRunMetaNoLOOA <- try(runMetaAnalysis(wrongTinyMetaObject,runLeaveOneOutAnalysis = FALSE))
checkTrue(is(wrongRunMetaNoLOOA,"try-error"))

#test with LeaveOneOut analysis 
wrongRunMetaLOOA   <- try(runMetaAnalysis(wrongTinyMetaObject,runLeaveOneOutAnalysis = TRUE))
checkTrue(is(wrongRunMetaLOOA,"try-error"))
