#This file was written by Winn Haynes
#Last updated: March 30, 2015

library(RUnit)

#Load the RData object we are testing on
# load("data//tinyMetaObject.RData")

# source("R/forestPlot.R")

#Generate some test plots
#This output should be compared to the same file in "goldStandardOutput"
#pdf(file.path("..","MetaIntegrator","unitTests","testOutput","forestPlots.pdf"))

#############
#Should generate warning
#############

tt <- try(forestPlot(tinyMetaObject, geneName="12341546576877564534"))
checkTrue(is(tt,"try-error"))

####################
#Should not generate warning
####################
tt <- try(forestPlot(tinyMetaObject, geneName="Gene1"))
checkTrue(!is(tt,"try-error"))

tt <- try(forestPlot(tinyMetaObject, geneName="Gene2"))
checkTrue(!is(tt,"try-error"))

#dev.off()
