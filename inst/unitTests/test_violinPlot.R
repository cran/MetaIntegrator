#This file was written by Winn Haynes
#Last updated: March 30, 2015

library(RUnit)

# #Load the RData object we are testing on
# load("data/tinyMetaObject.RData")
# 
# source("R/violinPlot.R")

#Generate some test plots
#This output should be compared to the same file in "goldStandardOutput"
#pdf(file.path("..","MetaIntegrator","unitTests","testOutput","violinPlots.pdf"))

#################
#These plots should all work
#################

tt <- try(violinPlot(tinyMetaObject$filterResults$pValueFDR0.05_es0_nStudies1_looaTRUE_hetero0, 
										 tinyMetaObject$originalData$Whole.Blood.Study.1, labelColumn="group"))
checkTrue(!is(tt,"try-error"))

tt <- try(violinPlot(tinyMetaObject$filterResults$pValueFDR0.05_es0_nStudies1_looaTRUE_hetero0, 
										 tinyMetaObject$originalData$Whole.Blood.Study.2, labelColumn="group"))
checkTrue(!is(tt,"try-error"))

tt <- try(violinPlot(tinyMetaObject$filterResults$pValueFDR0.05_es0_nStudies1_looaTRUE_hetero0, 
										 tinyMetaObject$originalData$PBMC.Study.1, labelColumn="group"))
checkTrue(!is(tt,"try-error"))

##################
# None of these plots should work
##################
tt <- try(violinPlot(tinyMetaObject$filterResults$pValueFDR0.05_es0_nStudies1_looaTRUE_hetero0, 
										 tinyMetaObject$originalData$Whole.Blood.Study.1, labelColumn="fake"))
checkTrue(is(tt,"try-error"))

#dev.off()
