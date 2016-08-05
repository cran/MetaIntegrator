#This file was written by Winn Haynes
#Last updated: April 8, 2015

library(RUnit)

################
# Test on simulated data
################

set.seed(1)
labels <- c(rep(0, 500), rep(1, 500))
scores <- runif(1000)
rocRes <- calculateROC(labels, scores)

#AUC should be near 0.5, but only want to fail test if way off
checkEquals(rocRes$auc[[1]], 0.507664, tolerance= 1e-3)

################
# Test on our data
################
scoreResults <- calculateScore(tinyMetaObject$filterResults[[1]], tinyMetaObject$originalData[[1]]) 
rocRes <- calculateROC(predictions=scoreResults, labels=tinyMetaObject$originalData[[1]]$class)
checkEquals(rocRes$auc[[1]], 0.85606, tolerance= 1e-3)

scoreResults <- calculateScore(tinyMetaObject$filterResults[[1]], tinyMetaObject$originalData[[2]]) 
rocRes <- calculateROC(predictions=scoreResults, labels=tinyMetaObject$originalData[[2]]$class)
checkEquals(rocRes$auc[[1]], 0.86738, tolerance= 1e-3)

scoreResults <- calculateScore(tinyMetaObject$filterResults[[1]], tinyMetaObject$originalData[[3]]) 
rocRes <- calculateROC(predictions=scoreResults, labels=tinyMetaObject$originalData[[3]]$class)
checkEquals(rocRes$auc[[1]], 0.95469, tolerance= 1e-3)
