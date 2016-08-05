#######################################################################################
#metaIntegrator Functions: runMetaAnalysis
#2015/03/02 3:48pm @ Stanford
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#the purpose of this code is to contain the runMetaAnalysis function
#######################################################################################

#######################################################################################
#runMetaAnalysis
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#this function runs the metanalysis pipeline from metaObject and returns a new
#metaObject that contains results from metaAnalysis and leaveOneOutAnalysis
#######################################################################################
runMetaAnalysis <- function(metaObject,
                            runLeaveOneOutAnalysis = TRUE,
                            maxCores=Inf){
	old=FALSE
  
  #insert check for names here:-> take care of this with Erika
  metaObject <- .originalDataNameConverter(metaObject)
  
  #Run checkDataObject function here to check for the object
  myObjectCheck <- checkDataObject(metaObject,'Meta','Pre-Analysis')
  
  #in case something did not pass just 
  if(myObjectCheck==FALSE){
    stop("Error in the input object!")
  }
  
  #Insert check here to ensure metaObject has at least
  #two cases and two controls
  classCheck <- sapply(metaObject$originalData,function(i) .classVectorChecker(i))
  
  #if this does not pass the check then crash
  if(all(classCheck)==FALSE){
    warningString <- paste("In ",names(which(classCheck!=TRUE))," there are not @ least two cases and/or two controls")
    stop(warningString)
  }
  
  #run runMetaAnalysisCore function on the all dataset and assign leaveOneOutAnalysis 
  #to a default value of NULL [rename function]
  metaObject$metaAnalysis        <- .runMetaAnalysisCore(metaObject$originalData,
                                                         old)
  metaObject$leaveOneOutAnalysis <- NULL
  
  #plot the ES histograms here
  es_plot <- ggplot(melt(metaObject$metaAnalysis$datasetEffectSizes, varnames=c("Gene", "Study")),
                    aes_string(x      = "value",
                        colour = "Study")) + 
    geom_density(size = 1.1)            + 
    theme_bw()                          + 
    scale_color_discrete(name = 'Dataset')
  
  #plot the ES histograms here [for the user to be collected]
  suppressWarnings(print(es_plot))
  
  #run leaveOneOutAnalysis if wanted and if possible 
  #[# of datasets > 2] and add tag to descriptor file
  if(length(metaObject$originalData) > 2 && runLeaveOneOutAnalysis == TRUE){
    #run meta-analysis for LOOA
    metaObject$leaveOneOutAnalysis <- .leaveOneOutMetaAnalysisWrapper(metaObject$originalData,old, maxCores=maxCores)
    
    #name objects for LOOA
    names(metaObject$leaveOneOutAnalysis) <- paste("removed",names(metaObject$originalData),sep = "_")
  }
  
  #return joint outputs in a single list
  return(metaObject)
}
