#######################################################################################
#metaIntegrator Functions: filterGenes
#2015/03/02 3:48pm @ Stanford
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#the purpose of this code is to contain the filterGenes function
#######################################################################################

#######################################################################################
#filterGenes
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#this function filters a metaAnalysis run on a metaObject
#######################################################################################
filterGenes <- function(metaObject,
                        isLeaveOneOut           = TRUE,
                        effectSizeThresh        = 0,
                        FDRThresh               = 0.05,
                        numberStudiesThresh     = 1,
                        heterogeneityPvalThresh = 0){
  
  #check metaObject here at this step
  myObjectCheck <- checkDataObject(metaObject,'Meta','Pre-Filter')
  
  #in case something did not pass just 
  if(myObjectCheck==FALSE){
    stop("Error in the input object!")
  }
  
  #filter metaObject
  meta_filter_out <- NULL
  
  #make leaveOneOut FALSE if it hasn't been run
  if(is.null(metaObject$leaveOneOutAnalysis)){
    isLeaveOneOut <- FALSE
  }
  
  #if this is leaveOneOut and there are more than 2 objects then go for it
  if(isLeaveOneOut == TRUE & length(metaObject$originalData) > 2){
    
    #iterate over all the leve-one-out metafilters
    loout_mf <- lapply(metaObject$leaveOneOutAnalysis,
                       function(i)
                         .filterMetaRun(i,
                                        effectSizeThresh        = effectSizeThresh,
                                        effectFDRSizeThresh     = FDRThresh,
                                        fisherFDRThresh         = FDRThresh,
                                        numberStudiesThresh     = numberStudiesThresh-1,
                                        heterogeneityPvalThresh = heterogeneityPvalThresh))
    
    #create final list of genes
    posGeneNames <- loout_mf[[1]]$posGeneNames
    negGeneNames <- loout_mf[[1]]$negGeneNames
    
    #perform sequential intersections
    for(i in 2:length(loout_mf)){
      #get only the common genes
      posGeneNames <- intersect(posGeneNames,loout_mf[[i]]$posGeneNames)
      negGeneNames <- intersect(negGeneNames,loout_mf[[i]]$negGeneNames)
    }
    
    #compose final meta_filter_out list
    meta_filter_out <-  list(posGeneNames            = posGeneNames,
                             negGeneNames            = negGeneNames,
                             effectSizeThresh        = effectSizeThresh,
                             FDRThresh               = FDRThresh,
                             numberStudiesThresh     = numberStudiesThresh,
                             heterogeneityPvalThresh = heterogeneityPvalThresh)
    
    #add labels for LOO
    meta_filter_out$isLeaveOneOut <- TRUE    
  }else{
    
    #just run vanilla metaFiltering
    meta_filter_out <- .filterMetaRun(metaObject$metaAnalysis,
                                      effectSizeThresh        = effectSizeThresh,
                                      effectFDRSizeThresh     = FDRThresh,
                                      fisherFDRThresh         = FDRThresh,
                                      numberStudiesThresh     = numberStudiesThresh,
                                      heterogeneityPvalThresh = heterogeneityPvalThresh)
    
    #add labels for LOO
    meta_filter_out$isLeaveOneOut <- FALSE
  }
  
  #get filter description
  meta_filter_out$filterDescription <- "MetaFilter from version 1.0"
  
  #get systems time
  meta_filter_out$timestamp         <- Sys.time()
  
  #compose filterName
  filter_name <- paste("FDR",      FDRThresh,
                       "_es",      effectSizeThresh,
                       "_nStudies",numberStudiesThresh,
                       "_looa",    isLeaveOneOut,
                       "_hetero",  heterogeneityPvalThresh,
                       sep = "")
  
  #add it to metaObject
  metaObject$filterResults[[filter_name]] <- meta_filter_out
  
  #return metaObject  
  return(metaObject)
}
