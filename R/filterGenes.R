#' Filter out significant genes from meta-analysis results
#' 
#' @description After the Meta-Analysis results have been written to the \code{metaObject}, 
#' the results can be examined using different gene filtering criteria. 
#' This function will use the given filterParameter to select genes that fulfill the filter conditions. 
#' The function returns a modified version of the \code{metaObject} with results stored in \code{metaObject$filterResults}
#' 
#' @param metaObject a Meta object which must have the \code{$originalData}, \code{$metaAnalysis} populated
#' @param isLeaveOneOut Do leave-one-out analysis on discovery datasets (default: TRUE). Needs at least 2 datasets for discovery.
#' @param FDRThresh FDR cutoff: a gene is selected, if it has a p-value less than or equal to the FDR cutoff (default: 0.05)
#' @param effectSizeThresh a gene is selected, if the absolute value of its effect size is above this threshold (default: 0)
#' @param numberStudiesThresh number of studies in which a selected gene has to be significantly up/down regulated (default: 1)
#' @param heterogeneityPvalThresh heterogeneity p-value cutoff (filter is off by default: \code{heterogeneityPvalThresh = 0}). 
#' Genes with significant heterogeneity and, thus a significant (low) heterogeneity p-value, can be filtered out by using e.g.: 
#' \code{heterogeneityPvalThresh = 0.05} (removes all genes with heterogeneity p-value < 0.05)
#' 
#' @return A modified version of the input metaObject with an additional filterObject stored within \code{metaObject$filterResults}
#' 
#' @author Francesco Vallania
#' @note Use \code{checkDataObject(metaObject, "Meta", "Pre-Filter")} to make sure your metaObject has the right format for filtering after running the meta-analysis with \code{runMetaAnalysis().}
#' @seealso \code{\link{checkDataObject}}
#' @examples
#' 	# filter genes with default settings 
#' 	#(false discovery rate cutoff of 5 percent and WITH leave-one-out analysis)
#' 	testMetaObject <- filterGenes(tinyMetaObject)
#' 	summarizeFilterResults(testMetaObject, getMostRecentFilter(testMetaObject))
#' 	
#' 	# filter genes with false discovery rate of 1 percent and WITHOUT leave-one-out analysis  
#' 	testMetaObject <- filterGenes(testMetaObject, FDRThresh = 0.01, isLeaveOneOut = FALSE)
#' 	summarizeFilterResults(testMetaObject, getMostRecentFilter(testMetaObject))
#' @keywords 
#' methods 
#'  classif 
#' @export
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
