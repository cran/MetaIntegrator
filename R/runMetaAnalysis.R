#' Run the meta-analysis algorithm 
#' @description
#' Given a \code{metaObject} with \code{$originalData} populated this function will run the meta-analysis algorithm.      
#' It returns a modified version of the \code{metaObject} with the meta-analysis results written into \code{metaObject$metaAnalysis} and the results of the leave-one-out analysis into \code{metaObject$leaveOneOutAnalysis}
#' @usage
#' 	runMetaAnalysis(metaObject, runLeaveOneOutAnalysis= TRUE, maxCores=Inf)
#' @param metaObject a metaObject which must have \code{metaObject$originalData} populated with a list of \code{datasetObjects} that will be used for discovery
#' @param runLeaveOneOutAnalysis TRUE to run leave one out analysis, FALSE otherwise (default: TRUE)
#' @param maxCores maximum number of cores to use during analysis (default: Inf)
#' @return modified version of the \code{metaObject} with \code{$metaAnalysis} and \code{$leaveOneOutAnalysis} populated
#' @author Francesco Vallania
#' @details
#' To make sure the input is correctly formatted, the input \code{metaObject} should be checked with \code{checkDataObject(metaObject, "Meta", "Pre-Analysis")} before starting the meta-analysis.
#' @seealso	\code{\link{checkDataObject}}
#' @examples
#'	#Run a meta analysis. 
#'	#		maxCores is set to 1 for package guideline compliance. 
#'	#		For personal purposes, leave parameter un-set.
#'	runMetaAnalysis(tinyMetaObject, maxCores=1)
#'@keywords 
#' methods 
#' @export
#' @import ggplot2 parallel data.table
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
  es_plot <- ggplot(reshape2::melt(metaObject$metaAnalysis$datasetEffectSizes, varnames=c("Gene", "Study")),
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
