#' Check for errors in objects used for analysis 
#' 
#' @description 
#' Given an object to check, its objectType and the objectStage, the function \code{checkDataObject} looks for errors within Meta, Dataset, MetaAnalyis, or MetaFilter objects. 
#' It returns TRUE if the object passed error checking, FALSE otherwise, and it prints warning messages explaining failed checks.
#' @usage checkDataObject(object, objectType, objectStage="")
#' @param object the object to be checked
#' @param objectType one type of "Meta", "Dataset", "MetaAnalysis", "MetaFilter"
#' @param objectStage if a \code{metaObject}, one of "Pre-Analysis", "Pre-Filter", or "Post-Filter". Otherwise: ""
#' @details For \code{metaAnalysisObject} and \code{filterObject}, it makes sure that each entry within the object is 1) not NULL and 2) the correct type.
#' For \code{datasetObjects}, it makes sure that:
#' 	1) the entries are not null (except \code{$class}, which is permitted to be NULL)  
#' 	2) the entries are the correct type and     
#' 	3) the sample names (within \code{$pheno}, \code{$expr}, and \code{$class}) match     
#' 	4) the probeIDs (within \code{$expr} and \code{$keys}) match.    
#' 	
#' 	For \code{metaObject}, it recursively checks the Dataset, MetaAnalysis, and MetaFilter objects contained within the \code{metaObject}. 
#'
#'The \code{objectStage} defines what entries a \code{metaObject} contains. Thus, "Pre-Analysis" \code{metaObjects} only contain \code{$originalData}. 
#'"Pre-Filter" \code{metaObjects} contain \code{$originalData}, \code{$metaAnalysis}, and \code{$leaveOneOutAnalysis}. 
#'"Post-Filter" \code{metaObjects} contain \code{$originalData}, \code{$metaAnalysis}, \code{$leaveOneOutAnalysis}, and \code{$filterResults}. 
#'
#' @return TRUE if passed error checking, FALSE otherwise
#' 	Prints warning messages explaining the portion of the error checking failed
#' @author Erika Bongen
#' @examples 
#'	# check a datasetObject
#'	checkDataObject(tinyMetaObject$originalData$Whole.Blood.Study.1, "Dataset")
#'	
#'	# check a metaObject before running the meta-analysis 
#'	checkDataObject(tinyMetaObject, "Meta", "Pre-Analysis")
#'	
#'	# check a metaObject after running the meta-analysis with runMetaAnalysis()
#'	checkDataObject(tinyMetaObject, "Meta", "Pre-Filter")
#'	
#'	# check a metaObject after filtering the meta-analysis results with filterGenes()
#'	checkDataObject(tinyMetaObject, "Meta", "Post-Filter")
#'	
#'	# check a metaAnalysisObject
#'	checkDataObject(tinyMetaObject$metaAnalysis, "MetaAnalysis")
#'
#'	# check a filterObject
#'	checkDataObject(tinyMetaObject$filterResults[[1]], "MetaFilter")
#' @keywords 
#'  utilities 
#' debugging 
#' @export
checkDataObject <- function(object, objectType, objectStage="") {
  if (objectType == "Meta") {
    result <- .metaCheckAll(object, objectStage)
  } else if (objectType == "Dataset") {
    result <- .datasetCheckAll(object)
  } else if (objectType == "MetaAnalysis") {
    result <- (.metaAnalysisCheckNull(object) && .metaAnalysisCheckType(object))
  } else if (objectType == "MetaFilter") {
    result <- .metaFilterCheckNull(object) && .metaFilterCheckType(object)
  } else {
    result <- FALSE
    warning("Invalid object type.")
  }
  
  return(result)
}
