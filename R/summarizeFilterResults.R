#' Summarize the filtered analysis results
#'@description
#'	Given a  \code{metaObject} and the name of the \code{filterObject} of interest, this function will print a summary style message about genes that passed the filtering step using the function \code{filterGenes()} and return a \code{dataFrame} that contains the \code{$pooledResults} information for each gene which passed the filter.
#'@usage
#'	summarizeFilterResults(metaObject, metaFilterLabel)
#'@param metaObject the metaObject that contains the \code{filterObject} of interest
#'@param metaFilterLabel the name of a \code{filterObject} generated with the function \code{filterGenes()}
#'@return Data frame, which contains \code{$pooledResults} information for each gene which passed the filter
#'@author Francesco Vallania
#'@seealso
#'	\code{\link{filterGenes}}
#'@examples
#'	# filter genes with default settings 
#'	#		false discovery rate cutoff of 5 percent and WITH leave-one-out analysis
#'	testMetaObject <- filterGenes(tinyMetaObject)
#'	summarizeFilterResults(testMetaObject, getMostRecentFilter(testMetaObject))
#'@keywords
#' utilities
#' methods
#' @export
summarizeFilterResults <- function(metaObject,metaFilterLabel){
  
  #check metaObject here at this step
  myObjectCheck <- checkDataObject(metaObject,'Meta','Post-Filter')
  
  #in case something did not pass just 
  if(myObjectCheck==FALSE){
    stop("Error in the input object!")
  }
  
  #in case label is missing
  if(is.null(metaObject$filterResults[[metaFilterLabel]])==TRUE){
    stop("Error! The input label does not match any filter run")
  }
  
  #extract genes from filteredObject
  pos_genes <- metaObject$filterResults[[metaFilterLabel]]$posGeneNames
  neg_genes <- metaObject$filterResults[[metaFilterLabel]]$negGeneNames
  
  #select genes
  posDF <- metaObject$metaAnalysis$pooledResults[match(pos_genes,row.names(metaObject$metaAnalysis$pooledResults)),]
  negDF <- metaObject$metaAnalysis$pooledResults[match(neg_genes,row.names(metaObject$metaAnalysis$pooledResults)),]
  
  #if empty set to NULL
  #re-sort the data.frames to rank genes by their FDRs
  if(nrow(posDF)==0){
    posDF <- NULL
  }else{
    posDF <- posDF[order(posDF$effectSizeFDR),]
  }
  if(nrow(negDF)==0){
    negDF <- NULL
  }else{
    negDF <- negDF[order(negDF$effectSizeFDR),]
  }
  
  #return data frames of positive and nagative genes
  return(list(pos=posDF,
              neg=negDF))
}
