#' Subset samples for a particular dataset
#' @details 
#' Subsets all relevant slots within the Dataset object to include only the 
#' desired samples.
#' @param datasetObject the Dataset object to subset
#' @param keepMe either a binary vector for whether each sample should be 
#'    in the subset or a list of names of samples to be in the subset
#' @return returns a Dataset object that has been subsetted to the desired samples
#' @export
#' @author Winston A. Haynes
#' @examples 
#' subsetObject <- subsetOriginalData(tinyMetaObject$originalData$Whole.Blood.Study.1, 
#'    keepMe= c("Sample 1", "Sample 13", "Sample 43"))
subsetOriginalData <- function(datasetObject, keepMe) {
  if(length(keepMe) < length(datasetObject$class)) {
    if(sum(keepMe %in% names(datasetObject$class)) == 0) {
      stop("No specified sample names are present in the dataset object.")
    } 
    missingSamples <- setdiff(keepMe, names(datasetObject$class))
    if(length(missingSamples > 0)) {
      warning(paste("The following samples were specified, but not present in the data", 
                    paste(missingSamples, collapse = ", ")))
    }
  }
  datasetObject$pheno <- datasetObject$pheno[keepMe, ]
  datasetObject$expr <- datasetObject$expr[, keepMe]
  datasetObject$class <- datasetObject$class[keepMe]
  return(datasetObject)
}