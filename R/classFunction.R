#' Helper function to build the class vector
#' @details 
#' Based on a defined set of disease terms, builds a class vector.
#' @param datasetObject the Dataset object to build a class vector for
#' @param column column from the $pheno slot to look for the disease terms
#' @param  diseaseTerms a list of terms identifying the disease samples
#' 
#' @return returns a Dataset object that has a class vector inserted
#' @export
#' @author Winston A. Haynes
#' @examples 
#' classObj <- classFunction(tinyMetaObject$originalData$Whole.Blood.Study.1, 
#'    column="group", diseaseTerms=c("Disease"))
classFunction<-function(datasetObject, column, diseaseTerms){
  currGroup<- trimws(as.character(datasetObject$pheno[,column]))
  
  #Healthy as 0
  classVec<-rep(0,length(currGroup))
  names(classVec) <- rownames(datasetObject$pheno)
  
  #Disease as 1
  classVec[currGroup %in% diseaseTerms] <- 1
  datasetObject$class <- classVec
  return(datasetObject)
}