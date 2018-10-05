#' Calculate a signature Z-score for a set of genes in a single dataset
#' 

#' @param filterObject	a MetaFilter object generated with \code{filterGenes()} containing the signature genes that will be used for Z-score calculation. 
#' @param datasetObject A Dataset object for which the signature score (Z-score) will be calculated. This vector would typically be added as \code{$score} column in \code{datasetObject$pheno}.
#' @param suppressMessages Boolean value (TRUE/FALSE) about whether to display verbose output. Default: FALSE.
#'
#' @return 	A vector of Z-scores, of length \code{ncols(datasetObject$expr)} (and in the same order).
#' 	
#' @description 	
#' Given a gene set of interest, it is often desirable to summarize the expression of that gene set using a single integrated score.
#' 	 The \code{calculateScore} method calculates the geometric mean of the expression level of all positive genes, 
#' 	 minus the geometric mean of the expression level of all negative genes. The resulting scores are then standardized within the given dataset, such that the output Z-score has mean=0 and std. dev=1.
#' 	  Such a Z-score can then be used for classification, etc. 
#' @usage calculateScore(filterObject, datasetObject, suppressMessages=FALSE)
#' 	 	  
#' @details 
#' 		The Z-score is based off of the geometric mean of expression. As such, negative expression values are not allowed. A dataset is thus always scaled by its minimum value + 1, such that the lowest value = 1. Any individual NANs or NAs are also set to 1. If a dataset does not have any information on a given gene, the entire gene is simply left out of the score. When run, the function will print to command line the number of genes used, and the number passed in. 
#' 		Although mostly used internally, the function has been exported in case users want to compare multiple classes, etc., using the same Z-score as is used for producing two-class comparisons. 
#' @author Timothy E. Sweeney, Winston A. Haynes
#' @seealso 	\code{\link{filterGenes}}
#' @export
#' @examples calculateScore(tinyMetaObject$filterResults[[1]], tinyMetaObject$originalData[[1]]) 

calculateScore <- function(filterObject, datasetObject, suppressMessages=FALSE) {
	#Check that our objects are the right type
	if(! checkDataObject(object = filterObject, objectType="MetaFilter")) {
		stop("filterObject that was passed to calculateScore was not appropriately formatted as a MetaFilter object")
	} 
	if(!checkDataObject(object =datasetObject, objectType =  "Dataset")){
		stop("datasetObject that was passed to calculateScore was not appropriately formatted as a Dataset object")
	}
	
	#Get the positive and negative genes, remove the empty gene name
	pos.genes <- setdiff(filterObject$posGeneNames,c("",NA))
	neg.genes <- setdiff(filterObject$negGeneNames,c("",NA))
	
	## Adjust Expression values so that they are all positive
	datasetObjectmin <- min(datasetObject$expr, na.rm=TRUE)
	if (datasetObjectmin <0) {datasetObject$expr <- datasetObject$expr + abs(datasetObjectmin) + 1}
	
	#Calculate the score of the positive genes from our signature
	posScore <- 0; posGenes <- NULL;
	if(sum(pos.genes %in% datasetObject$keys) > 0){
		posGenes <- getSampleLevelGeneData(datasetObject, pos.genes)
		posScore <- apply(posGenes, 2, .geomMean)
	}  
	
	#Calculate the score of the negative genes from our signature
	negScore <- 0; negGenes <- NULL;
	if(sum(neg.genes %in% datasetObject$keys) > 0){
		negGenes <- getSampleLevelGeneData(datasetObject, neg.genes)
		negScore <- apply(negGenes, 2, .geomMean)  
	}
	
	
	if(is.null(negGenes) & is.null(posGenes)) {
    cat("Used 0 of ", length(pos.genes)," pos genes, and 0 of ", length(neg.genes)," neg genes \n")
    return(rep(0, ncol(datasetObject$expr)))
	} 
  if(is.null(negGenes)){
		totalScore <- scale(posScore)
	} else { 
		if(is.null(posGenes)) {
			totalScore <- scale(-negScore)
		} else {
			totalScore <- scale(posScore -negScore)
		}
	}
	
  if(!suppressMessages) {
	  cat("Used ", ifelse(is.null(dim(posGenes)[1]),0,dim(posGenes)[1]),"of ", length(pos.genes)," pos genes, and ",
	      ifelse(is.null(dim(negGenes)[1]),0,dim(negGenes)[1]) ," of ", length(neg.genes)," neg genes \n")
  }
	return(as.vector(totalScore))
}

#Calculated geometric mean from a vector
#Args:
#	x: numeric vector
#Returns:
#	number
.geomMean <- function (x) {
	if (!is.numeric(x) && !is.complex(x) && !is.logical(x)) {
		stop("argument is not numeric or logical: returning NA")
	}
	if (any(x < 0)) {
		stop("'x' contains negative value(s)")
	}
	return(exp(sum(log(x))/length(x)))
}
