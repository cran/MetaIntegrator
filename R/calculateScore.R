calculateScore <- function(filterObject, datasetObject, suppressMessages=FALSE) {
	#Check that our objects are the right type
	if(! checkDataObject(object = filterObject, objectType="MetaFilter")) {
		stop("filterObject that was passed to calculateScore was not appropriately formatted as a MetaFilter object")
	} 
	if(!checkDataObject(object =datasetObject, objectType =  "Dataset")){
		stop("datasetObject that was passed to calculateScore was not appropriately formatted as a Dataset object")
	}
	
	#Get the positive and negative genes
	pos.genes <- filterObject$posGeneNames
	neg.genes <- filterObject$negGeneNames
	
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
	  cat("Used ", dim(posGenes)[1],"of ", length(pos.genes)," pos genes, and ",
			dim(negGenes)[1] ," of ", length(neg.genes)," neg genes \n")
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
