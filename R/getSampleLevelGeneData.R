#' Extract gene-level data from a given data object
#' 
#' @param  datasetObject a Dataset object that is used to extract sample level data (At least, must have a \code{$expr} of probe-level data, and \code{$keys} of probe:gene mappings).
#' @param geneNames		A vector of geneNames
#' 	
#' @return Returns a data frame with expression levels of only the genes of interest, for each sample in the dataset. 
#' 	Mostly used internally, but has been exposed to the user to allow advanced functionality on external datasets if desired. 
#' 	
#' @description 
#' 	Given a \code{datsetObject}, and a set of target genes, this function will summarize probe-level data to gene-level data for the target genes. Returns a data frame with only the genes of interest, for each sample in the dataset. 
#' 	
#' @details Summarizes probe-level data to gene-level data, using the mean of the probes, according to the probe:gene mapping in the \code{$keys} item in the dataset object. This is done only for the genes in the filter object.
#' @author Timothy E. Sweeney, Winston A. Haynes
#' @export
#' @examples 
#'	sampleResults <- getSampleLevelGeneData(datasetObject=tinyMetaObject$originalData[[1]], 
#'	geneNames=c(tinyMetaObject$filterResults[[1]]$posGeneNames, 
#'	  tinyMetaObject$filterResults[[1]]$negGeneNames))

getSampleLevelGeneData <- function(datasetObject, geneNames){
	#Subfunctions all contained within for cleanliness
	GenesMtx <- .extractDataFromGEM(datasetObject, geneNames)
	GenesMtx <- .replaceValues(GenesMtx, 0, 1)
	GenesMtx <- .replaceNaNs(GenesMtx, 1)
	GenesMtx <- .replaceNAs(GenesMtx, 1)
	return(GenesMtx)
}
.extractDataFromGEM <- function(datasetObject, geneNames, keys.sep=",") {
  #Set a null object for later rbind call
  tempExprs <- NULL
  
  #Expand the expression object for keys which map to multiple genes
	tempExprs2 <- .expand.df.sampleLevel(datasetObject$expr, datasetObject$keys)
	setkey(tempExprs2, keys)
	tempExprs2 <- tempExprs2[unique(geneNames)]
	tempSplit <- split(tempExprs2, by="keys")
	
	tempUnsplit <- lapply(tempSplit, .summarizeSplit)
	tempExprs <- matrix(unlist(tempUnsplit), nrow=length(tempUnsplit), byrow=T)
	rownames(tempExprs) <- names(tempUnsplit)
	colnames(tempExprs) <- make.names(colnames(tempUnsplit[[1]]))
	return(data.frame(tempExprs))
}

.summarizeSplit <- function(curSplit) {
  if(length(curSplit) == 0) {
    return()
  }
  if(!is.vector(curSplit)) {
    return(t(as.matrix(colMeans(curSplit[,1:(ncol(curSplit)-1)]) )))
  } else {
    return(t(as.matrix(curSplit[1:(length(curSplit)-1)])))
  }
}

.expand.df.sampleLevel <- function (df, keys=keys, keys.sep = ",") {
  ## every row in df may have multiple identities
  ## these identities are defined in keys, with individual keys separated by keys.sep
  ## this function expands every row that has 1:n mapping to n x 1:1 rows
  
  skey           <- strsplit( keys, split = keys.sep)
  df             <- df[rep(1:nrow(df), sapply(skey, length)), ]
  expr=data.table(df)
  expr$keys <- unlist(skey)
  return(expr)
}


.replaceValues <- function(x, thresholdValue = 0, replaceValue = 1) {
	for(i in 1:dim(x)[2]) {
		indices = which(x[,i] <= thresholdValue)
		if(length(indices) == 0)
			next
		x[indices,i] = replaceValue
	}
	return(x)
}

.replaceNaNs <- function(x, replaceValue=1) {
	for(i in 1:dim(x)[2]) {
		indices = which(is.nan(x[,i]) == TRUE)
		if(length(indices) == 0)
			next
		x[indices,i] = replaceValue
	}
	return(x)  
}

.replaceNAs <- function(x, replaceValue=1) {
	for(i in 1:dim(x)[2]) {
		indices = which(is.na(x[,i]) == TRUE)
		if(length(indices) == 0)
			next
		x[indices,i] = replaceValue
	}
	return(x)  
}


#declare global variables for variables in data.table/with notation to avoid R CMD CHECK notes
utils::globalVariables(c("keys"))
