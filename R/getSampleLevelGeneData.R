
getSampleLevelGeneData <- function(datasetObject, geneNames){
	#Subfunctions all contained within for cleanliness
	GenesMtx <- .extractDataFromGEM(datasetObject, geneNames)
	GenesMtx <- .replaceValues(GenesMtx, 0, 1)
	GenesMtx <- .replaceNaNs(GenesMtx, 1)
	GenesMtx <- .replaceNAs(GenesMtx, 1)
	return(GenesMtx)
}
.extractDataFromGEM <- function(datasetObject, geneNames) {
# 	genes <- filterObject
	tempExprs = NULL
	junk = lapply(as.matrix(geneNames), function(x, keys) which(keys == x), keys=datasetObject$keys)
	for(j in 1:length(junk)) {
		if(length(junk[[j]]) == 0) {
			next
		}
		temp = datasetObject$expr[junk[[j]],]
		if(!is.vector(temp)) {
			temp = t(as.matrix(colMeans(temp) ))
		} else {
			temp = t(as.matrix(temp))
		}
		rownames(temp) = geneNames[j]
		tempExprs = rbind(tempExprs, temp)
	}
	tempExprs = data.frame(tempExprs)
	
	return(tempExprs)
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
