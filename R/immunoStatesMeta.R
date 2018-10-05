#' immunoStates deconvolution analysis on MetaIntegrator object(s)
#' 
#' @param metaObject a MetaIntegrator formatted Meta object. 
#' 
#' @return Results from immunoStates stored in $originalData
#' 
#' @export
#' @import data.table 
#' @examples 
#' \dontrun{
#' # Example won't work on tinyMetaObject because it requires real gene names
#' # Download the needed datasets for processing. 
#' sleData <- getGEOData(c("GSE11909","GSE50635", "GSE39088"))
#' 
#' # Run immunoStates
#' immunoStatesEstimates <- immunoStateMeta(sleData)
#' }
#' @description 
#' Run immunoStates and load the results into $originalData for 
#' running meta-analysis on the cell proportion estimates.
immunoStatesMeta <- function(metaObject){
	
  #perform deconvolution
	immunoRes <- immunoStatesDecov(metaObject)
	
	#select for datasets that were deconvolved accurately
	iSnames <- names(which(sapply(immunoRes$immunoStates,function(i) length(i)>1)))
	
	#add check here
	if(length(setdiff(names(immunoRes$immunoStates),iSnames)) > 0){
	  #spit out a warning
	  excldNames <- setdiff(names(immunoRes$immunoStates),iSnames)
	  warning("These datasets were not deconvolved and will therefore be excluded:\n",
	          paste(excldNames,collapse="\n"))	  
	}else{
	  if(length(setdiff(names(immunoRes$originalData),names(immunoRes$immunoStates)))>0){
	    excldNames <- setdiff(names(immunoRes$originalData),names(immunoRes$immunoStates))
	    warning("These datasets were not deconvolved and will therefore be excluded:\n",
	            paste(excldNames,collapse="\n"))	 
	  }
	}
	
	#remove datasets that are not included
	metaObject$originalData <- metaObject$originalData[iSnames]
	
	#prepare for meta-analysis 
	for(datasetName in iSnames) {
		immunoExpr <- immunoRes$immunoStates[[datasetName]] 
		immunoExpr <- as.data.frame(immunoExpr)
		rownames(immunoExpr) <- immunoExpr$rn
		immunoExpr <- immunoExpr[,!(colnames(immunoExpr) %in% 
															c("rn", "P-value", "Correlation", "RMSE"))]
		immunoExpr <- t(immunoExpr)
		metaObject$originalData[[datasetName]]$expr        <- as.matrix(immunoExpr)
		metaObject$originalData[[datasetName]]$keys        <- rownames(immunoExpr)
		names(metaObject$originalData[[datasetName]]$keys) <- rownames(immunoExpr)
	}
	
	#remove everything that is not originalData
	for(objName in names(metaObject)){
		if(!(objName =="originalData")) {
			metaObject[[objName]] <- NULL
		}
	}
	return(metaObject)
}
