#' A wrapper function to run COCONUT on the MetaIntegrator objects.
#' 
#' @param metaObject a MetaIntegrator formatted Meta object. 
#' @param ... pass along arguments to COCONUT
#' 
#' @return Results from COCONUT analysis on the MetaIntegrator object
#' 
#' @export
#' @author Winston A. Haynes
#'

coconutMetaIntegrator <- function(metaObject, ...) {
	return(COCONUT::COCONUT(GSEs=geneLevelMeta(metaObject), control.0.col="class", ...))
}

geneLevelMeta <- function(metaObject) {
  geneLevelMeta <- list()
  for(datasetName in names(metaObject$originalData)) {
    datasetObject <- metaObject$originalData[[datasetName]]
    datasetPheno <- datasetObject$pheno
    rownames(datasetPheno) <- make.names(rownames(datasetPheno))
    datasetPheno$class <- datasetObject$class
    geneLevelData <- getSampleLevelGeneData(datasetObject, unique(datasetObject$keys))
    geneLevelDataset <- list(pheno=datasetPheno, genes=geneLevelData)
    geneLevelMeta[[length(geneLevelMeta)+1]] <- geneLevelDataset
    names(geneLevelMeta)[length(geneLevelMeta)] <- datasetName
  }
  return(geneLevelMeta)
}
