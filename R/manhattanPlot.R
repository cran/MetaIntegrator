#' Generates a Manhattan plot with effect size FDR as y-axis
#' @param metaObject a Meta object which must have meta-analysis run
#' 
#' @return Generates a Manhattan plot with effect size FDR as y-axis
#' @author Winston A. Haynes
#' @export

manhattanPlot <- function(metaObject) {
  if(!checkDataObject(metaObject, objectType="Meta", objectStage = "Pre-Filter")) {
    stop("Need properly formatted Meta object")
  }
  pooledRes <- metaObject$metaAnalysis$pooledResults
  geneInfo <- snplist::getBioMartData(rownames(pooledRes), biomart="ENSEMBL_MART_ENSEMBL",
                                      host="grch37.ensembl.org",
                                      path="/biomart/martservice",
                                      dataset="hsapiens_gene_ensembl")
  pooledRes$gene <- rownames(pooledRes)
  mergedLocation <- merge(geneInfo,pooledRes, by="gene")
  
  chrFac <- as.factor(mergedLocation$chr)
  chrFac <- dplyr::recode(chrFac,
                          "MT" = "25",
                          "X" = "23",
                          "Y" = "24")
  chrFac <- as.character(chrFac)
  mergedLocation$chrFac <- as.numeric(chrFac)
  return(manhattanly::manhattanly(manhattanly::manhattanr(mergedLocation, chr="chrFac", bp="start", p="effectSizeFDR", gene="gene")))
}