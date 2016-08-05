#' Generates a heatmap with effect sizes for all genes which pass a filter in all measured diseases
#' @param metaObject a Meta object which must have the $originalData, $metaAnalysis populated
#' @param filterObject a MetaFilter object containing the signature genes that will be used for the heatmap
#' @param colorRange a vector of length two with the minimum and maximum values for the heatmap colors. (default: c(-1,1))
#' 
#' @importFrom  gplots heatmap.2
#' @return Generates a heatmap with effect sizes for all genes which pass a filter
#' 
#' @author Winston A. Haynes
#' 
#' @examples
#' heatmapPlot(tinyMetaObject, tinyMetaObject$filterResults[[1]])
heatmapPlot <- function(metaObject, filterObject, colorRange=c(-1,1)) {
  geneNames <- c(filterObject$posGeneNames, filterObject$negGeneNames)
  geneNames <- intersect(geneNames, rownames(metaObject$metaAnalysis$datasetEffectSizes))
  my_palette <- colorRampPalette(c("purple", "white", "orange"))(n = 1000)
  
  if(length(geneNames)<=1) {
    warning("Must have at least 2 genes for heatmap") 
    return()
  }
  
  #Combine the pooled effect sizes
  heatmatrix <- cbind(data.frame(Pooled=metaObject$metaAnalysis$pooledResults[geneNames, "effectSize"]), 
                      metaObject$metaAnalysis$datasetEffectSizes[geneNames, ])
  
  
  #Use the formatted names in the plot
  nameFun<- function(obj) { 
    obj$formattedName }
  formattedNames <-c("Pooled", sapply(metaObject$originalData, nameFun))
  colnames(heatmatrix) <- formattedNames
  
  #Order the rows by pooled effect sizes 
  heatmatrix <- as.matrix(heatmatrix[order(heatmatrix$Pooled),c(order(as.character(colnames(heatmatrix))))])
  heatmatrix <- heatmatrix[,c(which(colnames(heatmatrix)=="Pooled"), which(colnames(heatmatrix)!="Pooled"))]
  heatmatrix[which(heatmatrix < colorRange[[1]])] <- colorRange[[1]]
  heatmatrix[which(heatmatrix > colorRange[[2]])] <- colorRange[[2]]
  heatmatrix <- t(heatmatrix)
  
  #Calculate margins based on the lengths of the labels we are displaying
  labCol <- if(ncol(heatmatrix)<40) colnames(heatmatrix) else FALSE
  marCol <- if(length(labCol)>1) max(sapply(colnames(heatmatrix), nchar))/1.5 else 2
  marRow <-  max(sapply(rownames(heatmatrix), nchar))/2
  
  heatmap.2(heatmatrix, na.rm=TRUE, col=my_palette,symbreaks=TRUE, density.info="none", trace="none",     
            dendrogram="none",
            key.title="", key.xlab="Effect Size", keysize=1,
            Colv=FALSE, Rowv=FALSE,
            na.color="white", sepwidth=c(0,0), labCol=labCol, labRow=rownames(heatmatrix), margins=c(marCol,marRow), cexRow=1, cexCol=1)
}
