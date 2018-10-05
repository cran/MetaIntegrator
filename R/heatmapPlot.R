#' Generates a heatmap with effect sizes for all genes which pass a filter in all measured diseases
#' @param metaObject a Meta object which must have the $originalData, $metaAnalysis populated
#' @param filterObject a MetaFilter object containing the signature genes that will be used for the heatmap
#' @param colorRange a vector of length two with the minimum and maximum values for the heatmap colors. (default: c(-1,1))
#' @param geneOrder FALSE if the genes should be ordered by pooled effect size in this datasets. Otherwise, the ordered names of the genes. (default: FALSE)
#' @param datasetOrder FALSE if the datasets should be ordered alphabetically. Otherwise, the ordered names of the datasets (default: FALSE)
#' @param displayPooled TRUE if the pooled effect sizes should be displayed. (default: TRUE)
#' @param useFormattedNames TRUE if the formatted datasetNames should be displayed. (default: TRUE)
#' 
#' @import  ggplot2 
#' @return Generates a heatmap with effect sizes for all genes which pass a filter
#' 
#' @author Winston A. Haynes
#' 
#' @examples
#' heatmapPlot(tinyMetaObject, tinyMetaObject$filterResults[[1]])
#' @export
heatmapPlot <- function(metaObject, filterObject, colorRange=c(-1,1), geneOrder=FALSE, datasetOrder=FALSE, 
												displayPooled=TRUE, useFormattedNames=TRUE) {
  geneNames <- c(filterObject$posGeneNames, filterObject$negGeneNames)
  geneNames <- intersect(geneNames, rownames(metaObject$metaAnalysis$datasetEffectSizes))
  my_palette <- grDevices::colorRampPalette(c("purple", "white", "orange"))(n = 1000)
  
  if(length(geneNames)<=1) {
    warning("Must have at least 2 genes for heatmap") 
    return()
  }
  
  #If dataset order was provided, use that here
  datasetEffectSizes <- metaObject$metaAnalysis$datasetEffectSizes[geneNames, ]
  if(length(datasetOrder) > 1) {
    datasetEffectSizes <- datasetEffectSizes[,datasetOrder]
  }
  
  heatFrame <-as.data.frame(datasetEffectSizes)
  if(displayPooled) {
    #Combine the pooled effect sizes with the data
    heatFrame <- cbind(data.frame(Pooled=metaObject$metaAnalysis$pooledResults[geneNames, "effectSize"]), 
                       heatFrame)
  }
  
  #Use the formatted names in the plot
  nameFun<- function(obj) { obj$formattedName }
  
  if(length(datasetOrder)==1) {
    formattedNames <-sapply(metaObject$originalData, nameFun)
  } else {
    formattedNames <- sapply(metaObject$originalData[datasetOrder], nameFun)
  }
  if(displayPooled) {
    formattedNames <- c("Pooled", formattedNames)
  }
  if(useFormattedNames) {
  	colnames(heatFrame) <- formattedNames
  }
  if(length(datasetOrder) == 1) {
    datasetOrder <- colnames(heatFrame)[order(as.character(colnames(heatFrame)))]
    if(displayPooled) {
      datasetOrder <- c("Pooled", setdiff(datasetOrder, c("Pooled")))
    }
  } else {
    datasetOrder <- colnames(heatFrame)
  }
  
  
  if(length(geneOrder) ==1 ) {
  	if(displayPooled) {
  		#Order the rows by pooled effect sizes 
      geneOrder <- rownames(heatFrame[order(heatFrame$Pooled),])
  	} else {
    	geneOrder <- rownames(heatFrame)
    }
  }
  
  heatFrame$geneName <- row.names(heatFrame)
  heatMelt <- reshape2::melt(heatFrame, id.vars=c("geneName"))
  
  #Threshold the maximum values
  heatMelt$value[which(heatMelt$value < colorRange[[1]])] <- colorRange[[1]]
  heatMelt$value[which(heatMelt$value > colorRange[[2]])] <- colorRange[[2]]
  
  #Set up factors so that rows order appropriately
  heatMelt$geneName <- factor(heatMelt$geneName, levels=geneOrder)
  heatMelt$variable <- factor(heatMelt$variable, levels=rev(datasetOrder))
  
  heatmapObj <- ggplot(heatMelt, aes_string(x='geneName', y='variable', fill='value')) + geom_tile() + scale_fill_gradient2(low = "purple", mid="white", high="orange", na.value = "white") + 
    xlab("") + ylab("") + theme(axis.line=element_blank(), axis.ticks = element_blank())
  if(nrow(heatFrame)>40) {
    heatmapObj <- heatmapObj + theme(axis.text.x = element_blank())
  } else {
    heatmapObj <- heatmapObj + theme(axis.text.x = element_text(angle=20, hjust = 1, vjust=1))
  }
  return(heatmapObj)
}
