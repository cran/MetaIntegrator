#' Generate a plot with multiple PRC curves
#' @description
#' for each dataset in the metaObject, prcPlot will return a ggplot of a Precision-Recall curve (and return the AUPRC) that describes how well a gene signature
#' (as defined in a \code{filterObject}) classifies groups in a dataset (in the form of a \code{datasetObject}).
#' 
#' @param metaObject a metaObject which must have \code{metaObject$originalData} populated with a list of \code{datasetObjects} that will be used for discovery
#' @param filterObject a metaFilter object containing the signature genes that will be used for calculating the score
#' @param title title of the plot 
#' @param legend.names the name listed for each dataset in the legend (default: the \code{datasetObject$formattedName} for each dataset)
#' @param curveColors \emph{Graphical:} vector of colors for the PRC curves
#' @param size use this to easily increase or decrease the size of all the text in the plot
#' @details Each PRC plot evaluates the ability of a given gene set to separate two classes. As opposed to ROC curves, PRC curves are more sensitive to class imbalances.
#' The gene set is evaluated as a Z-score of the difference in means between the positive genes and the negative genes (see calculateScore). 
#' @return Returns a ggplot PRC plot for all datasets
#' @seealso	
#'         \code{\link{prcPlot}}, \code{\link{multipleROCPlot}}
#' @examples 
#' multiplePRCPlot(tinyMetaObject, filterObject = 
#'    tinyMetaObject$filterResults$pValueFDR0.05_es0_nStudies1_looaTRUE_hetero0)
#' @export
#' @import ggplot2
#' @author Aditya Rao, Andrew Bo Liu
multiplePRCPlot <- function(metaObject, filterObject, title = NULL, legend.names=NULL, curveColors=NULL, size=22){
  
  plotData = plyr::llply(metaObject$originalData, function(dataset) {.prcData(dataset, filterObject)})
  plotData = list(perf = plyr::ldply(plotData, function(res) data.frame(x=res$perf@x.values[[1]], y=res$perf@y.values[[1]])),
                  auprcInfo = plyr::ldply(plotData, function(res) data.frame(auprc=res$auprcInfo$auprc, auprc.CI.lo=res$auprcInfo$auprc.CI[1], auprc.CI.hi=res$auprcInfo$auprc.CI[2])))
  
  legendNames = if (is.null(legend.names) || length(legend.names) != length(metaObject$originalData)) names(metaObject$originalData) else legend.names
  annotation = with(plotData$auprcInfo, sprintf("AUPRC=%s (95%% CI %s-%s)", signif(auprc, 2), signif(auprc.CI.lo, 2), signif(auprc.CI.hi, 2)))
  
  if (is.null(curveColors) || length(curveColors) != length(metaObject$originalData)) {
    if (!is.null(curveColors)) {
      warning("Using brewer color scale because curveColors has length different from number of datasets")
    }
    colorScale = scale_colour_brewer(palette="Set1", breaks = legendNames, 
                                     labels = paste(legendNames, ": ", annotation, sep = ""))
  } else {
    colorScale = scale_colour_manual(breaks = legendNames, labels = paste(legendNames, ": ", annotation, sep = ""),
                                     values = curveColors)
  }

  p = ggplot(plotData$perf, aes_string(x = "x", y = "y")) +
    geom_line(aes_string(colour = ".id")) +
    #geom_segment(x=0, y=0, xend=1, yend=1, linetype=2) +
    scale_x_continuous("Recall") +
    scale_y_continuous("Precision") +
    colorScale +
    ggtitle(title) +
    theme(text = element_text(size=size)) +
    theme(plot.title = element_text(size=size), 
          axis.title.x = element_text(size=size),
          axis.title.y = element_text(size=size, angle=90),
          legend.justification=c(1,0), 
          legend.position=c(1,0),
          legend.title = element_blank(),
          legend.key = element_blank(),
          legend.text=element_text(size=size),
          axis.text.x = element_text(size=size),
          axis.text.y = element_text(size=size)) +
    guides(colour=guide_legend(override.aes = list(size=3)))
  
  return(p)
}

.prcData = function(dataset, filterObject) {
  scores = calculateScore(filterObject,dataset,TRUE)
  pred_PRC = ROCR::prediction(scores, dataset$class)
  perf = ROCR::performance(pred_PRC, "prec", "rec")
  #check for NaN bug
  if(is.nan(perf@y.values[[1]][1])){
    perf@y.values[[1]][1] = 1
  }
  #add extra points at (1,0) and (0,1) so that it looks nicer
  perf@x.values[[1]] = c(0,perf@x.values[[1]],1)
  perf@y.values[[1]] = c(1,perf@y.values[[1]],0)
  auprcInfo = .calcauprc(perf@x.values[[1]],perf@y.values[[1]],dataset$class)
  return(list(auprcInfo=auprcInfo,perf=perf))
}
