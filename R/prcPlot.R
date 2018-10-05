#' Plot the PRC Curve for a Dataset
#' @description
#' prcPlot will plot a Precision-Recall curve (and return the AUPRC) that describes how well a gene signature (as defined in a \code{filterObject}) classifies groups in a dataset (in the form of a \code{datasetObject}).
#' 
#' @param filterObject a metaFilter object containing the signature genes that will be used for calculating the score
#' @param datasetObject a Dataset object for group comparison in the PRC plot. (At least, must have a \code{$expr} of probe-level data, \code{$keys} of probe:gene mappings, and \code{$class} of two-class labels.)
#' @param title title of the plot (default: \code{datasetObject$formattedName})
#' @param subtitle subtitle of the figure
#' @param textSize use this to easily increase or decrease the size of all the text in the plot
#' @param rounding how many digits to round the AUPRC and CI to (default: 3)
#' @param curveColors \emph{Graphical:} the color for the PRC curves (default: "red")
#' @param legend \emph{Graphical:} if TRUE, a legend will be included
#' @param PRC.lty \emph{Graphical:} PRC curve line type
#' @param PRC.lwd \emph{Graphical:} PRC curve line width
#' @param backgroundColor \emph{Graphical:} background color of the plot
#' @param grid.marks \emph{Graphical:} increment between grid lines
#' @param grid.color \emph{Graphical:} grid line color
#' @param grid.lty \emph{Graphical:} grid line type
#' @param grid.lwd \emph{Graphical:} grid line width
#' @param legend.lty \emph{Graphical:} legend style (0 is no box, 1 is boxed legend)
#' @param cex.main \emph{Graphical:} title size
#' @param cex.subtitle \emph{Graphical:} subtitle size
#' @details Evaluates the ability of a given gene set to separate two classes. As opposed to ROC curves, PRC curves are more sensitive to class imbalances.
#' The gene set is evaluated as a Z-score of the difference in means between the positive genes and the negative genes (see calculateScore). 
#' @return Returns a standard PRC plot, plus AUPRC with 95\% CI (calculated with the trapezoid method).
#' @seealso	
#'         \code{\link{multiplePRCPlot}}, \code{\link{rocPlot}}
#' @examples 
#' prcPlot(tinyMetaObject$filterResults[[1]], tinyMetaObject$originalData[[1]])
#' @export
#' @import ROCR
#' @author Aditya Rao & Jiaying Toh
prcPlot <- function(filterObject, datasetObject,title = datasetObject$formattedName,
                    subtitle=NULL, textSize = NULL, rounding=3, curveColors="red", legend=TRUE, PRC.lty=1,
                    PRC.lwd=1, backgroundColor="gray93", grid.marks=0.1, grid.color="white",
                    grid.lty=1, grid.lwd=0.9, legend.lty=0, cex.main=1, cex.subtitle=0.9){
  scores = calculateScore(filterObject,datasetObject,TRUE)
  pred_PRC = prediction(scores, datasetObject$class)
  perf = performance(pred_PRC, "prec", "rec")
  #add extra points at (1,0) and (0,1) to complete the curve
  perf@x.values[[1]] = c(0,perf@x.values[[1]],1)
  perf@y.values[[1]] = c(1,perf@y.values[[1]],0)
  auprcInfo=.calcauprc(perf@x.values[[1]],perf@y.values[[1]],datasetObject$class)
  
  legendText = sprintf("%s AUPRC=%s (95%% CI %s-%s)",datasetObject$formattedName,round(auprcInfo$auprc,rounding),
                       round(auprcInfo$auprc.CI[1],rounding),round(auprcInfo$auprc.CI[2],rounding))
  
  if(is.null(textSize)){
    plot(perf, main = title, cex.main = cex.main, ylim=c(0,1),xlim=c(0,1))
    graphics::mtext(subtitle, line = 0.5, cex=cex.subtitle)
  }else{
    plot(perf, main = title, ylim = c(0, 1), xlim = c(0,1),
         cex.lab = textSize, cex.axis = textSize, cex.main = textSize, cex.sub = textSize)
    graphics::mtext(subtitle, line = 0.5, cex=textSize)
  }
  graphics::rect(graphics::par("usr")[1],graphics::par("usr")[3],graphics::par("usr")[2],graphics::par("usr")[4],col = backgroundColor)
  graphics::abline(h=seq(0, 1, grid.marks), v=seq(0, 1, grid.marks), col=grid.color, lty=grid.lty, lwd=grid.lwd)
  #lines(c(1,0),c(0,1),col = "black", lty = diag.lty)
  plot(perf, col = curveColors, add = T, lty=PRC.lty, lwd=PRC.lwd)
  
  if(legend){
    if(is.null(textSize)){
      legend("bottomleft", inset = 0.03, legend = legendText, fill = curveColors, box.lty = legend.lty)
    }else{
      legend("bottomleft", inset = 0.03, legend = legendText, fill = curveColors, box.lty = legend.lty, cex = textSize * 0.75)
    }
  }
}



###-###-###-###-###-###
###  .calcauprc()   ###
###-###-###-###-###-###

#DESCRIPTION
#function for calculating the auprc of a precision-recall curve, as well as the
#confidence interval

#PARAMETERS
#precision - vector of precision values
#recall - vector of recall values
#class - vector of classifications for your samples, with 0 representing your negative class (controls) and 1 representing your positive class (cases)
#type.CI - method used to calculate the confidence interval. binomial is more balanced, but logit is guaranteed to be between 0 and 1
#conflevel - confidence level (percentage) of the confidence interval, e.g. 0.95 will produce a 95% confidence interval. Must be between 0 and 1.

#REQUIRED PACKAGES: zoo, pracma

.calcauprc <- function(precision, recall, class, type.CI="binomial", conflevel=0.95){

  ## calcauprc math
  NA.vals=union(which(is.na(precision)),which(is.na(recall)))
  if(length(NA.vals)==0){
    auprc = sum(diff(precision)*zoo::rollmean(recall,2))
  }else{
    auprc = sum(diff(precision[-NA.vals])*zoo::rollmean(recall[-NA.vals],2))
  }
  
  if(conflevel > 0 && conflevel < 1){
    conflevel = 1-((1-conflevel)/2)
    phi = stats::qnorm(conflevel,0,1)
  }else{
    warning("conflevel must be between 0 and 1\nJust the AUPRC will be returned, with no confidence interval\n")
    return(auprc)
  }
  
  if(type.CI=="binomial"){
    n.pos = sum(class==1)
    ci.lower = auprc-phi*sqrt((auprc*(1-auprc))/n.pos)
    ci.upper = auprc+phi*sqrt((auprc*(1-auprc))/n.pos)
    if(ci.upper>1){ci.upper=1}
    if(ci.lower<0){ci.lower=0}
  }else if(type.CI=="logit"){
    mu = pracma::logit(auprc)
    n.pos = sum(class==1)
    t = 1/sqrt(n.pos*auprc*(1-auprc))
    e=exp(1)
    ci.lower = (e^(mu-phi*t))/(1+e^(mu-phi*t))
    ci.upper = (e^(mu+phi*t))/(1+e^(mu+phi*t))
  }else{
    return(auprc)
  }
  
  return(list(auprc=auprc,auprc.CI=c(ci.lower,ci.upper)))
}
