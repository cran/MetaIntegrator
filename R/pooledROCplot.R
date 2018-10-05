#' Generate a plot with a pooled ROC curve
#' @description
#' Given a \code{metaObject} with \code{$originalData} populated, this function calculates and
#' plots a "pooled" ROC curve that represents the average of all the individual ROC
#' curves. This version of the function is for use with MetaIntegrator.
#'
#' @param metaObject a metaObject which must have \code{metaObject$originalData} populated with a list of \code{datasetObjects} that will be used for discovery
#' @param filterObject a metaFilter object containing the signature genes that will be used for calculating the score
#' @param points number of points to simulate for the approximated ROC curves during the linear interpolation (default: 1000)
#' @param weighting when calculating the mean AUC, if \code{weighting}=TRUE then the weighted mean AUC is calculated (default: TRUE)
#' @param title title of the plot
#' @param size size of the text/legend/etc (default: 14)
#' @param rounding how many digits to round the AUC and CI to (default: 3)
#' @param smoothed if TRUE, then a smoothed ROC curve is estimated using a modified version of the Kester and Buntinx Method
#' @param auc1.thresh (if \code{smoothed}=TRUE) if the AUC of a dataset is above this threshold, then it is treated as if the AUC were 1 (default: 0.99)
#' @param bootReps (if \code{smoothed}=TRUE) number of bootstrap iterations (default: 1000)
#' @param minPoints (if \code{smoothed}=TRUE) minimum number of points required for bootstrap to be used (default: 5)
#' @param numCores (if \code{smoothed}=TRUE) number of CPUs to use if parallel computing is desired (default: 1)
#' @param method (if \code{smoothed}=TRUE) method used to compute summary meta-statistics (default: "random")
#' @details
#' To make sure the input is correctly formatted, the input \code{metaObject} should be checked with
#' \code{checkDataObject(metaObject, "Meta", "Pre-Analysis")} before starting the meta-analysis.
#'
#' By default, this average ROC curve is calculated by first using linear
#' interpolation to create approximated versions of each given ROC curve that all
#' have the same set of FPR values. A pooled ROC curve is then calculated by
#' taking the weighted mean of the corresponding TPR values (weighting corresponds
#' to the number of samples in each dataset). This pooled curve is represented as
#' a black curve. In addition, the weighted standard deviation is calculated for
#' each TPR, which is represented by a grey area on the plot. The pooled AUC is
#' calculated by using the trapezoid method on the pooled ROC curve, and the 95\%
#' confidence interval of the pooled AUC is calculated using the pooled standard
#' error of the individual ROC curves.
#'
#' If \code{smoothed}=TRUE, then a smoothed version of the pooled ROC curve will be
#' plotted instead, with the surrounding gray area representing the weighted
#' standard deviation of the pooled ROC curve. The statistics for this smoothed
#' curve are based on the Kester and Buntinx Method, from (Kester and Buntinx,
#' \emph{Med Decis Making}, 2000). Methods have been added by Tim Sweeney (2015) for better
#' estimates in cases with low numbers of tpr/fpr values. Methods have also been
#' added by Aditya Rao (2018) to predict the curve's alpha parameter for a given
#' beta parameter and AUC, as well as to calculate the weighted standard deviation
#' of the given ROC curves.
#' @return Generates a plot with each individual ROC curve as well as the pooled ROC curve
#' @seealso	\code{\link{summaryROCPlot}}
#' @examples
#' pooledROCPlot(tinyMetaObject, filterObject = 
#'    tinyMetaObject$filterResults$pValueFDR0.05_es0_nStudies1_looaTRUE_hetero0)
#' @export
#' @import data.table ggplot2 ROCR
#' @references Kester and Buntinx, \emph{Med Decis Making}, 2000
#' @author Aditya Rao (with help from Hayley Warsinske and Francesco Vallania, original idea from Madeleine Scott, and some code adapted from Tim Sweeney)
pooledROCPlot <- function(metaObject,filterObject,points=1000,weighting=TRUE,title=NULL,size=14,rounding=3,
                          smoothed=FALSE,auc1.thresh=0.99,bootReps=1000,minPoints=5,numCores=1,method="random"){

  #calculate scores/ROC once at the beginning to save time
  metaObject$originalData = lapply(metaObject$originalData, function(dataset){
    dataset$scores = calculateScore(filterObject,dataset,TRUE)
    dataset$ROCinfo = calculateROC(as.numeric(as.character(dataset$class)), as.numeric(dataset$scores))
    return(dataset)
  })

  pooledStats = do.call(rbind,lapply(metaObject$originalData, function(dataset){
    auc = dataset$ROCinfo$auc
    auc.lo = max(dataset$ROCinfo$auc.CI[1],0)
    auc.up = min(dataset$ROCinfo$auc.CI[2],1)
    auc.se = .aucSE(auc,dataset$class)
    return(data.frame(name=dataset$formattedName,AUC=auc,CI.lower=auc.lo,CI.upper=auc.up,AUC.SE=auc.se,N=length(dataset$class)))
  }))

  legendText = sprintf("%s AUC=%s (95%% CI %s-%s)",pooledStats$name,round(pooledStats$AUC,rounding),
                       round(pooledStats$CI.lower,rounding),round(pooledStats$CI.upper,rounding))

  curves = data.table(do.call(rbind,lapply(metaObject$originalData, function(dataset){
    interp = stats::approx(dataset$ROCinfo$roc$x,dataset$ROCinfo$roc$y,n=points)
    new.curve = data.frame(name = rep(dataset$formattedName,points),N = length(dataset$class), FPR = interp$x, TPR = interp$y)
    new.curve = new.curve[points:1,]
    return(new.curve)
  })))

  if(weighting){
    roc.summ = curves[,.(TPR=sum(TPR*(pooledStats$N/sum(pooledStats$N)))),by=FPR]
  }else{
    roc.summ = curves[,.(TPR=mean(TPR)),by=FPR]
  }

  #weighted standard deviation w/ unbiased estimator (see https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Weighted_sample_variance)
  #this is with reliability weights (frequency would mean treating every TPR value as if it were repeated N_i times)
  weights = pooledStats$N/sum(pooledStats$N)
  v1 = sum(weights) #should sum to 1
  v2 = sum(weights^2)
  denom = v1-(v2/v1)
  roc.SD = curves[,.(TPR.SD=sqrt(sum(weights*(TPR-(sum(TPR*(weights/v1))))^2)/denom)),by=FPR]$TPR.SD

  roc.upper = data.frame(FPR = roc.summ$FPR,TPR = roc.summ$TPR + roc.SD)
  roc.lower = data.frame(FPR = roc.summ$FPR,TPR = roc.summ$TPR - roc.SD)
  roc.upper[roc.upper>1]=1
  roc.lower[roc.lower<0]=0

  # Alterative measures of spread/uncertainty/etc.
  #
  # standard error
  # roc.SE = curves[,.(TPR.SE=sd(TPR)/sqrt(mean(pooledStats$N))),by=FPR]$TPR.SE
  #
  # calculating standard error using number of datasets
  # roc.SE = curves[,.(TPR.SE=sd(TPR)/sqrt(length(metaObject$originalData))),by=FPR]$TPR.SE
  #
  # calculating standard error using total sample size
  # roc.SE = curves[,.(TPR.SE=sd(TPR)/sqrt(sum(pooledStats$N))),by=FPR]$TPR.SE
  #
  # 95% confidence interval
  # roc.upper = data.frame(FPR = roc.summ$FPR,TPR = roc.summ$TPR + (roc.SE*1.96))
  # roc.lower = data.frame(FPR = roc.summ$FPR,TPR = roc.summ$TPR - (roc.SE*1.96))
  #
  # standard deviation
  # roc.SD = curves[,.(TPR.SE=sd(TPR)),by=FPR]$TPR.SE
  # roc.upper = data.frame(FPR = roc.summ$FPR,TPR = roc.summ$TPR + roc.SD)
  # roc.lower = data.frame(FPR = roc.summ$FPR,TPR = roc.summ$TPR - roc.SD)

  auc.summ = .aucROCframe(data.frame(roc.summ))
  #calculate pooled standard error (https://stats.stackexchange.com/questions/55999/is-it-possible-to-find-the-combined-standard-deviation)
  auc.se = sqrt(sum((pooledStats$N-1)*pooledStats$AUC.SE^2+pooledStats$N*(pooledStats$AUC-as.vector(auc.summ))^2)/(sum(pooledStats$N)-1))
  auc.lower = max(0,auc.summ - (auc.se * 1.96))
  auc.upper = min(1,auc.summ + (auc.se * 1.96))

  if(smoothed){
    #use this if you want to calculate weighted standard deviation of AUCs
    AUC.wSD = sqrt(sum(weights*(pooledStats$AUC-(sum(pooledStats$AUC*(weights/v1))))^2)/denom)
    auc.lower2 = max(0,auc.summ - AUC.wSD)
    auc.upper2 = min(1,auc.summ + AUC.wSD)

    # use this if you want to use the AUCs of the current roc.lower and roc.upper instead
    # auc.lower2 = .aucROCframe(roc.lower)
    # auc.upper2 = .aucROCframe(roc.upper)

    #get beta parameter for the smoothed ROC curve
    ROC.stats = .getMetaROCStats(metaObject,filterObject,numCores=numCores,minPoints=minPoints,bootReps=bootReps,auc1.thresh=auc1.thresh)
    ROC.stats = stats::na.omit(ROC.stats)
    if(nrow(ROC.stats)>1){
      beta <- with(ROC.stats, rmeta::meta.summaries(d=tstar_beta, se=SE_beta*sqrt(N), method=method))
    } else {
      beta <- with(ROC.stats, data.frame(summary=tstar_beta, se.summary=SE_beta))
    }

    alpha.summ = .getSummROCalpha(auc.summ,beta$summary)
    alpha.lower = .getSummROCalpha(auc.lower2,beta$summary)
    alpha.upper = .getSummROCalpha(auc.upper2,beta$summary)

    roc.summ <- .newMetaROC(alpha=alpha.summ, beta=beta$summary)
    roc.lower <- .newMetaROC(alpha=alpha.lower, beta=beta$summary)
    roc.upper <- .newMetaROC(alpha=alpha.upper, beta=beta$summary)
  }

  #add points at (0,0) and (1,1)
  roc.summ = rbind(list(1,1),roc.summ,list(0,0))
  roc.lower = rbind(list(1,1),roc.lower,list(0,0))
  roc.upper = rbind(list(1,1),roc.upper,list(0,0))

  plotData = do.call(rbind,lapply(metaObject$originalData, function(dataset){
    return(data.frame(name=dataset$formattedName,FPR=dataset$ROCinfo$roc$x,TPR=dataset$ROCinfo$roc$y))
  }))
  plotData = rbind(plotData,data.frame(name="Pooled",FPR=roc.summ$FPR,TPR=roc.summ$TPR))

  legendText = c(legendText,sprintf("Pooled AUC=%s (95%% CI %s-%s)",round(auc.summ,rounding),
                                    round(auc.lower,rounding),round(auc.upper,rounding)))
  #force color to dark grey on palette
  hues <- seq(15, 375, length=length(metaObject$originalData)+1)
  dataPal = c(grDevices::hcl(h=hues, l=65, c=100)[1:length(metaObject$originalData)],"grey25")

  ggplot(plotData, aes_string(x = 'FPR', y = 'TPR')) +
    #Have this first to put it on bottom
    suppressWarnings(geom_ribbon(data=cbind(lower=roc.lower, upper=roc.upper),
                                 aes_string(x= 'lower.FPR', y='lower.TPR', ymin='lower.TPR', ymax='upper.TPR'), fill="gray75")) +
    geom_line(aes_string(colour = 'name')) +
    geom_segment(x=0, y=0, xend=1, yend=1, linetype=5, color="grey20") +
    scale_x_continuous("False Positive Rate (1-Specificity)",breaks=seq(0,1,0.2)) +
    scale_y_continuous("True Positive Rate (Sensitivity)",breaks=seq(0,1,0.2)) +
    scale_color_manual(values=dataPal, labels = legendText) +
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
    guides(colour=guide_legend(override.aes = list(size=2.8))) +
    geom_line(data=roc.summ[nrow(roc.summ):1,], aes_string(x='FPR', y='TPR'), size=1.3, color="gray20")
  ## If just want to define upper and lower lines without fill:
  #geom_line(data=roc.upper, aes(x=FPR, y=TPR), size=1.1, linetype=3, color="gray15") +
  #geom_line(data=roc.lower, aes(x=FPR, y=TPR), size=1.1, linetype=3, color="gray15")
}



#######################################################################################################################################

###-###-###-###-###-###-###-##
###   .getMetaROCStats()   ###
###-###-###-###-###-###-###-##

#DESCRIPTION
#This function calculates meta-statistics for multiple ROC curves, which can be
#used to make a summary ROC curve for those datasets

#Statistics are based on the Kester and Buntinx Method, from (Kester and
#Buntinx, Med Decis Making, 2000). Methods have also been added by Tim Sweeney
#(2015) for better estimates in cases with low numbers of tpr/fpr values

#PARAMETERS
#metaObject - a metaObject which must have metaObject$originalData populated with a list of datasetObjects that will be used for discovery
#filterObject - a MetaFilter object containing the signature genes that will be used for calculating the score
#auc1.thresh - if the AUC of a dataset is above this threshold, then it is treated as if the AUC were 1
#bootReps - number of bootstrap iterations
#minPoints - minimum number of points required for bootstrap to be used
#numCores - number of CPUs to use if parallel computing is desired

#RETURN VALUE
#data frame listing the sample size, alpha parameter value and standard error,
#beta parameter and standard error, auc value, 95% confidence interval upper and
#lower values, and the auc standard error for each of the ROC curves

#REQUIRED PACKAGES: MetaIntegrator, ROCR, boot

.getMetaROCStats <- function(metaObject, filterObject, auc1.thresh=0.99,bootReps=1000, minPoints=5, numCores=1){
  nameVec=rep("",length(metaObject$originalData))
  for(i in 1:length(metaObject$originalData)){
    nameVec[i]=metaObject$originalData[[i]]$formattedName
  }

  ROC.stats = do.call(rbind,mclapply(mc.cores=numCores, metaObject$originalData, function(dataset){
    auc = dataset$ROCinfo$auc
    auc.lo = max(dataset$ROCinfo$auc.CI[1],0)
    auc.up = min(dataset$ROCinfo$auc.CI[2],1)
    auc.se = .aucSE(auc,dataset$class)
    pred_ROC = prediction(dataset$scores, dataset$class)
    perf = performance(pred_ROC, "tpr", "fpr")
    P=sum(dataset$class==1)
    N=sum(dataset$class==0)
    if(auc>auc1.thresh){
      #as defined in the paper. alpha = log((P+0.5)*(N+0.5)/0.25)
      alpha = log((P+0.5)*(N+0.5)/0.25)
      #ASE = sqrt of asymptotic variance over sqrt n
      ASE = (1/0.5 + 1/0.5 + 1/(P+0.5) + 1/(N+0.5))^(-1/2)/sqrt(P+N)
      return(c(N+P, alpha, ASE , 0, 10, auc, auc.lo, auc.up, auc.se))
    }else if(auc==0){ #not sure why no threshold here
      #added for edge case of auc=0
      alpha = -log((P+0.5)*(N+0.5)/0.25)
      ASE = (1/0.5 + 1/0.5 + 1/(P+0.5) + 1/(N+0.5))^(-1/2)/sqrt(P+N)
      return(c(N+P, alpha, ASE , 0, 10, auc, auc.lo, auc.up, auc.se))
    }

    sens = unlist(perf@x.values); min1.spec = unlist(perf@y.values)
    sens[min(which(sens==1))] = 0.98 #not convinced this is necessary
    min1.spec[max(which(min1.spec==0))] = 0.02 #not convinced this is necessary
    #weights defined in paper appendix, sends inf/NAN to 0
    weights = (1/unlist(pred_ROC@tp) + 1/unlist(pred_ROC@tn) + 1/unlist(pred_ROC@fp) + 1/unlist(pred_ROC@fn))^(-1)
    weights[is.nan(weights) | is.na(weights)] = 0
    S = log(sens/(1-sens)) + log(min1.spec/(1-min1.spec))
    D = -1*(log(sens/(1-sens)) - log(min1.spec/(1-min1.spec)))

    #if there are enough non-0 weights, then bootstrap
    points = sum(weights != 0)
    #tim made some note about needing points that vary in x and y, not just non-0
    #weights, but then he commented out the corresponding code so it looks like
    #maybe non-0 weights is enough? If there's an issue, can try this code instead:
    #points = sum(!duplicated(sens) & !duplicated(min1.spec))

    if(points > minPoints){
      bs <- function(formula, data, indices){
        d = data[indices, ]
        fit = stats::lm(formula, data=d, weights=weights)
        return(stats::coef(fit))
      }
      boot.out = boot::boot(data=data.frame(min1.spec,sens,weights,S,D), statistic=bs, formula=D~S, R=bootReps, weights=weights)
      op <- NULL
      for (i in 1:2) op <- rbind(op, boot::imp.moments(boot.out, index = i)$rat)
      std.error <- sqrt(op[, 2]) #this used to be 2L
      return(c(N+P, op[1,1], std.error[1], op[2,1], std.error[2], auc, auc.lo, auc.up, auc.se))
    } else if (points >= 2){
      #else assign ASE as standard errors for both alpha and beta
      cat("\tCan't bootstrap, points < minSize...\t")
      model <- stats::lm(D~S, data=data.frame(min1.spec,sens,weights,S,D), weights=weights)
      alpha <- model$coefficients[1]
      beta <- 0
      ASE <- (1/0.5 + 1/0.5 + 1/(P+0.5) + 1/(N+0.5))^(-1/2)/sqrt(P+N)
      return(c(N+P, alpha, ASE , beta, ASE, auc, auc.lo, auc.up, auc.se))
    } else {
      #if only one point in ROC that is not on axis, can't use.
      cat(sprintf(" has <2 points in ROC curve; cannot compute summary stats\n", dataset$formattedName))
      return(c(N+P, rep(NA, 4), auc, auc.lo, auc.up))
    }
  }))

  colnames(ROC.stats)=c("N","tstar_alpha","SE_alpha","tstar_beta","SE_beta","AUC","CI.lower","CI.upper","AUC.SE")
  rownames(ROC.stats)=nameVec
  return(data.frame(ROC.stats))
}



###-###-###-###-##-#
###   .aucSE()   ###
###-###-###-###-##-#

#DESCRIPTION
#calculates the standard error given an AUC and vector of classification labels

#PARAMETERS
#auc - AUC value
#class - vector of classifications for your samples, with 0 representing your negative class (controls) and 1 representing your positive class (cases)

#RETURN VALUE
#Returns the standard error of the AUC

#REQUIRED PACKAGES: None

.aucSE <- function(auc,class){
  n.pos = sum(class == 1)
  n.neg = sum(class == 0)
  q1 = auc/(2 - auc)
  q2 = (2 * auc^2)/(1 + auc)
  se.auc = sqrt((((auc * (1 - auc))+(n.pos - 1)*(q1 - auc^2))+((n.neg - 1)*(q2 - auc^2)))/(n.pos * n.neg))
  return(se.auc)
}



###-###-###-###-###-###-###-##
###   .getSummROCalpha()   ###
###-###-###-###-###-###-###-##

#DESCRIPTION
#for the Kester and Buntinx method, calculates the alpha parameter given the AUC and beta parameter

#PARAMETERS
#auc - AUC value
#beta - beta parameter

#RETURN VALUE
#returns the predicted alpha parameter

#REQUIRED PACKAGES: None

#for the Kester and Buntinx method, get the alpha parameter given the AUC and beta parameter
.getSummROCalpha <- function(auc,beta){
  if(auc==0.5){
    alpha = 0
  }
  else if(beta>=0.95){
    alpha = (0.01-log((1-auc)/auc))/(0.51)
  }
  else{
    alpha = (0.01-log((1-auc)/auc))/(0.681-(beta^2)*0.197)
  }

  if(alpha == Inf){
    return(12) #for now, to avoid errors, I return a very high alpha if alpha was Inf
  }else if(alpha == -Inf || alpha < -12)
    return(-12)
  else{
    return(alpha)
  }
}


#declare global variables for variables in data.table/with notation to avoid R CMD CHECK notes
utils::globalVariables(c("TPR","FPR","."))
