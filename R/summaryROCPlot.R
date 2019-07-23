#' Generate a plot with a summary ROC curve
#' @param metaObject a Meta object which must have the $originalData populated
#' @param filterObject a MetaFilter object containing the signature genes that will be used for calculating the score
#' @param bootstrapReps number of bootstrap simulations to run for confidence interval on summary ROC
#' @param orderByAUC if TRUE, then order legend by summary AUC. Otherwise, use default ordering.
#' @param alphaBetaPlots if TRUE, then draw forest plots of alpha and beta. If false, suppress plotting. 
#' 
#' @author Timothy E. Sweeney
#' @import ggplot2
#' @return Generates a ROC plot for all datasets
#' @export
#' @examples 
#' \dontrun{
#' summaryROCPlot(tinyMetaObject,filterObject = 
#'    tinyMetaObject$filterResults$pValueFDR0.05_es0_nStudies1_looaTRUE_hetero0)
#' }
summaryROCPlot <- function(metaObject, filterObject, bootstrapReps=500, orderByAUC=TRUE, alphaBetaPlots=TRUE) {
  rocObjList <- .summaryROCrocObjList(metaObject, filterObject, bootstrapReps)
  roc.stats <- .metaROC.KestBuntix(rocObjList, bootReps=bootstrapReps, print=F) 
  if(alphaBetaPlots) {
  	.metaROCsummaries(roc.stats, "")
  }
  
  ## Plot summary ROC curve
  summ.roc <- .rocplot.mult.SummROC(test.data.list = rocObjList, summStats = roc.stats,
                                   title = "Summary ROC", size=16, orderByAUC = orderByAUC)
  
  
  print(.clinParamMultiRoc(rocObjList))
  return(summ.roc)
  
}

#' Calculate the summaryROC statistics
#' @param metaObject a Meta object which must have the $originalData populated
#' @param filterObject a MetaFilter object containing the signature genes that will be used for calculating the score
#' @param bootstrapReps number of bootstrap simulations to run for confidence interval on summary ROC
#' 
#' @author Timothy E. Sweeney
#' @import ggplot2
#' @return Summary AUC statistics
#' @export
#' @examples 
#' \dontrun{
#' summaryROCCalc(tinyMetaObject, filterObject = 
#'    tinyMetaObject$filterResults$pValueFDR0.05_es0_nStudies1_looaTRUE_hetero0)
#' }
summaryROCCalc <- function(metaObject, filterObject, bootstrapReps=500) {
  rocObjList <- .summaryROCrocObjList(metaObject, filterObject, bootstrapReps)
  roc.stats <- .metaROC.KestBuntix(rocObjList, bootReps=bootstrapReps, print=F) 
  summaryROC <- .summaryROCcalculate(roc.stats)
  
  return(data.frame(AUC=summaryROC$auc.summ, AUCupper=summaryROC$auc.upper, AUClower=summaryROC$auc.lower))
  
}

.summaryROCrocObjList <- function(metaObject, filterObject, bootstrapReps=500) {
  
  #We have to have at least 5 of each class to include in the summary ROC computation
  whichIndices <- sapply(metaObject$originalData, function(x) { 
    return(length(which(x$class==0))>=5 && length(which(x$class==1))>=5)
  })
  metaObject$originalData <- metaObject$originalData[whichIndices]
  
  if(sum(!whichIndices) > 0) {
    print(paste("Removed",paste(names(whichIndices)[!whichIndices], collapse=", "),
                "from the summary ROC because there are not at least 5 samples of each class"))
  }
  #We need classes with different scores to calculate the ROC
  zeroScoreFunction <- function(x, curFilter)  {
    return(length(which(calculateScore(datasetObject = x, filterObject = curFilter, suppressMessages = TRUE)!=0))>0) 
  }
  whichScores <- sapply(metaObject$originalData,  zeroScoreFunction, curFilter=filterObject)
  metaObject$originalData <- metaObject$originalData[whichScores]  
  
  if(sum(!whichScores) > 0) {
    print(paste("Removed",paste(names(whichScores)[!whichScores], collapse=", "),
                "from the summary ROC because there are not different scores for each sample"))
  }
  rocObjList <- .rocObjList(metaObject, filterObject)
  roc.stats <- .metaROC.KestBuntix(rocObjList, bootReps=bootstrapReps, print=F) 
  
  return(rocObjList)
  
}

###################3#  ROC meta  ######################################33
#### Kester and Buntix method 
#### TES 1/22/2015

.metaROC.KestBuntix <- function(rocDataList, print=F, bootReps=1000,  minSize=5, treatAs1=NULL){
	#preallocate
	comparisons <- names(rocDataList)
	stats <- data.frame("N" = rep(0, length(comparisons)), 
											"tstar_alpha" = 0, "SE_alpha"=0, 
											"tstar_beta"= 0, "SE_beta" =0)
	rownames(stats) <- comparisons
	
	plotdata <- plyr::llply(rocDataList, function(x) with(x, .rocdata(grp = group, pred = score)))
	plotdata <- list(roc = plyr::ldply(plotdata, function(x) x$roc),
									 stats = plyr::ldply(plotdata, function(x) x$stats) )
	
	## allows to treat specific instances as AUC=1. 
	## then alpha = log((P+0.5)*(N+0.5)/0.25), beta = 0
	if(is.null(treatAs1)){
		## if auc is > some thresh (0.99 arbitrarily)
		auc.1 <- which(plotdata$stats$auc > 0.99) 
		## or if no points off axes
		#    noPoints <- which(unlist(lapply(unique(plotdata$roc$.id), function(id){
		#        roc <- plotdata$roc[plotdata$roc$.id==id, ]
		#        all(roc$x==0 | roc$y==1)  
		#    })))
		auc.1 <- c(auc.1)            
	} else { 
		auc.1 <- which(names(rocDataList) %in% treatAs1)
	}
	
	auc.0 <- which(plotdata$stats$auc==0)
	
	for(i in 1:length(comparisons)) {
		comparison <- comparisons[i]   
		cat("\n", comparison, "\t")
		rocobj <- plotdata$roc[plotdata$roc$.id==comparison, c("x", "y")]
		colnames(rocobj) <- c("x1_spec", "sens")
		
		N <- table(rocDataList[[comparison]]$group)[1]
		P <- table(rocDataList[[comparison]]$group)[2]
		
		if(i %in% auc.1){
			cat("treated as AUC=1")
			#as defined in the paper. alpha = log((P+0.5)*(N+0.5)/0.25)
			alpha <- log((P+0.5)*(N+0.5)/0.25)
			#ASE = sqrt of asymptotic variance over sqrt n 
			ASE <- (1/0.5 + 1/0.5 + 1/(P+0.5) + 1/(N+0.5))^(-1/2)/sqrt(P+N)
			stats[comparison, ] <- c(N+P, alpha, ASE , 0, 10)
			next
		} else if(i %in% auc.0) {
			#added for edge case of auc=0
			alpha <- -log((P+0.5)*(N+0.5)/0.25)
			ASE <- (1/0.5 + 1/0.5 + 1/(P+0.5) + 1/(N+0.5))^(-1/2)/sqrt(P+N)
			stats[comparison, ] <- c(N+P, alpha, ASE , 0, 10)
			next
		} else {
			## added 12/30/15-- better estimates in few-point cases; however, not in original paper.
			##  Also, means always 3 points, so increased minSize default
			rocobj$sens[min(which(rocobj$sens==1))] <- 0.98
			rocobj$x1_spec[max(which(rocobj$x1_spec==0))] <- 0.02
			
			# weights defined in paper appendix, sends inf/NAN to 0
			getWeight <- function(sens, x1_spec, N, P){
				TP <- sens * P 
				TN <- (1-x1_spec) * N 
				FP <- P - TP
				FN <- N - TN
				weight <- (1/TP + 1/TN + 1/FP + 1/FN)^(-1)
				weight[is.nan(weight) | is.na(weight)] <- 0 
				weight 
			}
			
			rocobj$weights <- getWeight(sens=rocobj$sens, x1_spec=rocobj$x1_spec, N=N, P=P)
			
			### As defined in paper. Note x1_spec = 1-spec, so the fraction looks upside-down
			rocobj$S <- with(rocobj, log(sens/(1-sens)) + log(x1_spec/(1-x1_spec)) )
			rocobj$D <- with(rocobj, log(sens/(1-sens)) - log(x1_spec/(1-x1_spec)) )
			
			
			### this method changed to 'points' as below 12/8/15
			### if there are enough non-0 weights, then bootstrap
			points <- sum(rocobj$weights != 0) 
			### need points that vary in x and y, not just non-0 weights
			# points <- sum(!duplicated(rocobj$sens) & !duplicated(rocobj$x1_spec))
			
			if(points > minSize){
				cat("\tBootstrapping...\t")
				stats[comparison, ] <- .KestBuntBoot(rocobj, N, P, print, name=comparison, reps=bootReps)
				next
				## else assign ASE as errors for both alpha and beta
			} else if (points >= 2){
				cat("\tCan't bootstrap, points < minSize...\t")
				model <- stats::lm(D~S, data=rocobj, weights=weights)
				alpha <- model$coefficients[1]
				beta <- 0
				ASE <- (1/0.5 + 1/0.5 + 1/(P+0.5) + 1/(N+0.5))^(-1/2)/sqrt(P+N)            
				stats[comparison, ] <- c(N+P, alpha, ASE , beta, ASE)
				next
				## if only one point in ROC that is not on axis, can't use.
			} else {
				cat(sprintf(" has <2 points in ROC curve; cannot compute summary stats\n", comparison))
				stats[comparison, ] <- c(N+P, rep(NA, 4))
				next
			}
		}
	}
	cat("\n")
	print(stats)
	return(stats) 
}  

.KestBuntBoot <- function(rocobj, N, P, print, name, reps=10000) {
	
	#roc.lm <- lm(D~S, data=rocobj, weights=weights)
	#bootstrap coefficients
	bs <- function(formula, data, indices){
		d <- data[indices, ]
		fit <- stats::lm(formula, data=d, weights=weights)
		return(stats::coef(fit))
	}
	boot.out <- boot::boot(data=rocobj, statistic=bs, formula=D~S, R=reps, weights=rocobj$weights)
	if(print) print(boot.out)
	## from  boot:::print.boot():
	op <- NULL
	for (i in 1:2) op <- rbind(op, boot::imp.moments(boot.out, index = i)$rat)
	std.error <- sqrt(op[, 2L])
	
	##plot V, U, S, D if desired
	if(print){
	  print(with(rocobj, plot(log(x1_spec/(1-x1_spec)), log(sens/(1-sens)), 
	                          xlab="U", ylab="V", 
	                          main=paste0(name, "\n logit of ROC") )  ))
	  print(plot(rocobj$S, rocobj$D, xlab="V + U", ylab="V - U", 
	             main=paste0(name, "\n regression of D on S")), 
	        graphics::abline(a=op[1,1], b=op[2,1]) )
	}  
	
	##Bias-corrected bootstrap estimates:
	c(N+P, op[1,1], std.error[1], op[2,1], std.error[2]) 
}

.metaROCsummaries <- function(stats, title, method="random"){
	alpha <- with(stats, rmeta::meta.summaries(d=tstar_alpha, se=SE_alpha*sqrt(N), method=method))
	print(alpha)
	alpha$names <- rownames(stats)
	rmeta::metaplot(alpha$effects, alpha$stderrs, main=paste0("Summary of Alpha\n", title), labels = alpha$names,
	                colors = rmeta::meta.colors(box = "blue",lines = "lightblue",zero = "black", summary = "orange", text = "black"))
	
	beta <- with(stats, rmeta::meta.summaries(d=tstar_beta, se=SE_beta*sqrt(N), method=method))
	print(beta)
	beta$names <- rownames(stats)
	rmeta::metaplot(beta$effects, beta$stderrs, main=paste0("Summary of Beta\n", title), labels = beta$names,
	                colors = rmeta::meta.colors(box = "blue",lines = "lightblue",zero = "black", summary = "orange", text = "black"))
	# plot(.newMetaROC(alpha=alpha$summary, beta=beta$summary), type='l', main=title)
	# lines(.newMetaROC(alpha=alpha$summary + 1.96*alpha$se.summary, beta=beta$summary + 1.96*beta$se.summary), type='l', lty=2)
	# lines(.newMetaROC(alpha=alpha$summary - 1.96*alpha$se.summary, beta=beta$summary - 1.96*beta$se.summary), type='l', lty=2)
}

###-###-###-###-###-###-##
###   .aucROCframe()   ###
###-###-###-###-###-###-##

#DESCRIPTION
#Extracts the AUC from a matrix of FPR and TPR values

#PARAMETERS
#newroc - matrix where first column is FPR and second column is TPR

#RETURN VALUE
#the AUC value

#REQUIRED PACKAGES: None

.aucROCframe <- function(newroc){
  i <- 2:dim(newroc)[1]
  auc <- (newroc[i-1, "FPR"] - newroc[i, "FPR"]) %*% (newroc[i-1, "TPR"] + newroc[i, "TPR"])/2
  return(auc)
}


##-###-###-###-###-###-##
###   .newMetaROC()   ###
##-###-###-###-###-###-##

#DESCRIPTION
#Using the alpha and beta parameters from .getMetaROCStats(), this function makes
#a matrix of FPR and TPR values for the summary ROC curve

#PARAMETERS
#alpha - alpha parameter
#beta - beta parameter
#points - how many points to include for the fpr/tpr
#print - whether to print the summary AUC

#RETURN VALUE
#a matrix of FPR and TPR values

#REQUIRED PACKAGES: None

.newMetaROC <- function(alpha, beta, points=1000, print=F){
  ## in practice, beta > 0.95 becomes unbounded
  if(abs(beta)>0.95) beta <- 0.95 * sign(beta)
  A <- alpha/(1-beta)
  B <- (1+beta)/(1-beta)
  
  by <- 1/points
  newspec <- seq(by, 1-by, by)
  newsens <- exp(A+B*log((1-newspec)/newspec)) / (1 + exp(A+B*log((1-newspec)/newspec)) )
  newroc <- data.frame("FPR" = 1-newspec, "TPR" = newsens)
  
  auc <- .aucROCframe(newroc)
  if(print){
    cat("AUC: ", auc, "\n")
  }
  return(newroc)
}

.summaryROCcalculate <- function(summStats, method="random") {
  
  ## Find summary ROC parameters, remove any missing rows
  summStats <- stats::na.omit(summStats)
  ## Very low standard errors leads to NA values in the call to meta.summaries, so round them up to 1e-8
  summStats$SE_alpha[(summStats$SE_alpha *sqrt(summStats$N)) < 1e-8 ] <- 1e-8
  summStats$SE_beta[(summStats$SE_beta *sqrt(summStats$N)) < 1e-8 ] <- 1e-8
  if(nrow(summStats)>1){
    alpha <- with(summStats, rmeta::meta.summaries(d=tstar_alpha, se=SE_alpha*sqrt(N), method=method))
    beta <- with(summStats, rmeta::meta.summaries(d=tstar_beta, se=SE_beta*sqrt(N), method=method))
  } else {
    alpha <- with(summStats, data.frame(summary=tstar_alpha, se.summary=SE_alpha))
    beta <- with(summStats, data.frame(summary=tstar_beta, se.summary=SE_beta))
  }
  
  roc.summ <- .newMetaROC(alpha=alpha$summary, beta=beta$summary)
  roc.upper.b.hi <- .newMetaROC(alpha=alpha$summary + 1.96*alpha$se.summary, beta=beta$summary + 1.96*beta$se.summary, print=F)
  roc.upper.b.lo <- .newMetaROC(alpha=alpha$summary + 1.96*alpha$se.summary, beta=beta$summary - 1.96*beta$se.summary, print=F)
  roc.lower.b.hi <- .newMetaROC(alpha=alpha$summary - 1.96*alpha$se.summary, beta=beta$summary + 1.96*beta$se.summary, print=F)
  roc.lower.b.lo <- .newMetaROC(alpha=alpha$summary - 1.96*alpha$se.summary, beta=beta$summary - 1.96*beta$se.summary, print=F)
  
  roc.upper <- data.frame(FPR = roc.upper.b.hi$FPR, TPR = pmax(roc.upper.b.hi$TPR, roc.upper.b.lo$TPR, na.rm=T) )
  roc.lower <- data.frame(FPR = roc.lower.b.hi$FPR, TPR = pmin(roc.lower.b.hi$TPR, roc.lower.b.lo$TPR, na.rm=T) )
  
  auc.summ <- .aucROCframe(roc.summ)
  auc.upper <- .aucROCframe(roc.upper)
  auc.lower <- .aucROCframe(roc.lower)
  
  return(list(roc.summ=roc.summ, roc.upper=roc.upper, roc.lower=roc.lower, 
              auc.summ=auc.summ, auc.upper=auc.upper, auc.lower=auc.lower))
}

.rocplot.mult.SummROC <- function(test.data.list, summStats, groupName = "group", predName = "score", 
                                 title = "ROC Plot", size=16, method="random", printCI=T, orderByAUC=TRUE) {
	## Standard ROC stuff
	plotdata <- plyr::llply(test.data.list, function(x) with(x, .rocdata(grp = eval(parse(text = groupName)), pred = eval(parse(text = predName)))))
	names(plotdata) <- names(test.data.list)
	
	if(orderByAUC) {
		aucVals <- sapply(plotdata, function(x) x$stats$auc[[1]])
		plotdata <- plotdata[names(plotdata)[order(aucVals, decreasing=T)]]
	}
	plotdata <- list(roc = plyr::ldply(plotdata, function(x) x$roc),
									 stats = plyr::ldply(plotdata, function(x) x$stats))
	
	if(printCI){
		annotation <- with(plotdata$stats, paste0(plotdata$stats$.id, " AUC=", signif(auc, 2), " (95% CI ", signif(ci.lower, 2), " - ", signif(ci.upper, 2), ")"))
	} else {
		annotation <- with(plotdata$stats, paste0(plotdata$stats$.id, " AUC=", signif(auc, 2)))
	}
	
	sumRoc <- .summaryROCcalculate(summStats, method=method)
	
	## Add summary curve to labels, data list. 
	annotation <- c(annotation, paste0("Summary AUC=", signif(sumRoc$auc.summ, 2), " (95% CI ", signif(sumRoc$auc.lower, 2), " - ", signif(sumRoc$auc.upper, 2), ")") )
	names(annotation) <- c(plotdata$stats$.id, "Summary")
	
	colnames(sumRoc$roc.summ) <- c("x", "y")
	plotdata$roc <- rbind(plotdata$roc, cbind(.id="Summary", sumRoc$roc.summ))
	plotdata$roc$.id <- ordered(plotdata$roc$.id, levels=names(annotation))
	
	#force color to dark grey on palette
	gg_color_hue <- function(n){
		hues <- seq(15, 375, length=n+1)
		grDevices::hcl(h=hues, l=65, c=100)[1:n]
	}
	dataPal <- gg_color_hue(length(plotdata$stats$.id))
	mypal <- c(dataPal, "grey25")
	
	p <- ggplot(plotdata$roc, aes_string(x = 'x', y = 'y')) +
		#Have this first to put it on bottom
		geom_ribbon(data=cbind(summ=sumRoc$roc.summ, lower=sumRoc$roc.lower, upper=sumRoc$roc.upper), 
								aes_string(x= 'lower.FPR', y='summ.y', ymin='lower.TPR', ymax='upper.TPR'), 
								color="gray25", fill="gray75") +
		geom_line(aes_string(colour = '.id')) +
		geom_segment(x=0, y=0, xend=1, yend=1, linetype=5, color="grey20") +
		scale_x_continuous("False Positive Rate (1-Specificity)") +
		scale_y_continuous("True Positive Rate (Sensitivity)") +
		scale_color_manual(values=mypal, labels = annotation) +
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
		geom_line(data=sumRoc$roc.summ, aes_string(x='x', y='y'), size=1.3, color="gray20") 
	## If just want to define upper and lower lines without fill:
	#geom_line(data=roc.upper, aes(x=FPR, y=TPR), size=1.1, linetype=3, color="gray15") +
	#geom_line(data=roc.lower, aes(x=FPR, y=TPR), size=1.1, linetype=3, color="gray15") 
	
	cat(title, "\n")
	.youdenJ(sumRoc$roc.summ)
	
	return(p)
}

.youdenJ <- function(rocframe){
	#assume col 1 = FPR, col 2 = TPR
	spec <- (1 - rocframe[,1])
	sens <- rocframe[,2]
	tmp <- spec+sens-1
	maxsens <- signif(sens[which.max(tmp)]*100, 4)
	maxspec <- signif(spec[which.max(tmp)]*100, 4)
	cat(sprintf("At optimal cutoff, summary sensitivity = %s%%, summary specificity = %s%%", maxsens, maxspec))
	
	
	spec95 <- spec[which.min(abs(sens-0.95))]
	
	cat(sprintf("\nAt summary sensitivity = 95%%, summary specificity = %s%%", spec95))
}

.clinParamMultiRoc <- function(test.data.list){
	## assumes $score and $group columns
	out <- vector()
	
	for(i in 1:length(test.data.list)){
		data <- test.data.list[[i]]
		score.1 <- subset(data, group==min(group))$score
		score.2 <- subset(data, group!=min(group))$score
		thresh <- (mean(score.1)/stats::sd(score.1) +  mean(score.2)/stats::sd(score.2))/(1/stats::sd(score.1) + 1/stats::sd(score.2))
		
		TN <- with(data, sum(score<thresh & group==min(group)))
		FN <- with(data, sum(score<thresh & group!=min(group)))
		TP <- with(data, sum(score>=thresh & group!=min(group)))
		FP <- with(data, sum(score>=thresh & group==min(group)))
		
		sens <- TP/length(score.2)
		spec <- TN/length(score.1)
		PPV <- TP/(TP+FP)
		NPV <- TN/(TN+FN)
		acc <- (TP+TN)/(TP+TN+FP+FN)
		
		tmp <- c(sens, spec, PPV, NPV, acc)
		out <- rbind(out, tmp)
	}
	
	out <- apply(out, 2, signif, digits=3)
	
	colnames(out) <- c("sens", "spec", "PPV", "NPV", "acc")
	rownames(out) <- names(test.data.list)
	return(out)
}


#declare global variables for variables in data.table/with notation to avoid R CMD CHECK notes
utils::globalVariables(c("group","weights"))
