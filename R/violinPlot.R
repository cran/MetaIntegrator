#' Compare groups within a single dataset in a violin plot  
#'@description
#'	Given a \code{filterObject} and a \code{datasetObject} this function will use the selected genes of the \code{filterObject} to calculate and compare the z-scores of the groups (e.g. cases vs. controls) from the \code{datasetObject} by generating a violin plot. A violin plot is similar to a box plot, except the width of each violin is proportional to the density of points. \code{violinPlot()} is commonly used to validate a gene signature in an independent dataset.
#'@param filterObject a MetaFilter object containing the signature genes that will be used for the z-score calculation
#'@param datasetObject a Dataset object (typically independent validation dataset) for group comparison in a violin plot
#'@param labelColumn the label of the column in \code{$pheno} that specifies the groups to compare, typically case or control (default: 'label')
#'@param comparisonMethod statistical test that will be used (default="wilcox.test"). Other options include "t.test". 
#'@param pairwiseComparisons if TRUE, perform pairwise statistical comparisons against the first factor level. If FALSE, perform global statistical comparisons (default: TRUE).
#'@param autoLineBreak if TRUE, insert line breaks into labels on plots. If FALSE, don't insert line breaks (default: TRUE)
#'@details
#'	The z-score is based off of the geometric mean of expression. As such, negative expression values are not allowed. A dataset is thus always scaled by its minimum value + 1, such that the lowest value = 1. Any individual NANs or NAs are also set to 1. If a dataset does not have any information on a given gene, the entire gene is simply left out of the score. 
#'@return	Returns a violin plot as ggplot2 plot object
#'@author Winston A. Haynes
#'@import ggpubr ggplot2
#'@seealso
#'	\code{\link{filterGenes}},  \code{\link{runMetaAnalysis}}
#'@examples
#'	violinPlot(tinyMetaObject$filterResults$pValueFDR0.05_es0_nStudies1_looaTRUE_hetero0, 
#'	   tinyMetaObject$originalData$Whole.Blood.Study.1, 
#'	     labelColumn="group")
#'@keywords 
#'hplot 
#'graphs 
#'@export

violinPlot <- function(filterObject, datasetObject, labelColumn= "label", comparisonMethod="wilcox.test", 
                       pairwiseComparisons=TRUE, autoLineBreak=TRUE) {
	#Load relevant libraries
	
	#Check that our objects are the right type
	if(! checkDataObject(object = filterObject, objectType="MetaFilter")) {
		stop("filterObject that was passed to violinPlot was not appropriately formatted as a MetaFilter object")
	} 
	if(!checkDataObject(object =datasetObject, objectType =  "Dataset")){
		stop("datasetObject that was passed to violinPlot was not appropriately formatted as a Dataset object")
	}
	
	#Populate a pheno object
	datasetPheno <- datasetObject$pheno
	
	#Make sure the label column is present
	if(!labelColumn %in% colnames(datasetPheno)) {
		stop(paste("Column named \'",labelColumn, "\' not present in the datasetObject$pheno data frame", sep=""))
	}

	#Calculate the score
	datasetPheno$score <- calculateScore(filterObject, datasetObject, suppressMessages = TRUE)

	return(.createViolinPlot(datasetPheno, formattedName=datasetObject$formattedName, comparisonColumn=labelColumn, 
	                         pairwiseComparisons= pairwiseComparisons, autoLineBreak=autoLineBreak, comparisonMethod=comparisonMethod))
}

#Set some visual parameters. For some reason these have to be outside the function call....

#Draw the violin plot with:
#	standard error bars
#	Wilcox rank sum comparison between groups

#segmentWidth=0.05
#jitterWidth=0.025

.createViolinPlot <- function(datasetPheno, comparisonColumn="label", scoreColumn="score", formattedName="", 
                              comparisonMethod="wilcox.test", pairwiseComparisons=TRUE, autoLineBreak=TRUE) {
	
	#Set visualization parameters
	segmentWidth <- 0.05
	jitterWidth <- 0.025

	#Swap our column names to "x" and "y", so that we can work with them in the ggplot drawing
	datasetPheno$y <- datasetPheno[,scoreColumn]
	datasetPheno$x <- datasetPheno[,comparisonColumn]

  #Function to split names at spaces occurring every ~10 characters into new lines
  strSplitFun <- function(text) {
    prevPiece<- ""
    curPieces <- c()
    maxSplitLength <- 10
    for(textPiece in strsplit(text," ")[[1]]) {
      if(nchar(textPiece)+nchar(prevPiece) < maxSplitLength) {
        prevPiece <- paste(prevPiece, textPiece, sep=" ")
      } else {
        if(! prevPiece=="") {
          curPieces <- c(curPieces, prevPiece)
        }
        prevPiece <- textPiece
      }
    }
    curPieces <- c(curPieces, prevPiece)
    paste(curPieces, collapse="\n")
  }
  
	
	#Capitalize first letter for formatting
	simpleCap <- function(x) {
	  s <- strsplit(x, " ")[[1]]
	  paste(toupper(substring(s, 1,1)), substring(s, 2),
	        sep="", collapse=" ")
	}
	
	#Preserve the factor ordering
  if(class(datasetPheno$x) =="factor") {
  	factorLevels <- levels(datasetPheno$x)
	
		#Make sure that only levels which have contents are used
		factorLevels <- factorLevels[factorLevels %in% as.character(datasetPheno$x)]

  	factorLevels <- sapply(factorLevels, simpleCap)
  	if(autoLineBreak) {
  	  factorLevels <- sapply(factorLevels, strSplitFun)
  	}
  	
  	#If capitalization has lead to duplicated factors, remove the duplicates
  	factorLevels <- factorLevels[!duplicated(factorLevels)]
  }
  
  #Start off with datasetPheno as a character so we can split text
  datasetPheno$x <- as.character(datasetPheno$x)
  datasetPheno$x <- sapply(datasetPheno$x, simpleCap)
  
	#Turn into a factor so that we can use group reasoning
	if(exists("factorLevels")) {
	  if(autoLineBreak) {
	    datasetPheno$x <- factor(sapply(datasetPheno$x, strSplitFun), levels=factorLevels)
	  } else {
	    datasetPheno$x <- factor(datasetPheno$x, levels=factorLevels)
	  }
	} else {
	  if(autoLineBreak) {
	    datasetPheno$x <- as.factor(sapply(datasetPheno$x, strSplitFun))
	  } else {
	    datasetPheno$x <- as.factor(datasetPheno$x)
	  }
	}
  
	#Calculate standard errors
	standardError = Rmisc::summarySE(data.frame(x=datasetPheno$x, y=datasetPheno$y), measurevar="y", groupvars="x")[, c("x", "y", "se")]
	standardError$segmentWidth <- rep(segmentWidth, nrow(standardError))

	#Hack to supress NOTES when running CRAN checks:
	x <- NULL; y <- NULL; se <- NULL; label <- NULL; group <- NULL
	
	#Set up the comparisons for p-values
	comparisons <- split(matrix(c(rep(levels(datasetPheno$x)[1],(length(levels(datasetPheno$x))-1)),
	                              levels(datasetPheno$x)[2:(length(levels(datasetPheno$x)))]),ncol=2, byrow = F), 1:(length(levels(datasetPheno$x))-1))

	standardError$x <- match(standardError$x,levels(standardError$x))
	standardError$xPosSeg <- standardError$x + standardError$segmentWidth
	standardError$xNegSeg <- standardError$x - standardError$segmentWidth
	standardError$yPosSe <- standardError$y + standardError$se
	standardError$yNegSe <- standardError$y - standardError$se
	
	#Draw the ggplot
	plotObject <- ggplot(datasetPheno, aes_string('x','y')) + geom_violin(fill='grey',trim=FALSE) +
		geom_jitter(aes_string(color='x'), size = 3, position=position_jitter(width=jitterWidth)) +
		xlab(simpleCap(comparisonColumn)) +
		ylab("Meta Score") +
		theme(text = element_text(size=24)) +
		theme(legend.title=element_blank(), legend.text=element_text(size=24)) +
		theme(axis.text.x = element_text(size=24), axis.text.y = element_text(size=24)) +
		geom_segment(data=standardError, aes_string(x='xNegSeg',
								 		xend='xPosSeg', y='yNegSe', yend='yNegSe'), col='black') +
		geom_segment(data=standardError,  aes_string(x='xNegSeg',
								 		xend='xPosSeg', y='yPosSe', yend='yPosSe'), col='black') +
		geom_segment(data=standardError, aes_string(x='x',
								 		xend='x', y='yPosSe', yend='yNegSe'), col='black') +
		geom_point(data=standardError, aes(x=x, y=y), color="black", size=3) +
		ggtitle(formattedName) + guides(fill=FALSE, color=FALSE)
	if(! pairwiseComparisons) {
	  plotObject <- plotObject + stat_compare_means(label = "p.signif", method = comparisonMethod, ref.group = ".all.")
	} else {
	  plotObject <- plotObject + stat_compare_means(comparisons=comparisons, method=comparisonMethod)
	}
	return(plotObject)
}
