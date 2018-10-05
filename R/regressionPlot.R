# This file was written by Winn Haynes
# Last updated: January 17, 2017

#'Generate a plot which draws a regression line between the Meta Score and a continuous variable phenotype.
#'
#'@param filterObject a MetaFilter object containing the signature genes that will be used for the z-score calculation
#'@param datasetObject a Dataset object (typically independent validation dataset) for comparison in a regression plot
#'@param continuousVariableColumn the label of the column in $pheno that specifies the continuous variable to compare (default: 'continuousVariableColumn')
#'@param formattedVariableName label which will be used on the x-axis on the plot
#'@param corMethod method which will be passed to cor.test
#'@param correlationCorner one of topLeft, topRight, bottomLeft, bottomRight (default: bottomRight)
#'
#'@return Returns a regression plot as ggplot2 plot object
#'
#'@author Winston A. Haynes
#'
#'@keywords graph
#'@examples
#' regressionPlot(tinyMetaObject$filterResults[[1]], 
#'                tinyMetaObject$originalData$Whole.Blood.Study.1,
#'                continuousVariableColumn="age",
#'                formattedVariableName="Age")
#'@import ggplot2
#'@export
regressionPlot <- function(filterObject, datasetObject, continuousVariableColumn= "continuous", 
                           formattedVariableName="Continuous Variable", corMethod="pearson", 
                           correlationCorner="bottomRight") {
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
	if(!continuousVariableColumn %in% colnames(datasetPheno)) {
		stop(paste("Column named \'",continuousVariableColumn, "\' not present in the datasetObject$pheno data frame", sep=""))
	}
	
	#Calculate the score
	datasetPheno$score <- calculateScore(filterObject, datasetObject, suppressMessages = TRUE)
	
	#Populate the continuous column
	datasetPheno$continuous <- datasetPheno[,continuousVariableColumn]
	
	return(.createRegressionPlot(datasetPheno, formattedName=datasetObject$formattedName, 
	                             formattedVariableName=formattedVariableName, corMethod=corMethod, 
	                             correlationCorner=correlationCorner))
}

.createRegressionPlot <- function(datasetPheno, comparisonColumn="continuous", scoreColumn="score", formattedName="",
																	formattedVariableName="Continuous Variable", corMethod="pearson", 
																	correlationCorner="bottomRight") {
	
	#Set visualization parameters
	segmentWidth <- 0.05
	jitterWidth <- 0.025
	
	#Swap our column names to "x" and "y", so that we can work with them in the ggplot drawing
	datasetPheno$y <- datasetPheno[,scoreColumn]
	datasetPheno$x <- datasetPheno[,comparisonColumn]

	corRes <- stats::cor.test(datasetPheno$x,datasetPheno$y, use="complete.obs", method=corMethod)
	
	coordVal <- c(Inf, -Inf, 1, -0.1)
	if(correlationCorner == "bottomLeft") { coordVal <- c(-Inf, -Inf, 0, -0.1) }
	if(correlationCorner == "topLeft") { coordVal <- c(-Inf, Inf, 0, 1) }
	if(correlationCorner == "topRight") { coordVal <- c(Inf, Inf, 1, 1) }
	
	
	plotObject <- ggplot(datasetPheno, aes_string(x="x",y="y")) + 
		geom_point() +
		xlab("SLEDAI")  + 
		ylab("Meta Score")+
		ggtitle(formattedName) +
	  theme(text = element_text(size=24)) +
		annotate("text",x=coordVal[1], y=coordVal[2], hjust=coordVal[3], vjust=coordVal[4], size=6,
             label=paste("R=",signif(corRes$estimate,3), "\np=", signif(corRes$p.value,3))) + 
		geom_smooth(method="lm") +
		xlab(formattedVariableName)

	return(plotObject)
}
