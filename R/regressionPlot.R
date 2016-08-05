#This file was written by Winn Haynes
#Last updated: May 12, 2015

#'Generate a plot which draws a regrssion line between the Meta Score and a continuous variable phenotype.
#'
#'@param filterObject a MetaFilter object containing the signature genes that will be used for the z-score calculation
#'@param datasetObject a Dataset object (typically independent validation dataset) for comparison in a regression plot
#'@param continuousVariableColumn the label of the column in $pheno that sepecifies the continuous variable to compare (default: 'continuousVariableColumn')
#'@param formattedVariableName label which will be used on the x-axis on the plot
#'
#'@return Returns a regression plot as ggplot2 plot object
#'
#'@author Winston A. Haynes
#'
#'@keywords graph
#'
#'@import ggplot2
regressionPlot <- function(filterObject, datasetObject, continuousVariableColumn= "continuous", formattedVariableName="Continuous Variable") {
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
	
	return(.createRegressionPlot(datasetPheno, formattedName=datasetObject$formattedName, formattedVariableName=formattedVariableName))
}

.createRegressionPlot <- function(datasetPheno, comparisonColumn="continuous", scoreColumn="score", formattedName="",
																	formattedVariableName="Continuous Variable") {
	
	#Set visualization parameters
	segmentWidth <- 0.05
	jitterWidth <- 0.025
	
	#Swap our column names to "x" and "y", so that we can work with them in the ggplot drawing
	datasetPheno$y <- datasetPheno[,scoreColumn]
	datasetPheno$x <- datasetPheno[,comparisonColumn]
	
	plotObject <- ggplot(datasetPheno, aes_string(x="x",y="y")) + 
		geom_point() +
		xlab("SLEDAI")  + 
		ylab("Meta Score")+
		ggtitle(formattedName) +
	  theme(text = element_text(size=24)) +
		annotate("text",x=0.9*max(datasetPheno$x, na.rm=TRUE), y=1, size=6,
             label=paste("R=",signif(cor(datasetPheno$x,datasetPheno$y, use="complete.obs"),3))) + 
		geom_smooth(method="lm") +
		xlab(formattedVariableName)

	return(plotObject)
}
