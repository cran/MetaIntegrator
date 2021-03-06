% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/regressionPlot.R
\name{regressionPlot}
\alias{regressionPlot}
\title{Generate a plot which draws a regression line between the Meta Score and a continuous variable phenotype.}
\usage{
regressionPlot(filterObject, datasetObject,
  continuousVariableColumn = "continuous",
  formattedVariableName = "Continuous Variable", corMethod = "pearson",
  correlationCorner = "bottomRight")
}
\arguments{
\item{filterObject}{a MetaFilter object containing the signature genes that will be used for the z-score calculation}

\item{datasetObject}{a Dataset object (typically independent validation dataset) for comparison in a regression plot}

\item{continuousVariableColumn}{the label of the column in $pheno that specifies the continuous variable to compare (default: 'continuousVariableColumn')}

\item{formattedVariableName}{label which will be used on the x-axis on the plot}

\item{corMethod}{method which will be passed to cor.test}

\item{correlationCorner}{one of topLeft, topRight, bottomLeft, bottomRight (default: bottomRight)}
}
\value{
Returns a regression plot as ggplot2 plot object
}
\description{
Generate a plot which draws a regression line between the Meta Score and a continuous variable phenotype.
}
\examples{
regressionPlot(tinyMetaObject$filterResults[[1]], 
               tinyMetaObject$originalData$Whole.Blood.Study.1,
               continuousVariableColumn="age",
               formattedVariableName="Age")
}
\author{
Winston A. Haynes
}
\keyword{graph}
