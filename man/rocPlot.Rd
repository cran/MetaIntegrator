\name{rocPlot}
\alias{rocPlot}

\title{
Plot ROC Curve for a Dataset
}
\description{
rocPlot will plot an ROC curve (and return the AUC) that describes how well a gene signature (as defined in a \code{filterObject}) classifies groups in a dataset (in the form of a \code{datasetObject}).
}
\usage{
rocPlot(filterObject, datasetObject, title = datasetObject$formattedName)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filterObject}{
  a MetaFilter object containing the signature genes that will be used for calculation of the ROC plot.
}
  \item{datasetObject}{
a Dataset object for group comparison in the ROC plot. (At least, must have a \code{$expr} of probe-level data, \code{$keys} of probe:gene mappings, and \code{$class} of two-class labels.)
}
  \item{title}{
Title for the ROC plot.
}
}
\details{
Evaluates the ability of a given gene set to separate two classes. The gene set is evaluated as a Z-score of the difference in means between the positive genes and the negative genes (see calculateScore). Returns a standard ROC plot, plus AUC with 95\% CI (calculated according to Hanley method).
}
\value{
Returns a ggplot2 plot object
}
\author{
	Timothy E. Sweeney
}
\seealso{
\code{\link{calculateScore}}, \code{\link{calculateROC}}}
\examples{

rocPlot(tinyMetaObject$filterResults[[1]], tinyMetaObject$originalData[[1]])

}

\keyword{ graphs }

