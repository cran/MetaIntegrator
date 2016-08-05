\name{violinPlot}
\alias{violinPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compare groups within a single dataset in a violin plot  
}
\description{
Given a \code{filterObject} and a \code{datasetObject} this function will use the selected genes of the \code{filterObject} to calculate and compare the z-scores of the groups (e.g. cases vs. controls) from the \code{datasetObject} by generating a violin plot. A violin plot is similar to a box plot, except the width of each violin is proportional to the density of points. \code{violinPlot()} is commonly used to validate a gene signature in an independent dataset.
}
\usage{
violinPlot(filterObject, datasetObject, labelColumn="label")
}
\arguments{
  \item{filterObject}{a MetaFilter object containing the signature genes that will be used for the z-score calculation}
  \item{datasetObject}{a Dataset object (typically independent validation dataset) for group comparison in a violin plot}
  \item{labelColumn}{the label of the column in \code{$pheno} that sepecifies the groups to compare, typically case or control (default: 'label')}
}

\details{
The z-score is based off of the geometric mean of expression. As such, negative expression values are not allowed. A dataset is thus always scaled by its minimum value + 1, such that the lowest value = 1. Any individual NANs or NAs are also set to 1. If a dataset does not have any information on a given gene, the entire gene is simply left out of the score. 
}
\value{
Returns a violin plot as ggplot2 plot object

}
\author{
	Winston A. Haynes
}


\seealso{
 \code{\link{filterGenes}},  \code{\link{runMetaAnalysis}}
}
\examples{
violinPlot(tinyMetaObject$filterResults$pValueFDR0.05_es0_nStudies1_looaTRUE_hetero0, 
		tinyMetaObject$originalData$Whole.Blood.Study.1, labelColumn="group")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }   
\keyword{ graphs }
