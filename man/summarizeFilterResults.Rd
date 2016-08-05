\name{summarizeFilterResults}
\alias{summarizeFilterResults}
\title{
Summarize the filtered analysis results
}
\description{
Given a  \code{metaObject} and the name of the \code{filterObject} of interest, this function will print a summary style message about genes that passed the filtering step using the function \code{filterGenes()} and return a \code{dataFrame} that contains the \code{$pooledResults} information for each gene which passed the filter.
}
\usage{
summarizeFilterResults(metaObject, metaFilterLabel)
}
\arguments{
  \item{metaObject}{the metaObject that contains the \code{filterObject} of interest}
  \item{metaFilterLabel}{the name of a \code{filterObject} generated with the function \code{filterGenes()}}
}
\value{
  \item{dataFrame}{Data frame, which contains \code{$pooledResults} information for each gene which passed the filter}
}
\author{
	Francesco Vallania
}
\seealso{

\code{\link{filterGenes}}
}
\examples{
# filter genes with default settings 
#		false discovery rate cutoff of 5 percent and WITH leave-one-out analysis
testMetaObject <- filterGenes(tinyMetaObject)
summarizeFilterResults(testMetaObject, getMostRecentFilter(testMetaObject))
}
\keyword{ utilities }
\keyword{ methods }
