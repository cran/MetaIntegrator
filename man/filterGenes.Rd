\name{filterGenes}
\alias{filterGenes}
\title{
Filter out significant genes from meta-analysis results
}
\description{
After the Meta-Analysis results have been written to the \code{metaObject}, the results can be examined using different gene filtering criteria. This function will use the given filterParameter to select genes that fulfill the filter conditions. The function returns a modified version of the \code{metaObject} with results stored in \code{metaObject$filterResults}
    
}
\usage{
filterGenes(metaObject, isLeaveOneOut=TRUE, effectSizeThresh=0, FDRThresh=0.05,  
				numberStudiesThresh=1, heterogeneityPvalThresh=0)}

\arguments{
  \item{metaObject}{a Meta object which must have the \code{$originalData}, \code{$metaAnalysis} populated}
optional filterParameter:
  \item{isLeaveOneOut}{Do leave-one-out analysis on discovery datasets (default: TRUE). Needs at least 2 datasets for discovery.}
  \item{FDRThresh}{FDR cutoff: a gene is selected, if it has a p-value less than or equal to the FDR cutoff (default: 0.05)}
  \item{effectSizeThresh}{a gene is selected, if the absolute value of its effect size is above this threashold (default: 0)}
  \item{numberStudiesThresh}{number of studies in which a selected gene has to be significantly up/down regulated (default: 1)}
  \item{heterogeneityPvalThresh}{heterogeneity p-value cutoff (filter is off by default: \code{heterogeneityPvalThresh = 0}). Genes with significant heterogeneity and, thus a significant (low) heterogeneity p-value, can be filtered out by using e.g.: \code{heterogeneityPvalThresh = 0.05} (removes all genes with heterogeneity p-value < 0.05)}

}
\value{
  \item{metaObject}{A modified version of the input metaObject with an additional filterObject stored within \code{metaObject$filterResults}}
}
\author{
	Francesco Vallania
}
\note{
Use \code{checkDataObject(metaObject, "Meta", "Pre-Filter")} to make sure your metaObject has the right format for filtering after running the meta-analysis with \code{runMetaAnalysis().}

}

\seealso{
\code{\link{checkDataObject}}
}

\examples{
# filter genes with default settings 
#(false discovery rate cutoff of 5 percent and WITH leave-one-out analysis)
testMetaObject <- filterGenes(tinyMetaObject)
summarizeFilterResults(testMetaObject, getMostRecentFilter(testMetaObject))

# filter genes with false discovery rate of 1 percent and WITHOUT leave-one-out analysis  
testMetaObject <- filterGenes(testMetaObject, FDRThresh = 0.01, isLeaveOneOut = FALSE)
summarizeFilterResults(testMetaObject, getMostRecentFilter(testMetaObject))
}

\keyword{ methods }
\keyword{ classif }
