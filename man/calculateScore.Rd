\name{calculateScore}
\alias{calculateScore}

\title{
Calculate a signature Z-score for a set of genes in a single dataset
}
\description{
Given a gene set of interest, it is often desirable to summarize the expression of that gene set using a single integrated score. The \code{calculateScore} method calculates the geometric mean of the expression level of all positive genes, minus the geometric mean of the expression level of all negative genes. The resulting scores are then standardized within the given dataset, such that the output 'Z-score' has mean=0 and std. dev=1. Such a Z-score can then be used for classification, etc. 
}
\usage{
calculateScore(filterObject, datasetObject, suppressMessages=FALSE)
}
\arguments{
  \item{filterObject}{
a MetaFilter object generated with \code{filterGenes()} containing the signature genes that will be used for Z-score calculation. 
}
  \item{datasetObject}{
A Dataset object for which the signature score (Z-score) will be calculated. This vector would typically be added as \code{$score} column in \code{datasetObject$pheno}.
}
  \item{suppressMessages}{
	  Boolean value (TRUE/FALSE) about whether to display verbose output. Default: FALSE.}
}
\details{
The Z-score is based off of the geometric mean of expression. As such, negative expression values are not allowed. A dataset is thus always scaled by its minimum value + 1, such that the lowest value = 1. Any individual NANs or NAs are also set to 1. If a dataset does not have any information on a given gene, the entire gene is simply left out of the score. When run, the function will print to command line the number of genes used, and the number passed in. 

Although mostly used internally, the function has been exported in case users want to compare multiple classes, etc., using the same Z-score as is used for producing two-class comparisons. 
}
\value{
A vector of Z-scores, of length \code{ncols(datasetObject$expr)} (and in the same order).
}
\author{
	Timothy E. Sweeney, Winston A. Haynes
}
\seealso{

\code{\link{filterGenes}}
}
\examples{
calculateScore(tinyMetaObject$filterResults[[1]], tinyMetaObject$originalData[[1]]) 
}
