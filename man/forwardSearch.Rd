\name{forwardSearch}
\alias{forwardSearch}

\title{
Forward Search Function
}
\description{
Forward search is useful for reducing the size of the gene set in your filterObject. In general, forward search identifies a small set of genes with maximum ability to distinguish cases from controls.

forwardSearch is a method of optimizing a given set of significant genes to maximize discriminatory power, as measured by area under the ROC curve (AUC). The function works by taking a given set of genes (presumably a set that has been filtered for statistical significance), and iteratively adding one gene at a time, until the stopping threshold is reached. At each round, the gene whose addition contributes the greatest increase in weighted AUC is added. Weight AUC is defined as the sum of the AUC of each dataset, times the number of samples in that dataset. The stopping threshold is in units of weighted AUC. 
}
\usage{
forwardSearch(metaObject, filterObject, yes.pos = NULL, yes.neg = NULL, forwardThresh = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{metaObject}{
The metaObject from the main metaIntegrator function.
}
  \item{filterObject}{
  An object matching the specifications for Filter
}
  \item{yes.pos}{
Optional- if passed, the forwardSearch will start with the genes in yes.pos and yes.neg (instead of starting from zero genes).
}
  \item{yes.neg}{
Optional- if passed, the forwardSearch will start with the genes in yes.pos and yes.neg (instead of starting from zero genes).
}
  \item{forwardThresh}{
Stopping threshold for the forward search. Default=0.
}
}
\details{
The forwardSearch and backwardSearch functions are designed to assist in selection of gene sets optimized for discriminatory power. The selection of an optimized set is a non-convex problem, and hence both functions will yield gene sets that are only locally optimized (ie, they are not global optima). Both the forwardSearch and backwardSearch functions follow a greedy algorithm, either adding (or removing) genes that contribute the most (or the least) to the overall weighted AUC of the discovery datasets from the metaObject. 

Both search functions allow a user to set a stopping threshold; the fundamental traseoff here will be sparsity of the returned gene set vs. overall discriminatory power. The default threshold is 0, such the functions will return the set of genes at which no gene could be added or removed for the forward or backward functions, respectively, and increase the weighted AUC. 

Note that the weighted AUC returned during the function run is dependent on sample size; this was done (instead of a simple mean) so that the gene set discriminates the MOST SAMPLES, rather than being optimized for any particular dataset. 

}
\value{
A Filter  object which has results from forward search
}
\references{
Sweeney et al., Science Translational Medicine, 2015
}
\author{
Timothy E. Sweeney
}
\seealso{
\code{\link{backwardSearch}}
}
\examples{
#Run forward search to reduce the size of our filter results
forwardRes <- forwardSearch(tinyMetaObject, tinyMetaObject$filterResults[[1]], forwardThresh = 0) 

#See the results
print(forwardRes$posGeneNames)
print(forwardRes$negGeneNames)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ optimize }

