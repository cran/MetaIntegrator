\name{getSampleLevelGeneData}
\alias{getSampleLevelGeneData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract gene-level data from a given data object
}
\description{
Given a \code{datsetObject}, and a set of target genes, this function will summarize probe-level data to gene-level data for the target genes. Returns a data frame with only the genes of interest, for each sample in the dataset. 
}
\usage{
getSampleLevelGeneData(datasetObject, geneNames)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datasetObject}{
a Dataset object that is used to extract sample level data (At least, must have a \code{$expr} of probe-level data, and \code{$keys} of probe:gene mappings).
}
  \item{geneNames}{
A vector of geneNames
}
}
\details{
Summarizes probe-level data to gene-level data, using the mean of the probes, according to the probe:gene mapping in the \code{$keys} item in the dataset object. This is done only for the genes in the filter object.
}
\value{
Returns a data frame with expression levels of only the genes of interest, for each sample in the dataset. 

Mostly used internally, but has been exposed to the user to allow advanced functionality on external datasets if desired. 
}
\author{
	Timothy E. Sweeney
}

\examples{
sampleResults <- getSampleLevelGeneData(datasetObject=tinyMetaObject$originalData[[1]], 
					geneNames=c(tinyMetaObject$filterResults[[1]]$posGeneNames, 
					tinyMetaObject$filterResults[[1]]$negGeneNames))
}
