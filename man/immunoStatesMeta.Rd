% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/immunoStatesMeta.R
\name{immunoStatesMeta}
\alias{immunoStatesMeta}
\title{immunoStates deconvolution analysis on MetaIntegrator object(s)}
\usage{
immunoStatesMeta(metaObject)
}
\arguments{
\item{metaObject}{a MetaIntegrator formatted Meta object.}
}
\value{
Results from immunoStates stored in $originalData
}
\description{
Run immunoStates and load the results into $originalData for 
running meta-analysis on the cell proportion estimates.
}
\examples{
\dontrun{
# Example won't work on tinyMetaObject because it requires real gene names
# Download the needed datasets for processing. 
sleData <- getGEOData(c("GSE11909","GSE50635", "GSE39088"))

# Run immunoStates
immunoStatesEstimates <- immunoStateMeta(sleData)
}
}
\author{
Francesco Vallania
}
