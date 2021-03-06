% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cleanUpPheno.R
\name{cleanUpPheno}
\alias{cleanUpPheno}
\title{Automatic preprocessing of $pheno dataframe}
\usage{
cleanUpPheno(myDataset)
}
\arguments{
\item{myDataset}{a \code{datasetObject} that contains unprocessed \code{$pheno}}
}
\value{
myDataset a \code{datasetObject} that contains processed \code{$pheno} and original unprocessed \code{$rawPheno}
}
\description{
Takes a Dataset object and:
}
\examples{

\dontrun{
# Download and automatically preprocess pheno
gse53195 = getGEOData("GSE53195")
gse53195 = gse53195$originalData$GSE53195
View(gse53195$pheno) # Original $pheno
gse53195 = cleanUpPheno(gse53195)
View(gse53195$rawPheno) # Original $pheno
View(gse53195$pheno) # Preprocessed $Pheno
}

}
\author{
Erika Bongen
}
