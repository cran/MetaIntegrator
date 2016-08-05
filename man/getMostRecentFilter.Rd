\name{getMostRecentFilter}
\alias{getMostRecentFilter}
\title{
Get name of most recent filter 
}
\description{
Given a \code{metaObject} this function will look through \code{$filterResults} for the most recent filter used and return the filter name.
}
\usage{
getMostRecentFilter(metaObject)
}
\arguments{
  \item{metaObject}{A meta object}
}
\value{
  \item{FilterLabel}{Name of the most recent filter}
}
\author{
	Francesco Vallania
}
\examples{
getMostRecentFilter(tinyMetaObject)
}
\keyword{ attribute }
\keyword{ utilities }
