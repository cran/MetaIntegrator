% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ROC_and_search_functions.R
\name{calculateROC}
\alias{calculateROC}
\title{Calculate ROC Curve Statistics}
\usage{
calculateROC(labels, predictions, AUConly = FALSE)
}
\arguments{
\item{labels}{Vector of labels; must have exactly two unique values (ie, cases and controls).}

\item{predictions}{Vector of predictions (for instance, test scores) to be evaluated for ability to separate the two classes. Must be exactly the same length as labels.}

\item{AUConly}{Return all ROC values, or just the AUC.}
}
\value{
Assuming AUConly=F, returns a list of values:
 \item{roc}{dataframe consisting of two columns, FPR and TPR, meant for plotting}
 \item{auc}{area under the curve}
 \item{auc.CI}{95\% confidence interval for AUC}
}
\description{
Calculates receiver operating characteristic curve data, including AUC (using trapezoidal method). Takes only a vector of labels and a vector of predictions.
}
\details{
The code borrows its core ROC calculations from the ROCR package. AUC is calculated by the trapezoidal method. AUC standard errors are calculated according to Hanley's method.
}
\examples{
# expect an AUC near 0.5 with random test
labels <- c(rep(0, 500), rep(1, 500))
scores <- runif(1000)
calculateROC(labels, scores)
#With the real data, AUC should be around 0.85606
scoreResults <- calculateScore(tinyMetaObject$filterResults[[1]], tinyMetaObject$originalData[[1]]) 
rocRes <- calculateROC(predictions=scoreResults, labels=tinyMetaObject$originalData[[1]]$class)
print(rocRes$auc[[1]])
}
\references{
The code borrows its core ROC calculations from the ROCR package.
}
\seealso{
\code{\link{calculateScore}}, \code{\link{rocPlot}}
}
\author{
Timothy E. Sweeney
}
\keyword{classify}
