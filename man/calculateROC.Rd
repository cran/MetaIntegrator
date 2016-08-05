\name{calculateROC}
\alias{calculateROC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate ROC Curve Statistics
}
\description{
Calculates receiver operating characteristic curve data, including AUC (using trapezoidal method). Takes only a vector of labels and a vector of predictions. 
}
\usage{
calculateROC(labels, predictions, AUConly = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{labels}{
Vector of labels; must have exactly two unique values (ie, cases and controls).
}
  \item{predictions}{
Vector of predictions (for instance, test scores) to be evaluated for ability to separate the two classes. Must be exactly the same length as labels. 
}
  \item{AUConly}{
Return all ROC values, or just the AUC.
}
}
\details{
The code borrows its core ROC calculations from the ROCR package. AUC is calculated by the trapezoidal method. AUC standard errors are calculated according to Hanley's method.
}
\value{
Assuming AUConly=F, returns a list of values:
\item{roc}{dataframe consisting of two columns, FPR and TPR, meant for plotting}
\item{auc}{area under the curve}
\item{auc.CI}{95\% confidence interval for AUC}
}
\references{
The code borrows its core ROC calculations from the ROCR package.
}
\author{
	Timothy E. Sweeney
}

\seealso{
\code{\link{calculateScore}}, \code{\link{rocPlot}}
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
\keyword{ classify }
