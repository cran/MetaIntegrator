% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculatePower.R
\name{calcMetaPower}
\alias{calcMetaPower}
\title{Calculates the statistical power of a random effects meta-analysis}
\usage{
calcMetaPower(es, avg_n, nStudies, hg, tail=2)
}
\arguments{
\item{es}{effect size you're trying to detect (e.g. 0.6)}

\item{avg_n}{the average sample size of each GROUP in each STUDY (e.g. 10)}

\item{nStudies}{the number of studies you put in the meta-analysis (aka Discovery cohort) (e.g. 5)}

\item{hg}{heterogeneity, (".33" for small, "1" for moderate, & "3" for large) (e.g. 0.33)}

\item{tail}{whether you have a one tail or two tail p-value}
}
\value{
Statistic \code{Power} of the random effects meta-analysis described. 
Most statisticians want a statistical power of at least 0.8, which means that there is an 80% chance
that if there is a true effect, you will detect it.
}
\description{
Calculates the statistical power of a random effects meta-analysis
based on the methods described by Valentine et al. 2010, J of Educational and Behavioral Studies.
}
\details{
Based on the paper by Valentine et al.:
JC Valentine, TD Pigott, and HR Rothstein. 
How Many Studies Do You Need? A Primer on Statistical Power for Meta-Analysis
J of Educational and Behavioral Statistics
April 2010 Vol 35, No 2, pp 215-247

The code itself is adapted from a blog post by 
Dan Quintana, Researcher at Oslo University in Biological Psychiatry
On the website Towards Data Science, July 2017

https://towardsdatascience.com/how-to-calculate-statistical-power-for-your-meta-analysis-e108ee586ae8

\code{avg_n} is the average number people in each group in each study, so if you have
4 studies, and each study compared 10 cases and 10 controls, then \code{avg_n} = 10. 

NOTE: THIS CODE DOES NOT TAKE MULTIPLE HYPOTHESIS TESTING INTO ACCOUNT
      IT ASSUMES P< 0.05

For clarity, avg_n is the average number people in each group in each study, so if you have
4 studies, and each study compared 10 cases and 10 controls, then avg_n = 10.
}
\examples{
# effect size =0.7
# 10 samples on average in each group in each study
# 5 studies included in meta-analysis
# low heterogeneity (0.33)
calcMetaPower(es=0.7, avg_n=10, nStudies=5, hg=0.33)
}
