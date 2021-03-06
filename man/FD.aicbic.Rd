% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FDinterfaces.R
\name{FD.aicbic}
\alias{FD.aicbic}
\title{Information Criterions of a Flexible Dirichlet Model}
\usage{
FD.aicbic(x)
}
\arguments{
\item{x}{an object of class FDfitted, usually the result of \code{\link{FD.estimation}}.}
}
\description{
Akaike Information Criterion (AIC) and Bayesian Information Criterion (BIC) of a fitted Flexible Dirichlet model.
An Information Criterion for one fitted model object for which a log-likelihood value can be obtained is defined as
\eqn{ -2*log-likelihood + k*npar}, where \eqn{ npar} represents the number of parameters in the fitted model, and \eqn{ k = 2} for AIC, or \eqn{ k = log(n)} for BIC (\eqn{ n} being the number of observations).
}
\examples{
data <- FD.generate(n=20,a=c(12,7,15),p=c(0.3,0.4,0.3),t=8)
data
results <- FD.estimation(data, normalize=TRUE,iter.initial.SEM = 5,iter.final.EM = 10)
results
FD.aicbic(results)

}
\seealso{
\code{\link{FD.estimation}}, \code{\link{FD.stddev}}, \code{\link{FD.barycenters}}
}
