% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FDinterfaces.R
\name{FD.clusterdistances}
\alias{FD.clusterdistances}
\title{Flexible Dirichlet Cluster Distances}
\usage{
FD.clusterdistances(x)
}
\arguments{
\item{x}{an object of class FDfitted, usually the result of \code{\link{FD.estimation}}.}
}
\description{
Returns a measure of symmetrized Kullback-Leibler distance between mixture component densities of a fitted Flexible Dirichlet distribution.
}
\examples{
data <- FD.generate(n=20,a=c(12,7,15),p=c(0.3,0.4,0.3),t=8)
data
results <- FD.estimation(data, normalize=TRUE,iter.initial.SEM = 5,iter.final.EM = 10)
results
FD.clusterdistances(results)


}
\references{
{
 Ongaro, A. and Migliorati, S. (2013) A generalization of the Dirichlet distribution. Journal of Multivariate Analysis, \bold{114}, 412--426.\cr
\cr
 Migliorati, S., Ongaro, A. and Monti, G. S. (2016) A structured Dirichlet mixture model for compositional data: inferential and applicative issues. Statistics and Computing, doi:10.1007/s11222-016-9665-y.
}
}
\seealso{
\code{\link{FD.estimation}}, \code{\link{FD.barycenters}}, \code{\link{FD.moments}}
}
