% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FDinterfaces.R
\name{FD.subcomposition}
\alias{FD.subcomposition}
\title{Subcomposition}
\usage{
FD.subcomposition(data, columns)
}
\arguments{
\item{data}{a matrix or a dataframe containing only variables in the model.}

\item{columns}{numeric vector containing the position of the columns to keep in the new composition.}
}
\description{
Given a matrix or a numeric dataframe, this function returns a subcomposition made up of the specified columns.
}
\details{
Values must be positive. In case one row-entry (or more) is NA, the whole row will be returned as NA.
}
\examples{
data(oliveoil)
dataoil <- oliveoil
head(dataoil)
data <- FD.normalization(dataoil[,3:10])
head(data)
data.sub <- FD.subcomposition(data,c(1,3,4,5))
head(data.sub)
data.amalg <- FD.amalgamation(data,c(2,6,7,8),name='others')
head(data.amalg)

}
\seealso{
\code{\link{FD.generate}}, \code{\link{FD.amalgamation}}, \code{\link{FD.normalization}}
}
