\name{quantileFilter-methods}
\docType{methods}
\alias{quantileFilter-methods}
\alias{quantileFilter}
\alias{quantileFilter<-}
\title{Quantile Filter}
\description{
 This method initializes the Quantile Filter. \cr
 The Quantile Filter flags all rows with: \code{flag = (quantile[high]/quantile[low] >= cutoff)}

\emph{Usage}

\code{quantileFilter(object)} \cr
\code{quantileFilter(object, value)<-}
}
\arguments{
  \item{object}{object of class \code{PreFilter}.}
  \item{value}{numeric vector \code{c(cutoff, loquantile, hiquantile)}.}
}
\details{
 The method \code{quantileFilter} initializes the following parameters:
 \tabular{lll}{
  \tab \code{cutoff}:\tab the cutoff level for the filter. \cr
  \tab \code{loquantile}:\tab value for low  quantile (default is \code{loquantile=0.05}). \cr
  \tab \code{hiquantile}:\tab value for high quantile (default is \code{hiquantile=0.95}).
 }
}
\value{
  An initialized \code{\linkS4class{PreFilter}} object.
}
\author{Christian Stratowa}
\examples{
prefltr <- PreFilter()
quantileFilter(prefltr) <- c(3.0, 0.05, 0.95)
str(prefltr)
}
\keyword{methods}
