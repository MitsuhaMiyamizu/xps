\name{ratioFilter-methods}
\docType{methods}
\alias{ratioFilter-methods}
\alias{ratioFilter}
\alias{ratioFilter<-}
\title{Ratio Filter}
\description{
 This method initializes the Ratio Filter.
 The ratio is the maximum value divided by minimum value for each row of the expression dataframe. \cr
 The Ratio Filter flags all rows with: \code{flag = (max/min >= cutoff)}

\emph{Usage}

\code{ratioFilter(object)} \cr
\code{ratioFilter(object, value)<-}
}
\arguments{
  \item{object}{object of class \code{PreFilter}.}
  \item{value}{numeric value \code{c(cutoff)}.}
}
\details{
 The method \code{ratioFilter} initializes the following parameters:
 \tabular{lll}{
  \tab \code{cutoff}:\tab the cutoff level for the filter. \cr
 }
}
\value{
  An initialized \code{\linkS4class{PreFilter}} object.
}
\author{Christian Stratowa}
\examples{
prefltr <- PreFilter()
ratioFilter(prefltr) <- c(1.5)
str(prefltr)
}
\keyword{methods}
