\name{varFilter-methods}
\docType{methods}
\alias{varFilter-methods}
\alias{varFilter}
\alias{varFilter<-}
\title{Variance Filter}
\description{
 This method initializes the Variance Filter. \cr
 The Variance Filter flags all rows with: \code{flag = (var/mean >= cutoff)}

\emph{Usage}

\code{varFilter(object)} \cr
\code{varFilter(object, value)<-}
}
\arguments{
  \item{object}{object of class \code{PreFilter}.}
  \item{value}{numeric vector \code{c(cutoff, trim, epsilon)}.}
}
\details{
 The method \code{varFilter} initializes the following parameters:
 \tabular{lll}{
  \tab \code{cutoff}:\tab the cutoff level for the filter. \cr
  \tab \code{trim}:\tab the trim value for trimmed mean (default is \code{trim=0}). \cr
  \tab \code{epsilon}:\tab value to replace mean (default is \code{epsilon=0.01}): \cr
  \tab                \tab \code{epsilon > 0}: replace mean=0 with epsilon. \cr
  \tab                \tab \code{epsilon = 0}: always set mean=1.
 }
 Note, that for \code{epsilon = 0} the filter flags all rows with: variance >= cutoff
}
\value{
  An initialized \code{\linkS4class{PreFilter}} object.
}
\author{Christian Stratowa}
\examples{
prefltr <- PreFilter()
varFilter(prefltr) <- c(0.6,0.02,0.01)
str(prefltr)
}
\keyword{methods}
