\name{cvFilter-methods}
\docType{methods}
\alias{cvFilter-methods}
\alias{cvFilter}
\alias{cvFilter<-}
\title{Coefficient of Variation Filter}
\description{
 This method initializes the Coefficient of Variation Filter. \cr
 The coefficient of variation is the standard deviation divided by the absolute value of the mean. \cr
 The CV Filter flags all rows with: \code{flag = (cv >= cutoff)}

\emph{Usage}

\tabular{ll}{
 \tab \code{cvFilter(object)} \cr
 \tab \code{cvFilter(object, value)<-}
}
}
\arguments{
  \item{object}{object of class \code{PreFilter}.}
  \item{value}{numeric vector \code{c(cutoff, trim, epsilon)}.}
}
\details{
 The method \code{cvFilter} initializes the following parameters:
 \tabular{lll}{
  \tab \code{cutoff}:\tab the cutoff level for the filter. \cr
  \tab \code{trim}:\tab the trim value for trimmed mean (default is \code{trim=0}). \cr
  \tab \code{epsilon}:\tab value to replace mean (default is \code{epsilon=0.01}): \cr
  \tab                \tab \code{epsilon > 0}: replace mean=0 with epsilon. \cr
  \tab                \tab \code{epsilon = 0}: always set mean=1.
 }
 Note, that for \code{epsilon = 0} the filter flags all rows with: stdev >= cutoff
}
\value{
  An initialized \code{\linkS4class{PreFilter}} object.
}
\author{Christian Stratowa}
\examples{
prefltr <- PreFilter()
cvFilter(prefltr) <- c(0.3,0.0,0.01)
str(prefltr)
}
\keyword{methods}
