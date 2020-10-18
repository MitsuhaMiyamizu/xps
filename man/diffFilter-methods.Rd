\name{diffFilter-methods}
\docType{methods}
\alias{diffFilter-methods}
\alias{diffFilter}
\alias{diffFilter<-}
\title{Difference Filter}
\description{
 This method initializes the Difference Filter. \cr
 The difference is the maximum value minus minimum value for each row of the expression dataframe
 divided by the mean value of each row. \cr
 The Difference Filter flags all rows with: \code{flag = ((max - min)/mean >= cutoff)}

\emph{Usage}

\code{diffFilter(object)} \cr
\code{diffFilter(object, value)<-}
}
\arguments{
  \item{object}{object of class \code{PreFilter}.}
  \item{value}{numeric vector \code{c(cutoff, trim, epsilon)}.}
}
\details{
 The method \code{diffFilter} initializes the following parameters:
 \tabular{lll}{
  \tab \code{cutoff}:\tab the cutoff level for the filter. \cr
  \tab \code{trim}:\tab the trim value for trimmed mean (default is \code{trim=0}). \cr
  \tab \code{epsilon}:\tab value to replace mean (default is \code{epsilon=0.01}): \cr
  \tab                \tab \code{epsilon > 0}: replace mean=0 with epsilon. \cr
  \tab                \tab \code{epsilon = 0}: always set mean=1.
 }
 Note, that for \code{epsilon = 0} the filter flags all rows with: (max - min) >= cutoff
}
\value{
  An initialized \code{\linkS4class{PreFilter}} object.
}
\author{Christian Stratowa}
\examples{
prefltr <- PreFilter()
diffFilter(prefltr) <- c(2.2,0.0,0.01)
str(prefltr)
}
\keyword{methods}
