\name{gapFilter-methods}
\docType{methods}
\alias{gapFilter-methods}
\alias{gapFilter}
\alias{gapFilter<-}
\title{Gap Filter}
\description{
 This method initializes the Gap Filter. \cr
 The \code{gapFilter} looks for genes that might usefully discriminate between two groups.
 To do this we look for a gap in the ordered expression values. The gap should come in the
 central portion, thus a parameter \code{window} is defined to exclude jumps in the initial
 \code{window} values and the final \code{window} values. \cr
 The Gap Filter flags all rows with: \code{flag = ((gap[i+1] - gap[i])/mean >= cutoff)}

\code{gapFilter(object)} \cr
\code{gapFilter(object, value)<-}
}
\arguments{
  \item{object}{object of class \code{PreFilter}.}
  \item{value}{numeric vector \code{c(cutoff, window, trim, epsilon)}.}
}
\details{
 The method \code{gapFilter} initializes the following parameters:
 \tabular{lll}{
  \tab \code{cutoff}:\tab the cutoff level for the filter. \cr
  \tab \code{window}:\tab trim value for the ordered expression levels (default is \code{window=0.05}). \cr
  \tab \code{trim}:\tab the trim value for trimmed mean (default is \code{trim=0}). \cr
  \tab \code{epsilon}:\tab value to replace mean (default is \code{epsilon=0.01}): \cr
  \tab                \tab \code{epsilon > 0}: replace mean=0 with epsilon. \cr
  \tab                \tab \code{epsilon = 0}: always set mean=1.
 }
 Note, that for \code{epsilon = 0} the filter flags all rows with: (gap[i+1] - gap[i]) >= cutoff
}
\value{
  An initialized \code{PreFilter} object.
}
\author{Christian Stratowa}
\examples{
prefltr <- PreFilter()
gapFilter(prefltr) <- c(0.3,0.05,0.0,0.01)
str(prefltr)
}
\keyword{methods}
