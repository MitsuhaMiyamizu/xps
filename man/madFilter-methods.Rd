\name{madFilter-methods}
\docType{methods}
\alias{madFilter-methods}
\alias{madFilter}
\alias{madFilter<-}
\title{Median Absolute Deviation Filter}
\description{
 This method initializes the Median Absolute Deviation Filter. \cr
 The MAD Filter flags all rows with: \code{flag = (mad >= cutoff)}

\emph{Usage}

\code{madFilter(object)} \cr
\code{madFilter(object, value)<-}
}
\arguments{
  \item{object}{object of class \code{PreFilter}.}
  \item{value}{numeric vector \code{c(cutoff, epsilon)}.}
}
\details{
 The method \code{madFilter} initializes the following parameters:
 \tabular{lll}{
  \tab \code{cutoff}:\tab the cutoff level for the filter. \cr
  \tab \code{epsilon}:\tab value to replace mean (default is \code{epsilon=0.01}).
 }
 Note, that \code{epsilon} has no effect on \code{mad}.
}
\value{
  An initialized \code{\linkS4class{PreFilter}} object.
}
\author{Christian Stratowa}
\examples{
prefltr <- PreFilter()
madFilter(prefltr) <- c(0.5,0.01)
str(prefltr)
}
\keyword{methods}
