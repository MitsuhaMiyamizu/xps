\name{callFilter-methods}
\docType{methods}
\alias{callFilter-methods}
\alias{callFilter}
\alias{callFilter<-}
\title{Detection Call Filter}
\description{
 Detection Call Filter. \cr
 The \code{cutoff} value defines the upper threshold for allowed detection call p-values. If e.g.
 the number of samples exceeding this cutoff value is greater than \code{samples} then the
 corresponding expression dataframe row is flagged, i.e. \code{flag = 0}. \cr
 The Detection Call Filter flags all rows with: \code{flag = (sum(call[i] >= cutoff) >= samples)}

\emph{Usage}

\code{callFilter(object)} \cr
\code{callFilter(object, value)<-}
}
\arguments{
  \item{object}{object of class \code{PreFilter} or \code{UniFilter}.}
  \item{value}{character vector \code{c(cutoff, samples, condition)}.}
}
\details{
 The method \code{callFilter} initializes the following parameters:
 \tabular{lll}{
  \tab \code{cutoff}:\tab the cutoff value for the filter: \cr
  \tab               \tab \code{cutoff = 1.0}: present/absent call is used. \cr
  \tab               \tab \code{cutoff < 1.0}: detection p-value is used as cutoff. \cr
  \tab \code{samples}:\tab this value depends on the \code{condition} used: \cr
  \tab \code{condition}:\tab \code{condition="samples"}: number of samples (default): \cr
  \tab                  \tab \code{condition="percent"}: percent of samples.
 }
}
\value{
  An initialized \code{\linkS4class{PreFilter}} or \code{\linkS4class{UniFilter}} object.
}
\author{Christian Stratowa}
\examples{
## initialize PreFilter
prefltr <- PreFilter()
callFilter(prefltr) <- c(0.02,80.0,"percent")
str(prefltr)

## initialize UniFilter
unifltr <- UniFilter()
callFilter(unifltr) <- c(0.02,80.0,"percent")
str(unifltr)
}
\keyword{methods}
