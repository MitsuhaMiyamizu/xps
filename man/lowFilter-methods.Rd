\name{lowFilter-methods}
\docType{methods}
\alias{lowFilter-methods}
\alias{lowFilter}
\alias{lowFilter<-}
\title{Lower Threshold Filter}
\description{
 This method initializes the Lower Threshold Filter.
 The \code{cutoff} value defines the lower threshold for allowed expression levels. If e.g.
 the number of samples lower than this cutoff value is greater than \code{parameter} then the
 corresponding dataframe row is flagged, i.e. \code{flag = 0}. \cr
 The Lower Threshold Filter flags all rows with: \code{flag = (sum(expression[i] >= cutoff) >= parameter)}

\emph{Usage}

\code{lowFilter(object)} \cr
\code{lowFilter(object, value)<-}
}
\arguments{
  \item{object}{object of class \code{PreFilter}.}
  \item{value}{character vector \code{c(cutoff, parameter, condition)}.}
}
\details{
 The method \code{lowFilter} initializes the following parameters:
 \tabular{lll}{
  \tab \code{cutoff}:\tab the lower threshold level for the filter. \cr
  \tab \code{parameter}:\tab this value depends on the \code{condition} used: \cr
  \tab \code{condition}:\tab \code{condition="samples"}:    number of samples (default): \cr
  \tab                  \tab \code{condition="percent"}:    percent of samples. \cr
  \tab                  \tab \code{condition="mean"}:       mean value of samples. \cr
  \tab                  \tab \code{condition="percentile"}: percentile of samples.
 }
}
\value{
  An initialized \code{\linkS4class{PreFilter}} object.
}
\author{Christian Stratowa}
\examples{
prefltr <- PreFilter()
lowFilter(prefltr) <- c(4.0,3,"samples")
str(prefltr)
}
\keyword{methods}

