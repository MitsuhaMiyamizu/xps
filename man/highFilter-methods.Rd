\name{highFilter-methods}
\docType{methods}
\alias{highFilter-methods}
\alias{highFilter}
\alias{highFilter<-}
\title{Upper Threshold Filter}
\description{
 This method initializes the Upper Threshold Filter. \cr
 The \code{cutoff} value defines the upper threshold for allowed expression levels. If e.g.
 the number of samples exceeding this cutoff value is greater than \code{parameter} then the
 corresponding dataframe row is flagged, i.e. \code{flag = 0}. \cr
 The Upper Threshold Filter flags all rows with: \code{flag = (sum(expression[i] <= cutoff) >= parameter)}

\emph{Usage}

\code{highFilter(object)} \cr
\code{highFilter(object, value)<-}
}
\arguments{
  \item{object}{object of class \code{PreFilter}.}
  \item{value}{character vector \code{c(cutoff, parameter, condition)}.}
}
\details{
 The method \code{highFilter} initializes the following parameters:
 \tabular{lll}{
  \tab \code{cutoff}:\tab the upper threshold level for the filter. \cr
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
highFilter(prefltr) <- c(14.5,75.0,"percent")
str(prefltr)
}
\keyword{methods}
