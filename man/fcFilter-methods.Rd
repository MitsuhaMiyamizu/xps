\name{fcFilter-methods}
\docType{methods}
\alias{fcFilter-methods}
\alias{fcFilter}
\alias{fcFilter<-}
\title{Fold-Change Filter}
\description{
 This method initializes the Fold-Change Filter. \cr
 The fold-change is determined by the mean value of group 2 divided by the mean value of group 1. \cr
 The Fold-Change Filter flags all rows with: \code{flag = (fc >= cutoff)}

\emph{Usage}

\code{fcFilter(object)} \cr
\code{fcFilter(object, value)<-}
}
\arguments{
  \item{object}{object of class \code{UniFilter}.}
  \item{value}{numeric vector \code{c(cutoff, direction)}}
}
\details{
 The method \code{fcFilter} initializes the following parameters:
 \tabular{lll}{
  \tab \code{cutoff}:\tab the cutoff level for the filter. \cr
  \tab \code{direction}:\tab \code{direction="both"} (default): select up and downregulated genes. \cr
  \tab                  \tab \code{direction="up"}: select upregulated genes only. \cr
  \tab                  \tab \code{direction="down"}: select downregulated genes only.
 }
}
\value{
  An initialized \code{\linkS4class{UniFilter}} object.
}
\author{Christian Stratowa}
\examples{
unifltr <- UniFilter()
fcFilter(unifltr) <- c(1.5,"both")
str(unifltr)
}
\keyword{methods}
