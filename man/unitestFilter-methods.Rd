\name{unitestFilter-methods}
\docType{methods}
\alias{unitestFilter-methods}
\alias{unitestFilter}
\alias{unitestFilter<-}
\title{Unitest Filter}
\description{
 This method initializes the Unitest Filter. \cr
 Applying an unitest such as the \code{t.test} to two groups returns the p-value for the test
 and the value of the t-statistic. The Unitest Filter allows to select only rows satisfying
 e.g. a certain p-value as cutoff. \cr
 The Unitest Filter flags all rows with: \code{flag = (variable <= cutoff)}

\emph{Usage}

\code{unitestFilter(object)} \cr
\code{unitestFilter(object, value)<-}
}
\arguments{
  \item{object}{object of class \code{UniFilter}.}
  \item{value}{character vector \code{c(cutoff, variable)}.}
}
\details{
 The method \code{unitestFilter} initializes the following parameters:
 \tabular{lll}{
  \tab \code{cutoff}:\tab the cutoff level for the filter. \cr
  \tab \code{variable}:\tab \code{variable="pval"} (default): p-value. \cr
  \tab                 \tab \code{variable="stat"}: univariate statistic. \cr
  \tab                 \tab \code{variable="padj"}: optional adjusted p-value. \cr
  \tab                 \tab \code{variable="pcha"}: optional p-value obtained by permutations.
 }
}
\value{
  An initialized \code{\linkS4class{UniFilter}} object.
}
\author{Christian Stratowa}
\examples{
unifltr <- UniFilter()
unitestFilter(unifltr) <- c(0.01,"pval")
str(unifltr)
}
\keyword{methods}
