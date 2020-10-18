\name{Filter-class}
\docType{class}
\alias{Filter-class}
\alias{numberFilters,Filter-method}
\alias{numberFilters}

\title{Base Class Filter}
\description{Base class for classes \code{\linkS4class{PreFilter}} and \code{\linkS4class{UniFilter}}.}
\section{Slots}{
	\describe{
    \item{\code{numfilters}:}{Object of class \code{"numeric"} giving the number of filters applied.}
  }
}
\section{Methods}{
\describe{
    \item{numberFilters}{\code{signature(object = "Filter")}: number of filters applied.}
}
}
\author{Christian Stratowa}
\seealso{
related classes \code{\linkS4class{PreFilter}}, \code{\linkS4class{UniFilter}}.
}
\examples{
showClass("Filter")
}
\keyword{classes}
