\name{UniFilter-class}
\docType{class}
\alias{UniFilter-class}
\alias{callFilter,UniFilter-method}
\alias{callFilter<-,UniFilter,character-method}
\alias{fcFilter,UniFilter-method}
\alias{fcFilter<-,UniFilter,character-method}
\alias{uniTest,UniFilter-method}
\alias{uniTest<-,UniFilter,character-method}
\alias{unitestFilter,UniFilter-method}
\alias{unitestFilter<-,UniFilter,character-method}

\title{Class UniFilter}
\description{Class UniFilter allows to apply different unitest filters to class \code{\link{ExprTreeSet}},
 i.e. to the expression level data.frame \code{data}.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("UniFilter", ...)}. \cr
Alternatively, the contructor \code{\link{UniFilter}} can be used.
}
\section{Slots}{
	 \describe{
    \item{\code{foldchange}:}{Object of class \code{"list"} describing parameters for \code{fcFilter}.}
    \item{\code{prescall}:}{Object of class \code{"list"} describing parameters for \code{callFilter}.}
    \item{\code{unifilter}:}{Object of class \code{"list"} describing parameters for \code{unitestFilter}.}
    \item{\code{unitest}:}{Object of class \code{"list"} describing parameters for \code{uniTest}.}
    \item{\code{numfilters}:}{Object of class \code{"numeric"} giving the number of filters applied.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{Filter}"}, directly.
}
\section{Methods}{
  \describe{
    \item{callFilter}{\code{signature(object = "UniFilter")}: extracts slot \code{prescall}.}
    \item{callFilter<-}{\code{signature(object = "UniFilter", value = "character")}: replaces 
     slot \code{prescall} with character vector c(cutoff, samples, condition).}
    \item{fcFilter}{\code{signature(object = "UniFilter")}: extracts slot \code{foldchange}.}
    \item{fcFilter<-}{\code{signature(object = "UniFilter", value = "numeric")}: replaces 
     slot \code{foldchange} with numeric vector c(cutoff, direction).}
    \item{uniTest}{\code{signature(object = "UniFilter")}: extracts slot \code{unitest}.}
    \item{uniTest<-}{\code{signature(object = "UniFilter", value = "character")}: replaces 
     slot \code{unitest} with character vector c(type, alternative, correction, numperm, mu, paired, conflevel, varequ).}
    \item{unitestFilter}{\code{signature(object = "UniFilter")}: extracts slot \code{unifilter}.}
    \item{unitestFilter<-}{\code{signature(object = "UniFilter", value = "character")}: replaces 
     slot \code{unifilter} with character vector c(cutoff, variable).}
	 }
}
\author{Christian Stratowa}
\seealso{
related classes \code{\linkS4class{Filter}}, \code{\linkS4class{PreFilter}}.
}
\examples{
unifltr <- new("UniFilter", unitest=list("t.test"))
fcFilter(unifltr) <- c(1.5,"both")
unitestFilter(unifltr) <- c(0.01,"pval")
str(unifltr)
}
\keyword{classes}
