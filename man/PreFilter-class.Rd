\name{PreFilter-class}
\docType{class}
\alias{PreFilter-class}
\alias{callFilter,PreFilter-method}
\alias{callFilter<-,PreFilter,character-method}
\alias{cvFilter,PreFilter-method}
\alias{cvFilter<-,PreFilter,numeric-method}
\alias{diffFilter,PreFilter-method}
\alias{diffFilter<-,PreFilter,numeric-method}
\alias{gapFilter,PreFilter-method}
\alias{gapFilter<-,PreFilter,numeric-method}
\alias{highFilter,PreFilter-method}
\alias{highFilter<-,PreFilter,character-method}
\alias{lowFilter,PreFilter-method}
\alias{lowFilter<-,PreFilter,character-method}
\alias{madFilter,PreFilter-method}
\alias{madFilter<-,PreFilter,numeric-method}
\alias{quantileFilter,PreFilter-method}
\alias{quantileFilter<-,PreFilter,numeric-method}
\alias{ratioFilter,PreFilter-method}
\alias{ratioFilter<-,PreFilter,numeric-method}
\alias{varFilter,PreFilter-method}
\alias{varFilter<-,PreFilter,numeric-method}

\title{Class PreFilter}
\description{Class PreFilter allows to apply different filters to class \code{\link{ExprTreeSet}},
 i.e. to the expression level data.frame \code{data}.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("PreFilter", ...)}. \cr
Alternatively, the contructor \code{\link{PreFilter}} can be used.
}
\section{Slots}{
	 \describe{
    \item{\code{mad}:}{Object of class \code{"list"} describing parameters for \code{madFilter}.}
    \item{\code{cv}:}{Object of class \code{"list"} describing parameters for \code{cvFilter}.}
    \item{\code{variance}:}{Object of class \code{"list"} describing parameters for \code{varFilter}.}
    \item{\code{difference}:}{Object of class \code{"list"} describing parameters for \code{diffFilter}.}
    \item{\code{ratio}:}{Object of class \code{"list"} describing parameters for \code{ratioFilter}.}
    \item{\code{gap}:}{Object of class \code{"list"} describing parameters for \code{gapFilter}.}
    \item{\code{hithreshold}:}{Object of class \code{"list"} describing parameters for \code{highFilter}.}
    \item{\code{lothreshold}:}{Object of class \code{"list"} describing parameters for \code{lowFilter}.}
    \item{\code{quantile}:}{Object of class \code{"list"} describing parameters for \code{quantileFilter}.}
    \item{\code{prescall}:}{Object of class \code{"list"} describing parameters for \code{callFilter}.}
    \item{\code{numfilters}:}{Object of class \code{"numeric"} giving the number of filters applied.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{Filter}"}, directly.
}
\section{Methods}{
  \describe{
    \item{callFilter}{\code{signature(object = "PreFilter")}: extracts slot \code{prescall}.}
    \item{callFilter<-}{\code{signature(object = "PreFilter", value = "character")}: replaces 
     slot \code{prescall} with character vector c(cutoff, samples, condition).}
    \item{cvFilter}{\code{signature(object = "PreFilter")}: extracts slot \code{cv}.}
    \item{cvFilter<-}{\code{signature(object = "PreFilter", value = "numeric")}: replaces 
     slot \code{cv} with numeric vector c(cutoff, trim, epsilon).}
    \item{diffFilter}{\code{signature(object = "PreFilter")}: extracts slot \code{difference}.}
    \item{diffFilter<-}{\code{signature(object = "PreFilter", value = "numeric")}: replaces 
     slot \code{difference} with numeric vector c(cutoff, trim, epsilon).}
    \item{gapFilter}{\code{signature(object = "PreFilter")}: extracts slot \code{gap}.}
    \item{gapFilter<-}{\code{signature(object = "PreFilter", value = "numeric")}: replaces 
     slot \code{gap} with numeric vector c(cutoff, window, trim, epsilon).}
    \item{highFilter}{\code{signature(object = "PreFilter")}: extracts slot \code{hithreshold}.}
    \item{highFilter<-}{\code{signature(object = "PreFilter", value = "character")}: replaces 
     slot \code{hithreshold} with character vector c(cutoff, parameter, condition).}
    \item{lowFilter}{\code{signature(object = "PreFilter")}: extracts slot \code{lothreshold}.}
    \item{lowFilter<-}{\code{signature(object = "PreFilter", value = "character")}: replaces 
     slot \code{lothreshold} with character vector c(cutoff, parameter, condition).}
    \item{madFilter}{\code{signature(object = "PreFilter")}: extracts slot \code{mad}.}
    \item{madFilter<-}{\code{signature(object = "PreFilter", value = "numeric")}: replaces 
     slot \code{mad} with numeric vector c(cutoff, epsilon).}
    \item{quantileFilter}{\code{signature(object = "PreFilter")}: extracts slot \code{quantile}.}
    \item{quantileFilter<-}{\code{signature(object = "PreFilter", value = "numeric")}: replaces 
     slot \code{quantile} with numeric vector c(cutoff, loquantile, hiquantile).}
    \item{ratioFilter}{\code{signature(object = "PreFilter")}: extracts slot \code{ratio}.}
    \item{ratioFilter<-}{\code{signature(object = "PreFilter", value = "numeric")}: replaces 
     slot \code{ratio} with numeric vector c(cutoff).}
    \item{varFilter}{\code{signature(object = "PreFilter")}: extracts slot \code{variance}.}
    \item{varFilter<-}{\code{signature(object = "PreFilter", value = "numeric")}: replaces 
     slot \code{variance} with numeric vector c(cutoff, trim, epsilon).}
	 }
}
\author{Christian Stratowa}
\seealso{
related classes \code{\linkS4class{Filter}}, \code{\linkS4class{UniFilter}}.
}
\examples{
## for demonstration purposes only:  initialize all pre-filters
prefltr <- new("PreFilter")
madFilter(prefltr) <- c(0.5,0.01)
cvFilter(prefltr) <- c(0.3,0.0,0.01)
varFilter(prefltr) <- c(0.6,0.02,0.01)
diffFilter(prefltr) <- c(2.2,0.0,0.01)
ratioFilter(prefltr) <- c(1.5)
gapFilter(prefltr) <- c(0.3,0.05,0.0,0.01)
lowFilter(prefltr) <- c(4.0,3,"samples")
highFilter(prefltr) <- c(14.5,75.0,"percent")
quantileFilter(prefltr) <- c(3.0, 0.05, 0.95)
callFilter(prefltr) <- c(0.02,80.0,"percent")
str(prefltr)
}
\keyword{classes}
