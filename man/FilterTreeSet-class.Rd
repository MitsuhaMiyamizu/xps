\name{FilterTreeSet-class}
\docType{class}
\alias{FilterTreeSet-class}
\alias{FilterTreeSet}
\alias{callTreeset,FilterTreeSet-method}
\alias{callTreeset}
\alias{exprTreeset,FilterTreeSet-method}
\alias{exprTreeset}
\alias{getTreeData,FilterTreeSet-method}
\alias{validData,FilterTreeSet-method}

\title{Class FilterTreeSet}
\description{This class provides the link to the \code{\link{ROOT}} filter file and 
the \code{\link{ROOT}} trees contained therein. It extends class \code{\link{ProcesSet}}.}
\section{Objects from the Class}{
Objects are currently created using function \code{\link{prefilter}}.
}
\section{Slots}{
	 \describe{
    \item{\code{filter}:}{Object of class \code{"Filter"} currently providing access to 
     the \code{\linkS4class{PreFilter}} settings.}
    \item{\code{exprset}:}{Object of class \code{"ExprTreeSet"} providing direct access to 
     the \code{\linkS4class{ExprTreeSet}} used for filtering.}
    \item{\code{callset}:}{Object of class \code{"CallTreeSet"} providing direct access to 
     the optional \code{\linkS4class{CallTreeSet}} used for filtering.}
    \item{\code{scheme}:}{Object of class \code{"SchemeTreeSet"} providing access to 
     \code{\link{ROOT}} scheme file.}
    \item{\code{data}:}{Object of class \code{"data.frame"}. The data.frame contains
     the data of the filter stored in \code{\link{ROOT}} filter trees.}
    \item{\code{params}:}{Object of class \code{"list"} representing relevant parameters.}
    \item{\code{setname}:}{Object of class \code{"character"} representing the name to the
     \code{\link{ROOT}} file subdirectoy where the \code{\link{ROOT}} trees are stored,
     currently \sQuote{PreFilterSet}.}
    \item{\code{settype}:}{Object of class \code{"character"} describing the type of treeset
     stored in \code{setname}, currently \sQuote{prefilter}.}
    \item{\code{rootfile}:}{Object of class \code{"character"} representing the name 
     of the \code{\link{ROOT}} file, including full path.}
    \item{\code{filedir}:}{Object of class \code{"character"} describing the full path to
     the system directory where \code{rootfile} is stored.}
    \item{\code{numtrees}:}{Object of class \code{"numeric"} representing the number of 
     \code{\link{ROOT}} trees stored in subdirectoy \code{setname}.}
    \item{\code{treenames}:}{Object of class \code{"list"} representing the names of the 
     \code{\link{ROOT}} trees stored in subdirectoy \code{setname}.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{ProcesSet}"}, directly.
Class \code{"\linkS4class{TreeSet}"}, by class "ProcesSet", distance 2.
}
\section{Methods}{
  \describe{
    \item{callTreeset}{\code{signature(object = "FilterTreeSet")}: extracts slot \code{callset}.}
    \item{exprTreeset}{\code{signature(object = "FilterTreeSet")}: extracts slot \code{exprset}.}
    \item{getTreeData}{\code{signature(object = "FilterTreeSet")}: exports tree data and returns
     a data.frame.}
    \item{validData}{\code{signature(object = "FilterTreeSet")}: extracts data.frame \code{data}.}
	 }
}
\author{Christian Stratowa}
\seealso{
related classes \code{\linkS4class{AnalysisTreeSet}}.
}
\examples{
showClass("FilterTreeSet")
}
\keyword{classes}
