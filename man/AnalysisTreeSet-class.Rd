\name{AnalysisTreeSet-class}
\docType{class}
\alias{AnalysisTreeSet-class}
\alias{AnalysisTreeSet}
\alias{filterTreeset,AnalysisTreeSet-method}
\alias{filterTreeset}
\alias{getTreeData,AnalysisTreeSet-method}
\alias{validData,AnalysisTreeSet-method}
\alias{validFilter,AnalysisTreeSet-method}
\alias{validFilter}
\alias{volcanoplot,AnalysisTreeSet-method}

\title{Class AnalysisTreeSet}
\description{This class provides the link to the \code{\link{ROOT}} analysis file and 
the \code{\link{ROOT}} trees contained therein. It extends class \code{\link{ProcesSet}}.}
\section{Objects from the Class}{
Objects are currently created using function \code{\link{unifilter}}.
}
\section{Slots}{
	 \describe{
    \item{\code{fltrset}:}{Object of class \code{"FilterTreeSet"} providing indirect access to 
     the \code{\linkS4class{ExprTreeSet}} used and the \code{\linkS4class{UniFilter}} settings.}
    \item{\code{scheme}:}{Object of class \code{"SchemeTreeSet"} providing access to 
     \code{\link{ROOT}} scheme file.}
    \item{\code{data}:}{Object of class \code{"data.frame"}. The data.frame contains
     the data of the unitest stored in \code{\link{ROOT}} data trees.}
    \item{\code{params}:}{Object of class \code{"list"} representing relevant parameters.}
    \item{\code{setname}:}{Object of class \code{"character"} representing the name to the
     \code{\link{ROOT}} file subdirectoy where the \code{\link{ROOT}} trees are stored,
     currently \sQuote{UniFilterSet}.}
    \item{\code{settype}:}{Object of class \code{"character"} describing the type of treeset
     stored in \code{setname}, currently \sQuote{unifilter}.}
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
    \item{filterTreeset}{\code{signature(object = "AnalysisTreeSet")}: extracts slot \code{fltrset}.}
    \item{getTreeData}{\code{signature(object = "AnalysisTreeSet")}: exports tree data and returns
     a data.frame.}
    \item{validData}{\code{signature(object = "AnalysisTreeSet")}: extracts data.frame \code{data}.}
    \item{validFilter}{\code{signature(object = "AnalysisTreeSet")}: extracts data.frame \code{data}
     from \code{fltrset}.}
    \item{volcanoplot}{\code{signature(x = "AnalysisTreeSet")}: creates a volcano-plot.}
	 }
}
\author{Christian Stratowa}
\seealso{
related classes \code{\linkS4class{FilterTreeSet}}.
}
\examples{
showClass("AnalysisTreeSet")
}
\keyword{classes}
