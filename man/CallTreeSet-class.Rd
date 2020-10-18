\name{CallTreeSet-class}
\docType{class}
\alias{CallTreeSet-class}
\alias{CallTreeSet}
\alias{attachCall,CallTreeSet-method}
\alias{attachPVal,CallTreeSet-method}
\alias{callplot,CallTreeSet-method}
\alias{presCall,CallTreeSet-method}
\alias{presCall<-,CallTreeSet,data.frame-method}
\alias{pvalData,CallTreeSet-method}
\alias{pvalData<-,CallTreeSet,data.frame-method}
\alias{removeCall,CallTreeSet-method}
\alias{removePVal,CallTreeSet-method}
\alias{validCall,CallTreeSet-method}
\alias{validPVal,CallTreeSet-method}

\title{Class CallTreeSet}
\description{This class provides the link to the \code{\link{ROOT}} call file and 
the \code{\link{ROOT}} trees contained therein. It extends class \code{\link{ProcesSet}}.}
\section{Objects from the Class}{
Objects are created using functions \code{\link{mas5.call}} or \code{\link{dabg.call}}, respectively.
}
\section{Slots}{
	 \describe{
    \item{\code{calltype}:}{Object of class \code{"character"} representing the call type, 
     i.e. \sQuote{mas5} or \sQuote{dabg}.}
    \item{\code{detcall}:}{Object of class \code{"data.frame"}. The data.frame can contain 
     the detection calls stored in \code{\link{ROOT}} call trees.}
    \item{\code{scheme}:}{Object of class \code{"SchemeTreeSet"} providing access to 
     \code{\link{ROOT}} scheme file.}
    \item{\code{data}:}{Object of class \code{"data.frame"}. The data.frame can contain 
     the data (i.e. p-values) stored in \code{\link{ROOT}} call trees.}
    \item{\code{params}:}{Object of class \code{"list"} representing relevant parameters.}
    \item{\code{setname}:}{Object of class \code{"character"} representing the name to the
     \code{\link{ROOT}} file subdirectoy where the \code{\link{ROOT}} call trees are stored,
     usually \sQuote{CallTreeSet}.}
    \item{\code{settype}:}{Object of class \code{"character"} describing the type of treeset
     stored in \code{setname}, usually \sQuote{preprocess}.}
    \item{\code{rootfile}:}{Object of class \code{"character"} representing the name 
     of the \code{\link{ROOT}} call file, including full path.}
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
    \item{attachCall}{\code{signature(object = "CallTreeSet")}: exports detection call data from 
     \code{\link{ROOT}} call file and and saves as data.frame \code{detcall}.}
    \item{attachPVal}{\code{signature(object = "CallTreeSet")}: exports call p-values from 
     \code{\link{ROOT}} call file and and saves as data.frame \code{data}.}
    \item{callplot}{\code{signature(x = "CallTreeSet")}: creates a barplot of percent present
     and absent calls.}
    \item{presCall}{\code{signature(object = "CallTreeSet")}: extracts the detection call data.frame.}
    \item{presCall<-}{\code{signature(object = "CallTreeSet", value = "data.frame")}: replaces the detection call data.frame.}
    \item{pvalData}{\code{signature(object = "CallTreeSet")}: extracts the detection p-value data.frame.}
    \item{pvalData<-}{\code{signature(object = "CallTreeSet", value = "data.frame")}: replaces the detection p-value data.frame.}
    \item{removeCall}{\code{signature(object = "CallTreeSet")}: replaces data.frame \code{detcall} 
     with an empty data.frame of dim(0,0).}
    \item{removePVal}{\code{signature(object = "CallTreeSet")}: replaces data.frame \code{data} 
     with an empty data.frame of dim(0,0).}
    \item{validCall}{\code{signature(object = "CallTreeSet")}: extracts a subset of columns 
     from data.frame \code{detcall}.}
    \item{validPVal}{\code{signature(object = "CallTreeSet")}: extracts a subset of columns 
     from data.frame \code{data}.}
	 }
}
\author{Christian Stratowa}
\seealso{
related classes \code{\linkS4class{DataTreeSet}}, \code{\linkS4class{ExprTreeSet}}.
}
\examples{
showClass("CallTreeSet")
}
\keyword{classes}
