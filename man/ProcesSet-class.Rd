\name{ProcesSet-class}
\docType{class}
\alias{ProcesSet-class}
\alias{ProcesSet}
\alias{attachData,ProcesSet-method}
\alias{boxplot,ProcesSet-method}
\alias{chipName,ProcesSet-method}
\alias{chipType,ProcesSet-method}
\alias{export,ProcesSet-method}
\alias{getTreeData,ProcesSet-method}
\alias{hist,ProcesSet-method}
\alias{image,ProcesSet-method}
\alias{mboxplot,ProcesSet-method}
\alias{removeData,ProcesSet-method}
\alias{schemeFile,ProcesSet-method}
\alias{schemeFile}
\alias{schemeFile<-,ProcesSet,character-method}
\alias{schemeFile<-}
\alias{schemeSet,ProcesSet-method}
\alias{schemeSet}
\alias{schemeSet<-,ProcesSet,SchemeTreeSet-method}
\alias{schemeSet<-}
\alias{treeData,ProcesSet-method}
\alias{treeData}
\alias{validData,ProcesSet-method}

\title{Class ProcesSet}
\description{This class provides access to class \code{\link{SchemeTreeSet}} for the derived 
classes \code{\link{DataTreeSet}}, \code{\link{ExprTreeSet}} and \code{\link{CallTreeSet}}. 
It extends class \code{\link{TreeSet}}.}
\section{Objects from the Class}{
Usually, no objects are created from it.
}
\section{Slots}{
	 \describe{
    \item{\code{scheme}:}{Object of class \code{"SchemeTreeSet"} providing access to 
     \code{\link{ROOT}} scheme file.}
    \item{\code{data}:}{Object of class \code{"data.frame"}. The data.frame can contain 
     the data stored in \code{\link{ROOT}} data trees.}
    \item{\code{params}:}{Object of class \code{"list"} representing relevant parameters.}
    \item{\code{setname}:}{Object of class \code{"character"} representing the name to the
     \code{\link{ROOT}} file subdirectoy where the \code{\link{ROOT}} trees are stored,
     usually one of \sQuote{DataTreeSet}, \sQuote{PreprocesSet}, \sQuote{CallTreeSet}.}
    \item{\code{settype}:}{Object of class \code{"character"} describing the type of treeset
     stored in \code{setname}, usually one of \sQuote{rawdata}, \sQuote{preprocess}.}
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
Class \code{"\linkS4class{TreeSet}"}, directly.
}
\section{Methods}{
  \describe{
    \item{attachData}{\code{signature(object = "ProcesSet")}: exports data from 
     \code{\link{ROOT}} data file and and saves as data.frame \code{data}.}
    \item{boxplot}{\code{signature(x = "ProcesSet")}: creates a \code{\link{boxplot}} of the
     data from data.frame \code{data}.}
    \item{chipName}{\code{signature(object = "ProcesSet")}: extracts slot \code{chipname} 
     from slot \code{scheme}.}
    \item{chipType}{\code{signature(object = "ProcesSet")}:  extracts slot \code{chiptype} 
     from slot \code{scheme}.}
    \item{export}{\code{signature(object = "ProcesSet")}: exports \code{\link{ROOT}} trees 
     as text file, see \code{\link{export-methods}}.}
    \item{getTreeData}{\code{signature(object = "ProcesSet")}: exports tree data from \code{\link{ROOT}} 
     file \code{rootfile}, and saves as data.frame \code{data}.}
    \item{hist}{\code{signature(x = "ProcesSet")}: creates a plot showing the histograms 
     for data.frame \code{data}.}
    \item{image}{\code{signature(x = "ProcesSet")}: creates an image for each column from
     data.frame \code{data} or \code{bgrd}, respectively.}
    \item{mboxplot}{\code{signature(x = "ProcesSet")}: creates an M-boxplot of the
     data from data.frame \code{data}.}
    \item{removeData}{\code{signature(object = "ProcesSet")}: replaces data.frame \code{data} 
     with an empty data.frame of dim(0,0).}
    \item{schemeFile}{\code{signature(object = "ProcesSet")}: extracts the \code{\link{ROOT}}  
     scheme file from slot \code{scheme}.}
    \item{schemeFile<-}{\code{signature(object = "ProcesSet"), value = "character")}: replaces the \code{\link{ROOT}}  
     scheme file from slot \code{scheme}.}
    \item{schemeSet}{\code{signature(object = "ProcesSet")}: extracts slot \code{scheme}.}
    \item{schemeSet<-}{\code{signature(object = "ProcesSet"), value = "SchemeTreeSet")}: replaces slot \code{scheme} 
     with a different \code{SchemeTreeSet}.}
    \item{treeData}{\code{signature(object = "ProcesSet")}: extracts all columns 
     from data.frame \code{data}.}
    \item{validData}{\code{signature(object = "ProcesSet")}: extracts a subset of columns 
     from data.frame \code{data}.}
	 }
}
\author{Christian Stratowa}
\seealso{
derived classes \code{\linkS4class{DataTreeSet}}, \code{\linkS4class{ExprTreeSet}}, 
\code{\linkS4class{CallTreeSet}}, \code{\linkS4class{QualTreeSet}}.
}
\examples{
showClass("ProcesSet")
}
\keyword{classes}
