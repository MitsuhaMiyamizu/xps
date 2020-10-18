\name{TreeSet-class}
\docType{class}
\alias{TreeSet-class}
\alias{TreeSet}
\alias{export,TreeSet-method}
\alias{fileDir,TreeSet-method}
\alias{fileDir}
\alias{fileDir<-,TreeSet,character-method}
\alias{fileDir<-}
\alias{rootFile,TreeSet-method}
\alias{rootFile}
\alias{rootFile<-,TreeSet,character-method}
\alias{rootFile<-}
\alias{root.browser,TreeSet-method}
\alias{setName,TreeSet-method}
\alias{setName}
\alias{setName<-,TreeSet,character-method}
\alias{setName<-}
\alias{setType,TreeSet-method}
\alias{setType}
\alias{setType<-,TreeSet,character-method}
\alias{setType<-}
\alias{treeInfo,TreeSet-method}
\alias{treeNames,TreeSet-method}
\alias{treeNames}

\title{Class TreeSet}
\description{
This is the virtual base class for all other classes providing the link to 
a \code{\link{ROOT}} file and the \code{\link{ROOT}} trees contained therein.
}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
	 \describe{
    \item{\code{setname}:}{Object of class \code{"character"} representing the name to the
     \code{\link{ROOT}} file subdirectory where the \code{\link{ROOT}} trees are stored,
     usually one of \sQuote{DataTreeSet}, \sQuote{PreprocesSet}, \sQuote{CallTreeSet}.}
    \item{\code{settype}:}{Object of class \code{"character"} describing the type of treeset
     stored in \code{setname}, usually one of \sQuote{scheme}, \sQuote{rawdata}, \sQuote{preprocess}.}
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
\section{Methods}{
  \describe{
    \item{export}{\code{signature(object = "TreeSet")}: exports \code{\link{ROOT}} trees 
     as text file, see \code{\link{export-methods}}.}
    \item{fileDir}{\code{signature(object = "TreeSet")}: extracts slot \code{filedir}.}
    \item{fileDir<-}{\code{signature(object = "TreeSet", value = "character")}: replaces slot \code{filedir}.}
    \item{root.browser}{\code{signature(object = "TreeSet")}: opens the \code{\link{ROOT}} file browser.}
    \item{rootFile}{\code{signature(object = "TreeSet")}: extracts slot \code{rootfile}.}
    \item{rootFile<-}{\code{signature(object = "TreeSet", value = "character")}: replaces slot \code{rootfile}.}
    \item{setName}{\code{signature(object = "TreeSet")}: extracts slot \code{setname}.}
    \item{setName<-}{\code{signature(object = "TreeSet", value = "character")}: replaces slot \code{setname}.}
    \item{setType}{\code{signature(object = "TreeSet")}: extracts slot \code{settype}.}
    \item{setType<-}{\code{signature(object = "TreeSet", value = "character")}: replaces slot \code{settype}.}
    \item{treeInfo}{\code{signature(object = "TreeSet")}: extracts UserInfo from \code{\link{ROOT}} trees.}
    \item{treeNames}{\code{signature(object = "TreeSet")}: extracts slot \code{treenames}.}
	 }
}
\author{Christian Stratowa}
\seealso{
derived classes \code{\linkS4class{SchemeTreeSet}}, \code{\linkS4class{DataTreeSet}},
 \code{\linkS4class{ExprTreeSet}}, \code{\linkS4class{CallTreeSet}}.
}
\examples{
showClass("TreeSet")
}
\keyword{classes}
