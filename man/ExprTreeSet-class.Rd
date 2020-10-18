\name{ExprTreeSet-class}
\docType{class}
\alias{ExprTreeSet-class}
\alias{ExprTreeSet}
\alias{attachExpr,ExprTreeSet-method}
\alias{corplot,ExprTreeSet-method}
\alias{exprType,ExprTreeSet-method}
\alias{exprType}
\alias{exprType<-,ExprTreeSet,character-method}
\alias{exprType<-}
\alias{exprs,ExprTreeSet-method}
\alias{exprs<-,ExprTreeSet,data.frame-method}
\alias{madplot,ExprTreeSet-method}
\alias{mvaplot,ExprTreeSet-method}
\alias{normType,ExprTreeSet-method}
\alias{normType}
\alias{normType<-,ExprTreeSet,character-method}
\alias{normType<-}
\alias{nuseplot,ExprTreeSet-method}
\alias{pcaplot,ExprTreeSet-method}
\alias{rleplot,ExprTreeSet-method}
\alias{removeExpr,ExprTreeSet-method}
\alias{se.exprs,ExprTreeSet-method}
\alias{se.exprs}
\alias{validExpr,ExprTreeSet-method}
\alias{validSE,ExprTreeSet-method}
\alias{xpsNormalize,ExprTreeSet-method}
\alias{xpsPreFilter,ExprTreeSet-method}
\alias{xpsUniFilter,ExprTreeSet-method}

\title{Class ExprTreeSet}
\description{This class provides the link to the \code{\link{ROOT}} expression file and 
the \code{\link{ROOT}} trees contained therein. It extends class \code{\link{ProcesSet}}.}
\section{Objects from the Class}{
Objects are created using functions \code{\link{express}}, \code{\link{summarize}} 
or \code{\link{normalize}}, or the specialized functions \code{\link{rma}}, \code{\link{mas5}} 
or \code{\link{mas4}}.
}
\section{Slots}{
	 \describe{
    \item{\code{exprtype}:}{Object of class \code{"character"} representing the exression type, 
     i.e. \sQuote{rma}, \sQuote{mas5}, \sQuote{mas4} or \sQuote{custom}.}
    \item{\code{normtype}:}{Object of class \code{"character"} representing the normalization 
     type, i.e. \sQuote{mean}, \sQuote{median}, \sQuote{lowess},\sQuote{supsmu}.}
    \item{\code{scheme}:}{Object of class \code{"SchemeTreeSet"} providing access to 
     \code{\link{ROOT}} scheme file.}
    \item{\code{data}:}{Object of class \code{"data.frame"}. The data.frame can contain 
     the data (e.g. expression levels) stored in \code{\link{ROOT}} data trees.}
    \item{\code{params}:}{Object of class \code{"list"} representing relevant parameters.}
    \item{\code{setname}:}{Object of class \code{"character"} representing the name to the
     \code{\link{ROOT}} file subdirectoy where the \code{\link{ROOT}} data trees are stored,
     usually \sQuote{PreprocesSet}.}
    \item{\code{settype}:}{Object of class \code{"character"} describing the type of treeset
     stored in \code{setname}, usually \sQuote{preprocess}.}
    \item{\code{rootfile}:}{Object of class \code{"character"} representing the name 
     of the \code{\link{ROOT}} data file, including full path.}
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
    \item{attachExpr}{\code{signature(object = "ExprTreeSet")}: exports expression trees from 
     \code{\link{ROOT}} expression file and and saves as data.frame \code{data}.}
    \item{corplot}{\code{signature(x = "ExprTreeSet")}: creates a correlation heat map.}
    \item{exprType}{\code{signature(object = "ExprTreeSet")}: extracts slot \code{exprtype}.}
    \item{exprType<-}{\code{signature(object = "ExprTreeSet", value = "character")}: replaces slot \code{exprtype}.}
    \item{exprs}{\code{signature(object = "ExprTreeSet")}: extracts the expression data.frame.}
    \item{exprs<-}{\code{signature(object = "ExprTreeSet", value = "data.frame")}: replaces the expression data.frame.}
    \item{madplot}{\code{signature(x = "ExprTreeSet")}: creates a false color display of between arrays distances.}
    \item{mvaplot}{\code{signature(x = "ExprTreeSet")}: creates an MvA-plot.}
    \item{normType}{\code{signature(object = "ExprTreeSet")}: extracts slot \code{normtype}.}
    \item{normType<-}{\code{signature(object = "ExprTreeSet", value = "character")}: replaces slot \code{normtype}.}
    \item{nuseplot}{\code{signature(x = "ExprTreeSet")}: creates a NUSE-plot.}
    \item{pcaplot}{\code{signature(x = "ExprTreeSet")}: plots first two principal components of PCA.}
    \item{rleplot}{\code{signature(x = "ExprTreeSet")}: creates a RLE-plot.}
    \item{removeExpr}{\code{signature(object = "ExprTreeSet")}: replaces data.frame \code{data} 
     with an empty data.frame of dim(0,0).}
    \item{se.exprs}{\code{signature(object = "ExprTreeSet")}: extracts the standard deviation data.frame.}
    \item{validExpr}{\code{signature(object = "ExprTreeSet")}: extracts a subset of columns 
     from data.frame \code{data}.}
    \item{validSE}{\code{signature(object = "ExprTreeSet")}: extracts data columns 
     from data.frame \code{se.exprs}.}
    \item{xpsNormalize}{\code{signature(object = "ExprTreeSet")}: applies normalization methods.}
    \item{xpsPreFilter}{\code{signature(object = "ExprTreeSet")}: applies prefiltering methods.}
    \item{xpsUniFilter}{\code{signature(object = "ExprTreeSet")}: applies unifiltering methods.}
	 }
}
\author{Christian Stratowa}
\seealso{
related classes \code{\linkS4class{DataTreeSet}}, \code{\linkS4class{CallTreeSet}}, \code{\linkS4class{QualTreeSet}}.
}
\examples{
showClass("ExprTreeSet")
}
\keyword{classes}
