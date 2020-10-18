\name{QualTreeSet-class}
\docType{class}
\alias{QualTreeSet-class}
\alias{QualTreeSet}
\alias{borders,QualTreeSet-method}
\alias{borders}
\alias{borderplot,QualTreeSet-method}
\alias{coiplot,QualTreeSet-method}
\alias{image,QualTreeSet-method}
\alias{NUSE,QualTreeSet-method}
\alias{nuseplot,QualTreeSet-method}
\alias{qualOption,QualTreeSet-method}
\alias{qualOption}
\alias{qualOption<-,QualTreeSet,character-method}
\alias{qualOption<-}
\alias{qualType,QualTreeSet-method}
\alias{qualType}
\alias{qualType<-,QualTreeSet,character-method}
\alias{qualType<-}
\alias{residuals,QualTreeSet-method}
\alias{residuals}
\alias{RLE,QualTreeSet-method}
\alias{rleplot,QualTreeSet-method}
\alias{weights,QualTreeSet-method}
\alias{weights}
\alias{xpsRNAdeg,QualTreeSet-method}

\title{Class QualTreeSet}
\description{This class provides the link to the \code{\link{ROOT}} quality control file and 
the \code{\link{ROOT}} trees contained therein. It extends class \code{\link{ProcesSet}}.}
\section{Objects from the Class}{
Objects are created using functions \code{\link{qualify}}, \code{\link{fitQC}}, or the specialized 
functions \code{\link{qualify.rlm}}, \code{\link{fitRLM}} or \code{\link{rmaPLM}}.
}
\section{Slots}{
	 \describe{
    \item{\code{qualopt}:}{Object of class \code{"character"} representing the quality control option, 
     i.e. \sQuote{raw}, \sQuote{adjusted}, \sQuote{normalized} or \sQuote{all}.}
    \item{\code{qualtype}:}{Object of class \code{"character"} representing the quality control type, 
     i.e. \sQuote{rlm}.}
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
    \item{borderplot}{\code{signature(x = "QualTreeSet")}: creates a boxplot of positive and 
     negative border elements.}
    \item{borders}{\code{signature(object = "QualTreeSet")}: exports border trees from 
     \code{\link{ROOT}} quality control file as data.frame \code{data}.}
    \item{coiplot}{\code{signature(x = "QualTreeSet")}: creates a Center-of-Intensity-plot for
     positive and negative feature intensities.}
    \item{image}{\code{signature(x = "QualTreeSet")}: creates a pseudo image for each quality 
     control tree, i.e. residual images.}
    \item{NUSE}{\code{signature(x = "QualTreeSet")}: plot Normalized Unscaled Standard Errors, or return stats, values.}
    \item{nuseplot}{\code{signature(x = "QualTreeSet")}: creates a NUSE-plot.}
    \item{qualOption}{\code{signature(object = "QualTreeSet")}: extracts slot \code{qualopt}.}
    \item{qualOption<-}{\code{signature(object = "QualTreeSet", value = "character")}: replaces slot \code{qualopt}.}
    \item{qualType}{\code{signature(object = "QualTreeSet")}: extracts slot \code{qualtype}.}
    \item{qualType<-}{\code{signature(object = "QualTreeSet", value = "character")}: replaces slot \code{qualtype}.}
    \item{residuals}{\code{signature(object = "QualTreeSet")}: exports residuals from the
     residuals trees of the \code{\link{ROOT}} quality control file as data.frame \code{data}.}
    \item{RLE}{\code{signature(x = "QualTreeSet")}: plot Relative Log Expression, or return stats, values.}
    \item{rleplot}{\code{signature(x = "QualTreeSet")}: creates a RLE-plot.}
    \item{weights}{\code{signature(object = "QualTreeSet")}: exports weights from the
     residuals trees of the \code{\link{ROOT}} quality control file as data.frame \code{data}.}
    \item{xpsRNAdeg}{\code{signature(x = "QualTreeSet")}: list with parameters for RNA degradation.}
	 }
}
\author{Christian Stratowa}
\seealso{
related classes \code{\linkS4class{DataTreeSet}}, \code{\linkS4class{CallTreeSet}}, \code{\linkS4class{ExprTreeSet}}.
}
\examples{
showClass("QualTreeSet")
}
\keyword{classes}
