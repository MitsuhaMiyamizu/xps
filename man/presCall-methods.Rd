\name{presCall-methods}
\docType{methods}
\alias{presCall-methods}
\alias{presCall}
\alias{presCall<-}
\alias{pvalData-methods}
\alias{pvalData}
\alias{pvalData<-}
\title{Get/Set Present Call Values}
\description{
  Get/set present call values from/for class \code{CallTreeSet}.

\emph{Usage}

\code{presCall(object)} \cr
\code{presCall(object, treenames = NULL) <- value} \cr
 \cr
\code{pvalData(object)} \cr
\code{pvalData(object, treenames = NULL) <- value}
}
\arguments{
  \item{object}{object of class \code{\linkS4class{CallTreeSet}}.}
  \item{treenames}{\code{character} vector containing optional tree names to be used as subset.}
  \item{value}{\code{data.frame} containing present call values.}
}
\details{
 Get the p-values from slot \code{data} or present calls from slot \code{detcall},
 or set slot \code{data} or \code{detcall}, respectively, to \code{value}.

 Method \code{presCall} returns the present calls from slot \code{detcall} as \code{data.frame},
 while replacement method \code{presCall<-} allows to replace slot \code{detcall} with a
 \code{data.frame}.

 Method \code{pvalData} returns the p-values from slot \code{data} as \code{data.frame},
 while replacement method \code{pvalData<-} allows to replace slot \code{data} with a
 \code{data.frame}.

 In order to create an \code{CallTreeSet} containing only a subset of e.g. slot \code{data}, first
 export slot \code{data} using method \code{pvalData}, create a \code{character} vector containing
 only \code{treenames} to be used in the subset, and then use replacement method \code{pvalData<-}
 to replace slot \code{data} with the subset. Slots \code{treenames} and \code{numtrees} will be
 updated automatically for \code{pvalData<-} but not for \code{presCall<-}.

 Note: When creating \code{character} vector \code{treenames} it is sufficient to use the name
 part of the tree name w/o the extension.

 Note: If you do not want to replace your current object, create first a copy of type
 \code{CallTreeSet} by simply writing \code{newobj <- oldobj}, and use \code{newobj} for
 replacement.
}
\author{Christian Stratowa}
\seealso{\code{\link{exprs}}}
\examples{
\dontrun{
## load ROOT scheme file and ROOT data file
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.test3 <- root.data(scheme.test3, paste(path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))

## create an CallTreeSet
call.mas5 <- mas5.call(data.test3,"tmp_TestMAS5Call",tmpdir="",verbose=FALSE)

## get p-values
value <- pvalData(call.mas5)

## selected treenames only
treenames <- c("TestA2", "TestB1")

## make a copy of your object if you do not want to replace it
subset.call <- call.mas5

## replace slot data with subset
exprs(subset.call, treenames) <- value
str(subset.call)
}
}
\keyword{methods}
