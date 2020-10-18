\name{exprs-methods}
\docType{methods}
\alias{exprs-methods}
\alias{exprs}
\alias{exprs<-}
\title{Get/Set Expression Values}
\description{
  Get/set expression values from/for class \code{ExprTreeSet}.

\emph{Usage}

\code{exprs(object)} \cr
\code{exprs(object, treenames = NULL) <- value}
}
\arguments{
  \item{object}{object of class \code{\linkS4class{ExprTreeSet}}.}
  \item{treenames}{\code{character} vector containing optional tree names to be used as subset.}
  \item{value}{\code{data.frame} containing expression values.}
}
\details{
 Get the expression values from slot \code{data} or set slot \code{data} to \code{value}.

 Method \code{exprs} returns the expression values from slot \code{data} as \code{data.frame},
 while replacement method \code{exprs<-} allows to replace slot \code{data} with a \code{data.frame}.

 In order to create an \code{ExprTreeSet} containing only a subset of slot \code{data}, first
 export slot \code{data} using method \code{exprs}, create a \code{character} vector containing
 only \code{treenames} to be used in the subset, and then use replacement method \code{exprs<-} to
 replace slot \code{data} with the subset. Slots \code{treenames} and \code{numtrees} will be
 updated automatically.

 Note: When creating \code{character} vector \code{treenames} it is sufficient to use the name
 part of the tree name w/o the extension.

 Note: If you do not want to replace your current object, create first a copy of type
 \code{ExprTreeSet} by simply writing \code{newobj <- oldobj}, and use \code{newobj} for
 replacement. This is important since \code{exprs<-} does also update slots \code{treenames}
 and \code{numtrees} as already mentioned.
}
\author{Christian Stratowa}
\seealso{\code{\link{pvalData}}, \code{\link{presCall}}}
\examples{
\dontrun{
## load ROOT scheme file and ROOT data file
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.test3 <- root.data(scheme.test3, paste(path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))

## create an ExprTreeSet
data.rma <- rma(data.test3,"tmp_TestRMA",tmpdir="",background="pmonly",normalize=TRUE,verbose=FALSE)

## get expression values
value <- exprs(data.rma)

## selected treenames only
treenames <- c("TestA2", "TestB1")

## make a copy of your object if you do not want to replace it
subset.rma <- data.rma

## replace slot data with subset
exprs(subset.rma, treenames) <- value
str(subset.rma)
}
}
\keyword{methods}
