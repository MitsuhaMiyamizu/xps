\name{attachExpr-methods}
\docType{methods}
\alias{attachExpr-methods}
\alias{attachExpr}
\alias{removeExpr-methods}
\alias{removeExpr}
\title{Attach/Remove Expression Measures}
\description{
 Attach/remove expression levels to/from \code{\link{ExprTreeSet}}.

\emph{Usage}

\code{attachExpr(object, treenames = "*")}

\code{removeExpr(object)}
}
\arguments{
  \item{object}{Object of class \code{"ExprTreeSet"}.}
  \item{treenames}{Object of class \code{"list"} representing the names of the 
   \code{\link{ROOT}} expression trees.}
}
\details{
 By default expression levels will be saved in class \code{\link{ExprTreeSet}} as slot
 \code{data}, since usually the \code{data.frame} obtained as result of e.g. \code{rma}
 normalization is of reasonable size. However, when normalizing many arrays, especially exon
 arrays at probeset levels, it may be better to compute \code{rma} with slot \code{add.data=FALSE}
 thus avoiding memory problems. In this case, function \code{attachExpr} allows to  fill
 slot \code{data} on demand.

 \code{attachExpr} exports expression levels from expression trees from \code{\link{ROOT}}
 expression file and and saves as data.frame \code{data}. \code{treenames} is a vector of
 tree names to attach; for \code{treenames="*"} all trees from slot \code{treenames}
 will be exported and expression levels attached as data.frame \code{data}.

 \code{removeExpr} removes expression levels from \code{\link{ExprTreeSet}} and replaces
 data.frame \code{data} with an empty data.frame of dim(0,0).
}
\value{
  A \code{\link{ExprTreeSet}} object.
}
\author{Christian Stratowa}
\note{
 Do not use \code{attachExpr} unless you know that your computer has sufficient RAM,
 especially when using exon arrays. It may be advisible to use a subset of \code{treenames}
 only.
}
\seealso{\code{\link{attachCall}}, \code{\link{removeCall}}}
\examples{
## first, load ROOT scheme file and ROOT data file
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.test3 <- root.data(scheme.test3, paste(path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))

data.rma <- rma(data.test3,"tmp_Test3RMA0",tmpdir="",background="pmonly",normalize=TRUE,add.data=FALSE,verbose=FALSE)

## attach data
data.rma <- attachExpr(data.rma)

## get data.frame
expr.rma <- validData(data.rma)
head(expr.rma)

## remove data
data.rma <- removeExpr(data.rma)

rm(scheme.test3, data.test3)
gc()
}
\keyword{methods}
