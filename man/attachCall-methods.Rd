\name{attachCall-methods}
\docType{methods}
\alias{attachCall-methods}
\alias{attachCall}
\alias{attachPVal-methods}
\alias{attachPVal}
\alias{removeCall-methods}
\alias{removeCall}
\alias{removePVal-methods}
\alias{removePVal}
\title{Attach/Remove Detecion Call Measures}
\description{
 Attach/remove detection call and detection p-value to/from \code{\link{CallTreeSet}}.

\emph{Usage}

\code{attachCall(object, treenames = "*")}

\code{attachPVal(object, treenames = "*")}

\code{removeCall(object)}

\code{removePVal(object)}
}
\arguments{
  \item{object}{Object of class \code{"CallTreeSet"}.}
  \item{treenames}{Object of class \code{"list"} representing the names of the 
   \code{\link{ROOT}} call trees.}
}
\details{
 By default detection calls will be saved in class \code{\link{CallTreeSet}} in slots
 \code{data} and \code{detcall}, respectively, since usually the \code{data.frame}s obtained as
 result of e.g. \code{\link{mas5.call}} are of reasonable size. However, when computing many
 arrays, especially exon arrays at probeset levels, it may be better to compute detection
 calls with slot \code{add.data=FALSE} thus avoiding memory problems. In this case, functions
 \code{attachCall} and  \code{attachPVal} allow to  fill slots \code{detcall} and \code{data},
 respectively, on demand.

 \code{attachCall} exports detection calls from call trees from \code{\link{ROOT}} call
 file and and saves as data.frame \code{detcall}. \code{treenames} is a vector of
 tree names to attach; for \code{treenames="*"} all trees from slot \code{treenames}
 will be exported and detection calls attached as data.frame \code{detcall}.

 \code{attachPVal} exports detection p-values from call trees from \code{\link{ROOT}} call
 file and and saves as data.frame \code{data}. \code{treenames} is a vector of
 tree names to attach; for \code{treenames="*"} all trees from slot \code{treenames}
 will be exported and detection p-values attached as data.frame \code{data}.

 \code{removeCall} removes detection calls from \code{\link{CallTreeSet}} and replaces
 data.frame \code{detcall} with an empty data.frame of dim(0,0).

 \code{removePVal} removes detection p-values from \code{\link{CallTreeSet}} and replaces
 data.frame \code{data} with an empty data.frame of dim(0,0).
}
\value{
  A \code{\link{CallTreeSet}} object.
}
\author{Christian Stratowa}
\note{
 Do not use \code{attachCall} and \code{attachPVal}  unless you know that your computer has
 sufficient RAM, especially when using exon arrays. It may be advisible to use a subset of 
\code{treenames} only.
}
\seealso{\code{\link{attachExpr}}, \code{\link{removeExpr}}}
\examples{
## first, load ROOT scheme file and ROOT data file
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.test3 <- root.data(scheme.test3, paste(path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))

## MAS5 detection call
call.mas5 <- mas5.call(data.test3,"tmp_Test3Call0",tmpdir="",add.data=FALSE,verbose=FALSE)

## attach data
call.mas5 <- attachPVal(call.mas5)
call.mas5 <- attachCall(call.mas5)

## get data.frames
pval.mas5 <- pvalData(call.mas5)
pres.mas5 <- presCall(call.mas5)
head(pval.mas5)
head(pres.mas5)

## remove data
call.mas5 <- removePVal(call.mas5)
call.mas5 <- removeCall(call.mas5)

rm(scheme.test3, data.test3)
gc()
}
\keyword{methods}
