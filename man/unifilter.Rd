\name{unifilter}
\alias{unifilter}
\alias{xpsUniFilter-methods}
\alias{xpsUniFilter}
\title{Function for Applying an UniFilter to an ExprTreeSet}
\description{
 This function applies an \code{\link{UniFilter}} to an \code{\link{ExprTreeSet}}.
}
\usage{
unifilter(xps.expr,
          filename   = character(0),
          filedir    = getwd(),
          filter     = NULL,
          minfilters = 999,
          logbase    = "log2",
          group      = character(0),
          treename   = "UniTest",
          xps.fltr   = NULL,
          xps.call   = NULL,
          update     = FALSE,
          verbose    = TRUE)


xpsUniFilter(object, ...)
}
\arguments{
  \item{xps.expr}{object of class \code{ExprTreeSet}.}
  \item{filename}{file name of ROOT filter file.}
  \item{filedir}{system directory where ROOT filter file should be stored.}
  \item{filter}{object of class \code{UniFilter}.}
  \item{minfilters}{minimum number of initialized filter methods to satisfy (default is all filters).}
  \item{logbase}{convert data to logarithm of base: \code{"0"}, \code{"log"}, \code{"log2"} (default), \code{"log10"}}
  \item{group}{a character vector assigning the trees of \code{xps.expr} to one of two groups.}
  \item{treename}{tree name to be used in ROOT filter file.}
  \item{xps.fltr}{optional object of class \code{FilterTreeSet}.}
  \item{xps.call}{optional object of class \code{CallTreeSet}.}
  \item{update}{logical. If \code{TRUE} the existing ROOT filter file \code{filename} will be updated.}
  \item{verbose}{logical, if \code{TRUE} print status information.}
  \item{object}{object of class \code{ExprTreeSet}.}
  \item{\dots}{same arguments as function \code{unifilter}.}
}
\details{
 This function applies the different filters initialized with constructor \code{\link{UniFilter}}
 to the \code{\link{ExprTreeSet}} \code{xps.expr}.

 Slot \code{minfilters} determines the minimum number of initialized filters, which must be satisfied
 so that the mask is set to \code{flag=1}. For \code{minfilters=1} at least one filter must be 
 satisfied, equivalent to logical \sQuote{OR}; for \code{minfilters=999} all filters must be 
 satisfied, equivalent to logical \sQuote{AND}.

 If pre-filtering should be done before applying function \code{unifilter} then a
 \code{\link{FilterTreeSet}} \code{xps.fltr} must be supplied, created with function
 \code{\link{prefilter}}.

 If method \code{\link{callFilter}} was initialized with constructor \code{\link{UniFilter}}
 then \code{\link{CallTreeSet}} \code{xps.call} must be supplied, usually created with function
 \code{\link{mas5.call}}.
}
\value{
 An \code{\linkS4class{AnalysisTreeSet}}
}
\author{Christian Stratowa}
\note{
 Internally, slot \code{group} will be converted to integer values using \code{as.integer(as.factor(group))},
 thus \code{group=c("GrpA","GrpA","GrpB","GrpB")} will result in a fold-change of \code{fc=mean(GrpB)/mean(GrpA)}.
}
\seealso{\code{\link{UniFilter}}, \code{\link{prefilter}}}
\examples{
\dontrun{
## first, load ROOT scheme file and ROOT data file
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.test3 <- root.data(scheme.test3, paste(path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))

## second, create an ExprTreeSet
data.rma <- rma(data.test3,"tmp_Test3_RMA",tmpdir="",background="pmonly",normalize=TRUE,verbose=FALSE)
## note: do not copy/paste this code, it is necessary only because R CMD check fails since it does not find tmp_Test3_RMA.root:
data.rma@rootfile <- paste(path.package("xps"),"rootdata/tmp_Test3RMA.root",sep="/")
data.rma@filedir  <- paste(path.package("xps"),"rootdata",sep="/")

## third, construct an UniFilter
unifltr <- UniFilter(unitest=c("t.test","two.sided","none",0,0.0,FALSE,0.95,TRUE),foldchange=c(1.3,"both"),unifilter=c(0.1,"pval"))

## finally, create an AnalysisTreeSet
rma.ufr <- unifilter(data.rma,"tmp_Test3Unifilter",getwd(),unifltr,group=c("GrpA","GrpA","GrpB","GrpB"),verbose=FALSE)
str(rma.ufr)
}
}
\keyword{manip}
