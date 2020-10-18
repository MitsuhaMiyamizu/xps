\name{prefilter}
\alias{prefilter}
\alias{xpsPreFilter-methods}
\alias{xpsPreFilter}
\title{Function for Applying a PreFilter to an ExprTreeSet}
\description{
 This function applies a \code{\link{PreFilter}} to an \code{\link{ExprTreeSet}}.
}
\usage{
prefilter(xps.expr,
          filename   = character(0),
          filedir    = getwd(),
          filter     = NULL,
          minfilters = 999,
          logbase    = "log2",
          treename   = "PreFilter",
          xps.call   = NULL,
          verbose    = TRUE)


xpsPreFilter(object, ...)
}
\arguments{
  \item{xps.expr}{object of class \code{ExprTreeSet}.}
  \item{filename}{file name of ROOT filter file.}
  \item{filedir}{system directory where ROOT filter file should be stored.}
  \item{filter}{object of class \code{PreFilter}.}
  \item{minfilters}{minimum number of initialized filter methods to satisfy (default is all filters).}
  \item{logbase}{convert data to logarithm of base: \code{"0"}, \code{"log"}, \code{"log2"} (default), \code{"log10"}}
  \item{treename}{tree name to be used in ROOT filter file.}
  \item{xps.call}{optional object of class \code{CallTreeSet}.}
  \item{verbose}{logical, if \code{TRUE} print status information.}
  \item{object}{object of class \code{ExprTreeSet}.}
  \item{\dots}{same arguments as function \code{prefilter}.}
}
\details{
 This function applies the different filters initialized with constructor \code{\link{PreFilter}}
 to the \code{\link{ExprTreeSet}} \code{xps.expr}.

 Slot \code{minfilters} determines the minimum number of initialized filters, which must be satisfied
 so that the mask is set to \code{flag=1}. For \code{minfilters=1} at least one filter must be 
 satisfied, equivalent to logical \sQuote{OR}; for \code{minfilters=999} all filters must be 
 satisfied, equivalent to logical \sQuote{AND}.

 If method \code{\link{callFilter}} was initialized with constructor \code{\link{PreFilter}}
 then \code{\link{CallTreeSet}} \code{xps.call} must be supplied, usually created with function
 \code{\link{mas5.call}}.
}
\value{
 A \code{\linkS4class{FilterTreeSet}}
}
\author{Christian Stratowa}
\seealso{\code{\link{PreFilter}}, \code{\link{unifilter}}}
\examples{
\dontrun{
## first, load ROOT scheme file and ROOT data file
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.test3 <- root.data(scheme.test3, paste(path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))

## second, create an ExprTreeSet
data.rma <- rma(data.test3,"tmp_TestRMA",tmpdir="",background="pmonly",normalize=TRUE,verbose=FALSE)
## note: do not copy/paste this code, it is necessary only because R CMD check fails since it does not find tmp_TestRMA.root:
data.rma@rootfile <- paste(path.package("xps"),"rootdata/tmp_Test3RMA.root",sep="/")
data.rma@filedir  <- paste(path.package("xps"),"rootdata",sep="/")

## third, construct a PreFilter
prefltr <- PreFilter(mad=c(0.5,0.01),lothreshold=c(6.0,0.02,"mean"),hithreshold=c(10.5,80.0,"percent"))

## finally, create a FilterTreeSet
rma.pfr <- prefilter(data.rma,"tmp_Test3Prefilter",getwd(),prefltr,2,verbose=FALSE)
str(rma.pfr)
}
}
\keyword{manip}
