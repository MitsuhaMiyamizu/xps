\name{symbol2unitID-methods}
\docType{methods}
\alias{symbol2unitID-methods}
\alias{symbol2unitID}
\alias{unitID2symbol-methods}
\alias{unitID2symbol}
\title{Conversion between Gene Symbols and UnitIDs}
\description{
 Convert gene symbols to internal UNIT_IDs and vice verse. 

\emph{Usage}

\code{symbol2unitID(object, symbol, unittype = "transcript", as.list = TRUE)}

\code{unitID2symbol(object, unitID, unittype = "transcript", as.list = TRUE)}
}
\arguments{
  \item{object}{Object of class \code{"SchemeTreeSet"} or \code{"DataTreeSet"}.}
  \item{symbol}{\code{character} vector of gene symbol(s).}
  \item{unitID}{vector of UNIT_IDs.}
  \item{unittype}{\code{character} vector, \dQuote{transcript} or \dQuote{probeset}.}
  \item{as.list}{if TRUE a \code{list} will be returned (default is \code{data.frame}).}
}
\details{
 Functions \code{symbol2unitID} and \code{unitID2symbol} returns the UNIT_ID(s) for selected gene symbols
 and vice verse.  

 For exon arrays the internal UNIT_ID(s) depend on \code{unittype}. 

 By default a \code{list} is returned, however for \code{as.list=FALSE} a character vector of IDs is returned.
}
\value{
  A \code{list} or character \code{vector}.
}
\author{Christian Stratowa}
\seealso{\code{\link{transcriptID2unitID}}, \code{\link{probesetID2unitID}}}
\examples{
## load ROOT scheme file
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))

## unitnames not attached
id <- symbol2unitID(scheme.test3, symbol="ACTB", as.list=TRUE)
id
id <- unitID2symbol(scheme.test3, unitID=274, as.list=TRUE)
id

## unitnames attached
scheme.test3 <- attachUnitNames(scheme.test3)
id <- symbol2unitID(scheme.test3, symbol="ACTB", as.list=TRUE)
id
id <- unitID2symbol(scheme.test3, unitID=274, as.list=TRUE)
id
scheme.test3 <- removeUnitNames(scheme.test3)

rm(scheme.test3)
gc()
}
\keyword{methods}
