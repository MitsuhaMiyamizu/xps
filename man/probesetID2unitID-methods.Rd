\name{probesetID2unitID-methods}
\docType{methods}
\alias{probesetID2unitID-methods}
\alias{probesetID2unitID}
\alias{transcriptID2unitID-methods}
\alias{transcriptID2unitID}
\alias{unitID2probesetID-methods}
\alias{unitID2probesetID}
\alias{unitID2transcriptID-methods}
\alias{unitID2transcriptID}
\title{Conversion between Probeset IDs and UnitIDs}
\description{
 Convert probeset IDs and transcript IDs to internal UNIT_IDs and vice verse. 

\emph{Usage}

\code{probesetID2unitID(object, probesetID = NULL, as.list = TRUE)}

\code{transcriptID2unitID(object, transcriptID = NULL, as.list = TRUE)}

\code{unitID2probesetID(object, unitID = NULL, as.list = TRUE)}

\code{unitID2transcriptID(object, unitID = NULL, as.list = TRUE)}
}
\arguments{
  \item{object}{Object of class \code{"SchemeTreeSet"} or \code{"DataTreeSet"}.}
  \item{probesetID}{optional vector of probeset IDs.}
  \item{transcriptID}{optional vector of transcript IDs.}
  \item{unitID}{optional vector of UNIT_IDs.}
  \item{as.list}{if TRUE a \code{list} will be returned (default is \code{data.frame}).}
}
\details{
 Functions \code{probesetID2unitID} and \code{transcriptID2unitID} return the UNIT_ID(s) for all or selected
 probeset IDs or transcript IDs, respectively.  

 Conversely, functions \code{unitID2probesetID} and \code{unitID2transcriptID} return the probeset IDs
 or transcript IDs, respectively, for all or selected UNIT_IDs.
.  
 For expression arrays the functions for probeset IDs and transcript IDs return identical IDs. 
.  
 For exon arrays the functions for probeset IDs and transcript IDs return the probeset_id(s) or
 transcript_cluster_id(s), respectivley. 

 By default a \code{list} is returned, however for \code{as.list=FALSE} a character vector of IDs is returned.
}
\value{
  A \code{list} or character \code{vector}.
}
\author{Christian Stratowa}
\seealso{\code{\link{unitID2transcriptID}}, \code{\link{unitID2probesetID}}}
\examples{
## load ROOT scheme file
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))

## unitNames not attached
id <- unitID2probesetID(scheme.test3, c(2,34,229))
id

## unitNames attached
scheme.test3 <- attachUnitNames(scheme.test3)
id <- probesetID2unitID(scheme.test3, c("PA1178_oprH_at","AFFX-Bt_eIF-4E_3_at","100084_at"))
id
scheme.test3 <- removeUnitNames(scheme.test3)

rm(scheme.test3)
gc()
}
\keyword{methods}
