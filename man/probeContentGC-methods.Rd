\name{probeContentGC-methods}
\docType{methods}
\alias{probeContentGC-methods}
\alias{probeContentGC}
\title{Get G/C Content for Probes}
\description{
 Get G/C content for all or selected UNIT_IDs. 

\emph{Usage}

\code{probeContentGC(object, which = "", unitID = NULL, unittype = "transcript")}
}
\arguments{
  \item{object}{Object of class \code{"SchemeTreeSet"} or \code{"DataTreeSet"}.}
  \item{which}{type of probes to be used, for details see \code{\link{validData}}.}
  \item{unitID}{optional vector of UNIT_IDs.}
  \item{unittype}{\code{character} vector, one of \dQuote{transcript}, \dQuote{probeset}.}
}
\details{
 Function \code{probeContentGC} returns a \code{data.frame} containing columns \dQuote{Mask}
 and \dQuote{ContentGC} for all or selected the UNIT_ID(s).  

 For exon arrays the type of UNIT_ID(s) depends on \code{unittype}. 
}
\value{
  A \code{data.frame}.
}
\author{Christian Stratowa}
\seealso{\code{\link{probeSequence}}}
\examples{
## load ROOT scheme file
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
scheme.test3 <- attachProbeContentGC(scheme.test3)

## get UNIT_ID for probeset IDs
id <- probesetID2unitID(scheme.test3, c("PA1178_oprH_at","AFFX-Bt_eIF-4E_3_at","100084_at"))

## get GC content
gc <- probeContentGC(scheme.test3, unitID=id)
head(gc)

scheme.test3 <- removeProbeContentGC(scheme.test3)

rm(scheme.test3)
gc()
}
\keyword{methods}
