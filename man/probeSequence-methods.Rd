\name{probeSequence-methods}
\docType{methods}
\alias{probeSequence-methods}
\alias{probeSequence}
\title{Get Probe Sequence}
\description{
 Get probe sequences for all or selected UNIT_IDs. 

\emph{Usage}

\code{probeSequence(object, unitID = NULL)}
}
\arguments{
  \item{object}{Object of class \code{"SchemeTreeSet"} or \code{"DataTreeSet"}.}
  \item{unitID}{optional vector of UNIT_IDs.}
}
\details{
 Function \code{probeSequence} returns a \code{data.frame} containing column \dQuote{ProbeSequence}
 for all or selected the UNIT_ID(s).  
}
\value{
  A \code{data.frame}.
}
\author{Christian Stratowa}
\seealso{\code{\link{probeContentGC}}}
\examples{
## load ROOT scheme file
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
scheme.test3 <- attachProbeSequence(scheme.test3)

## get UNIT_ID for probeset ID
id <- probesetID2unitID(scheme.test3, "100084_at")

## get GC content
seq <- probeSequence(scheme.test3, unitID=id)
head(seq)

scheme.test3 <- removeProbeSequence(scheme.test3)

rm(scheme.test3)
gc()
}
\keyword{methods}
