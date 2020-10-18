\name{attachProbe-methods}
\docType{methods}
\alias{attachProbe-methods}
\alias{attachProbe}
\alias{attachProbe-methods}
\alias{attachProbeContentGC}
\alias{attachProbeContentGC-methods}
\alias{attachProbeSequence}
\alias{removeProbeSequence-methods}
\alias{removeProbe}
\alias{removeProbeContentGC-methods}
\alias{removeProbeContentGC}
\alias{removeProbeSequence-methods}
\alias{removeProbeSequence}
\title{Attach/Remove Probe Sequence and/or GC Content}
\description{
 Attach/remove probe sequence and/or GC content to/from \code{\link{SchemeTreeSet}} or to slot \code{scheme}
 of \code{\link{DataTreeSet}}.

\emph{Usage}

\code{attachProbe(object, varlist)}

\code{attachProbeContentGC(object)}

\code{attachProbeSequence(object)}

\code{removeProbe(object)}

\code{removeProbeContentGC(object)}

\code{removeProbeSequence(object)}
}
\arguments{
  \item{object}{Object of class \code{"SchemeTreeSet"} or \code{"DataTreeSet"}.}
  \item{varlist}{names of probe tree leaves to import to slot \code{probe}.}
}
\details{
 Function \code{attachProbe} exports leaves from probe tree of \code{\link{ROOT}} scheme file
 and and saves the data as data.frame \code{probe} of slot \code{scheme}.

 Following \code{varlist} parameters are valid:
 \tabular{lll}{
  \tab \code{fPosition}:\tab probe interrogation position. \cr
  \tab \code{fSequence}:\tab probe sequence. \cr
  \tab \code{fNumberGC}:\tab number of G/C nucleotides in probe sequence. \cr
  \tab \code{fTm}:\tab probe melting temperature dependent on G/C number. \cr
  \tab \code{fIsAntisense}:\tab probe type (sense/antisense).
 }

 Function \code{attachProbeContentGC} saves \code{fNumberGC} in data.frame \code{probe} of
 \code{\link{SchemeTreeSet}} or in slot \code{scheme} of \code{\link{DataTreeSet}}.

 Function \code{attachProbeSequence} saves \code{fSequence} in data.frame \code{probe} of
 \code{\link{SchemeTreeSet}}.

 Function \code{removeProbe} removes \code{probe} data from \code{\link{SchemeTreeSet}} or from 
 slot \code{scheme} of \code{\link{DataTreeSet}} and replaces data.frame \code{probe} with an
 empty data.frame of dim(0,0). 
}
\value{
  A \code{\link{SchemeTreeSet}} object or \code{\link{DataTreeSet}}.
}
\author{Christian Stratowa}
\note{
 Do not use \code{attachProbe} unless you know that your computer has sufficient RAM,
 especially for exon array schemes.
}
\seealso{\code{\link{attachMask}}}
\examples{
## load existing ROOT scheme file
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
dim(chipProbe(scheme.test3))

scheme.test3 <- attachProbe(scheme.test3, varlist="fSequence:fNumberGC")
dim(chipProbe(scheme.test3))
head(chipProbe(scheme.test3))

scheme.test3 <- removeProbe(scheme.test3)
dim(chipProbe(scheme.test3))
}
\keyword{methods}
