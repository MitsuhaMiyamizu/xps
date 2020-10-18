\name{pm-methods}
\docType{methods}
\alias{pm-methods}
\alias{mm-methods}
\alias{pm}
\alias{mm}
\title{Methods for accessing perfect matches and mismatches}
\description{
  Methods for accessing perfect match (PM) and mismatch (MM) probes.

\emph{Usage}

\code{pm(object, which = "pm", unitID = NULL, unittype = "transcript")}

\code{mm(object, which = "mm", unitID = NULL, unittype = "transcript")}
}
\arguments{
  \item{object}{object of class \code{DataTreeSet}.}
  \item{which}{type of perfect match or mismatch probes to be returned.}
  \item{unitID}{optional vector of UNIT_IDs.}
  \item{unittype}{\code{character} vector, \dQuote{transcript} or \dQuote{probeset}.}
}
\details{
For expression arrays all the perfect match (pm) or mismatch (mm) probes on the arrays 
the object represents are returned as data.frame.

For exon arrays, \code{pm} returns the probes of the different exon levels as data.frame, 
i.e. \code{which} can have one of the following values:
\tabular{lll}{
  \tab \code{core}:\tab probesets supported by RefSeq and full-length GenBank transcripts. \cr
  \tab \code{metacore}:\tab core meta-probesets. \cr
  \tab \code{extended}:\tab probesets with other cDNA support. \cr
  \tab \code{metaextended}:\tab extended meta-probesets. \cr
  \tab \code{full}:\tab probesets supported by gene predictions only. \cr
  \tab \code{metafull}:\tab full meta-probesets. \cr
  \tab \code{affx}:\tab standard AFFX controls.
}

For whole genome arrays, \code{pm} returns the probes of the different exon levels as data.frame, 
i.e. \code{which} can have one of the following values:
\tabular{lll}{
  \tab \code{core}:\tab probesets with category \sQuote{unique} and \sQuote{mixed}. \cr
  \tab \code{metacore}:\tab probesets with category \sQuote{unique} only. \cr
  \tab \code{affx}:\tab standard AFFX controls.
}

For exon/genome arrays, \code{mm} returns the background probes as data.frame, 
i.e. \code{which} is either \dQuote{genomic} or \dQuote{antigenomic}.
}
\value{
  A \code{\link{data.frame}}.
}
\author{Christian Stratowa}
\seealso{\code{\link{validData}}}
\examples{
## load existing ROOT scheme file and ROOT data file
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.test3 <- root.data(scheme.test3, paste(path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))

## need to attach scheme mask and probe intensities
data.test3 <- attachMask(data.test3)
data.test3 <- attachInten(data.test3)

pm <- pm(data.test3)
mm <- mm(data.test3)
head(pm)
head(mm)

## need to convert Affy ID to UNIT_ID first
id <- transcriptID2unitID(schemeSet(data.test3), transcriptID="100084_at", as.list=FALSE)
pm <- pm(data.test3, unitID=id)
mm <- mm(data.test3, unitID=id)
head(pm)
head(mm)

## optionally remove mask and data to free memory
data.test3 <- removeInten(data.test3)
data.test3 <- removeMask(data.test3)
}
\keyword{methods}
