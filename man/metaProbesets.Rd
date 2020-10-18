\name{metaProbesets}
\alias{metaProbesets}
\title{Create MetaProbeset File for APT}
\description{
 Create MetaProbeset File for APT function \dQuote{apt-probeset-summarize}.
}
\usage{
metaProbesets(xps.scheme, infile = character(0), outfile = character(0), exonlevel="metacore")
}
\arguments{
  \item{xps.scheme}{exon \code{SchemeTreeSet}.}
  \item{infile}{Name of file containing exon transcript\_cluster\_ids.}
  \item{outfile}{Name of resulting file containing meta probeset definitions.}
  \item{exonlevel}{exon annotation level determining which probes should be used.}
}
\details{
 This function allows to create a metaprobeset file for APT function \dQuote{apt-probeset-summarize}
 to be used with option \dQuote{-m}. The \code{infile} must contain exon transcript\_cluster\_ids,
 one per line, e.g. one can export the \code{rownames(data.rma)} \cr
 The resulting file may be useful if you want to compare results created with \code{xps} to results
 created with APT function \dQuote{apt-probeset-summarize}.
}
\value{
 None.
}
\author{Christian Stratowa}
\examples{
\dontrun{
## first, load ROOT exon scheme file:
scmdir <- "/Volumes/GigaDrive/CRAN/Workspaces/Schemes"
scheme.exon <- root.scheme(paste(scmdir,"Scheme_HuEx10stv2r2_na25.root",sep="/"))

metaProbesets(scheme.exon,"metacore.txt","metacoreList.mps","metacore")
}
}
\keyword{manip}
