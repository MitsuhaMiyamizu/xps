\name{getProbeInfo}
\alias{getProbeInfo}
\title{Get Probe Information}
\description{
  Get GeneChip probe information from root scheme file.
}
\usage{
getProbeInfo(rootfile)
}
\arguments{
  \item{rootfile}{name of ROOT scheme file, including full path.}
}
\details{
 Extracts GeneChip probe information directly from \code{\link{ROOT}} scheme file \code{rootfile}.
}
\value{
a \code{list} with parameters:
  \item{nrows }{physical number of rows in the array.}
  \item{ncols }{physical number of columns in the array.}
  \item{nprobes }{number of probes on the array.}
  \item{ncontrols }{number of controls on the array.}
  \item{ngenes }{number of genes on the array.}
  \item{nunits }{number of units on the array.}
  \item{nprobesets }{umber of probesets on the array.}
  \item{naffx }{number of AFFX controls on the array.}
}
\author{Christian Stratowa}
\examples{
getProbeInfo(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
}
\keyword{manip}
