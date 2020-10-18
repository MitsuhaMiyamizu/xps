\name{intensity2GCplot-methods}
\docType{methods}
\alias{intensity2GCplot-methods}
\alias{intensity2GCplot}
\title{Boxlot of Probe Intensities Stratified by GC Content.}
\description{
 Creates a boxplot of probe intensities stratified by GC content.

\emph{Usage}

\code{
intensity2GCplot(x,
                 treename, 
                 which   = "",
                 transfo = log2,
                 range   = 0,
                 col     = c("lightblue", "darkblue"),
                 ...) 
}
}
\arguments{
  \item{x}{object of class \code{\link{DataTreeSet}}.}
  \item{treename}{\code{character} vector, tree name containing intensities.}
  \item{which}{type of probes to be used, for details see \code{\link{validData}}.}
  \item{transfo}{a valid function to transform the data, usually \dQuote{log2}, or \dQuote{0}.}
  \item{range}{determines how far the plot whiskers extend out from the box.}
  \item{col}{color pair to be used by function colorRampPalette.}
  \item{\dots}{optional arguments to be passed to \code{intensity2GCplot}.}
}
\details{
 Creates a boxplot of probe intensities for \code{treename} stratified by GC content for an object
 of class \code{\link{DataTreeSet}}.  
}
\author{Christian Stratowa}
\note{
 G/C content must first be attached to class \code{\link{DataTreeSet}} using method 
 \code{\link{attachProbeContentGC}}. It is also recommended to attach the probe mask using
 method \code{\link{attachMask}}. 
}
\seealso{\code{\link{plotIntensity2GC}}}
\examples{
## load existing ROOT scheme file and ROOT data file
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.test3 <- root.data(scheme.test3, paste(path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))

## need to attach probe G/C content and optionally mask
data.test3 <- attachProbeContentGC(data.test3)
data.test3 <- attachMask(data.test3)

if (interactive()) {
intensity2GCplot(data.test3, treename = "TestA1.cel", which="mm")
}

## optionally remove probe G/C content and mask to free memory
data.test3 <- removeMask(data.test3)
data.test3 <- removeProbeContentGC(data.test3)
}
\keyword{methods}
