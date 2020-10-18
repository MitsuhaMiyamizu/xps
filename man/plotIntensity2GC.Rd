\name{plotIntensity2GC}
\alias{plotIntensity2GC}
\title{Boxlot of Probe Intensities Stratified by GC Content for Device.}
\description{
 Creates a boxplot of probe intensities stratified by GC content for the selected device.
}
\usage{
plotIntensity2GC(x,
                 treename, 
                 which   = "",
                 transfo = log2,
                 range   = 0,
                 col     = c("lightblue", "darkblue"),
                 dev     = "screen",
                 outfile = "Intensity2GCPlot",
                 w       = 540,
                 h       = 540,
                 ...) 
}
\arguments{
  \item{x}{object of class \code{\link{DataTreeSet}}.}
  \item{treename}{\code{character} vector, tree name used for intensities.}
  \item{which}{type of probes to be used, for details see \code{\link{validData}}.}
  \item{transfo}{a valid function to transform the data, usually \dQuote{log2}, or \dQuote{0}.}
  \item{range}{determines how far the plot whiskers extend out from the box.}
  \item{col}{color pair to be used by function colorRampPalette.}
  \item{dev}{graphics device to plot to, i.e. one of \dQuote{screen}, \dQuote{jpeg},\dQuote{png},
   \dQuote{pdf} or \dQuote{ps}.}
  \item{outfile}{the name of the output file.}
  \item{w}{the width of the device in pixels.}
  \item{h}{the height of the device in pixels.}
  \item{\dots}{optional arguments to be passed to \code{plotIntensity2GC}.}
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
\seealso{\code{\link{intensity2GCplot}}}
\keyword{device}
