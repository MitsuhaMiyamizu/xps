\name{plotProbeset}
\alias{plotProbeset}
\title{Plot of Probe Intensities for a Probeset for Device.}
\description{
 Creates a line plot of probe intensities for a probeset for the selected device.
}
\usage{
plotProbeset(x,
             unitID,
             unittype   = "transcript",
             which      = "pm",
             transfo    = log2,
             names      = "namepart",
             ylim       = NULL,
             col        = 1:6,
             lty        = 1:5,
             add.legend = FALSE,
             dev        = "screen",
             outfile    = "ProbesetPlot",
             w          = 540,
             h          = 540,
             ...) 
}
\arguments{
  \item{x}{object of class \code{\link{DataTreeSet}}.}
  \item{unitID}{unit ID of probeset with type of ID determined by parameter \code{unittype}.}
  \item{unittype}{\code{character} vector, one of \dQuote{unit}, \dQuote{transcript}, \dQuote{probeset}.}
  \item{which}{type of probes to be used, for details see \code{\link{validData}}.}
  \item{transfo}{a valid function to transform the data, usually \dQuote{log2}, or \dQuote{0}.}
  \item{names}{optional vector of sample names.}
  \item{ylim}{range for the plotted y values.}
  \item{col}{color to use for the different samples.}
  \item{lty}{line types to use for the different samples.}
  \item{add.legend}{logical, if TRUE a legend of sample names will be drawn. Optionally, a character
   indicating the position of the legend, default is \dQuote{topleft}.}
  \item{dev}{graphics device to plot to, i.e. one of \dQuote{screen}, \dQuote{jpeg},\dQuote{png},
   \dQuote{pdf} or \dQuote{ps}.}
  \item{outfile}{the name of the output file.}
  \item{w}{the width of the device in pixels.}
  \item{h}{the height of the device in pixels.}
  \item{\dots}{optional arguments to be passed to \code{plotProbeset}.}
}
\details{
 Produces line plots of the probe intensities for probeset \code{unitID}. Probe intensities are 
 taken from slot \code{data}. 

 For \code{names=NULL} full column names of slot \code{data} will be displayed while for 
 \code{names="namepart"} column names will be displayed without name extension. If \code{names}
 is a vector of column names, line plots of probe intensities will only be drawn for these columns.
}
\author{Christian Stratowa}
\note{
 Data must first be attached to class \code{\link{DataTreeSet}} using method \code{\link{attachInten}}.
 Furthermore, unit names must be attached using method \code{\link{attachUnitNames}}.
}
\seealso{\code{\link{probesetplot}}}
\keyword{device}
