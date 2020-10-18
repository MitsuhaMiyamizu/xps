\name{plotPM}
\alias{plotPM}
\title{Barplot of PM and MM Intensities for Device}
\description{
 Creates a barplot of mean perfect match and mismatch intensities for the selected device.
}
\usage{
plotPM(x,
       which   = "",
       size    = 0,
       transfo = NULL,
       method  = mean,
       names   = "namepart",
       beside  = TRUE,
       col     = c("red", "blue"),
       legend  = c("PM","MM"),
       las     = 2,
       ylab    = "mean intensities",
       dev     = "screen",
       outfile = "PMPlot",
       w       = 540,
       h       = 540,
       ...) 
}
\arguments{
  \item{x}{object of class \code{\link{DataTreeSet}}.}
  \item{which}{type of probes to be used, for details see \code{\link{validData}}.}
  \item{size}{length of sequence to be generated as subset.}
  \item{transfo}{a valid function to transform the data, usually \dQuote{log2}, or \dQuote{0}.}
  \item{method}{method to compute average intensities, \dQuote{mean} or \dQuote{median}.}
  \item{names}{optional vector of sample names.}
  \item{beside}{logical. If \code{FALSE}, mean intensities are portrayed as stacked bars, 
   and if \code{TRUE} the columns are portrayed as juxtaposed bars.}
  \item{col}{color of PM, MM bars.}
  \item{legend}{a vector of text used to construct a legend for the plot, 
   or a logical indicating whether a legend should be included.}
  \item{las}{the style of axis labels.}
  \item{ylab}{a title for the y axis.}
  \item{dev}{graphics device to plot to, i.e. one of \dQuote{screen}, \dQuote{jpeg},\dQuote{png},
   \dQuote{pdf} or \dQuote{ps}.}
  \item{outfile}{the name of the output file.}
  \item{w}{the width of the device in pixels.}
  \item{h}{the height of the device in pixels.}
  \item{\dots}{optional arguments to be passed to \code{barplot}.}
}
\details{
 Produces barplots of mean perfect match and mismatch intensities for slot \code{data}
 for an object of class \code{\link{DataTreeSet}}.

 For \code{names=NULL} full column names of slot \code{data} will be displayed while for 
 \code{names="namepart"} column names will be displayed without name extension. If \code{names}
 is a vector of column names, only these columns will displayed as pmplot.
}
\author{Christian Stratowa}
\seealso{\code{\link{pmplot}}}
\keyword{device}
