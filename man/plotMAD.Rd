\name{plotMAD}
\alias{plotMAD}
\title{Array-Array Expression Level Distance Plot for Device}
\description{
 A false color display of between arrays distances, computed as the MAD of the M-values of
 each pair of arrays for the selected device.
}
\usage{
plotMAD(x,
        which      = "UnitName",
        transfo    = log2,
        col        = NULL,
        names      = "namepart",
        sort       = FALSE,
        bmar       = NULL,
        add.legend = FALSE,
        dev        = "screen",
        outfile    = "MADPlot",
        w          = 540,
        h          = 540,
        ...) 
}
\arguments{
  \item{x}{object of class \code{\link{ExprTreeSet}}.}
  \item{which}{type of probes to be used, for details see \code{\link{validData}}.}
  \item{transfo}{a valid function to transform the data, usually \dQuote{log2}, or \dQuote{0}.}
  \item{col}{vector of colors for plot, length is number of samples.}
  \item{names}{optional vector of sample names.}
  \item{sort}{logical, if TRUE the correlation matrix will be sorted decreasingly.}
  \item{bmar}{optional \code{list} for bottom margin and axis label magnification \code{cex.axis}.}
  \item{add.legend}{logical, if TRUE then a color bar will be drawn.}
  \item{dev}{graphics device to plot to, i.e. one of \dQuote{screen}, \dQuote{jpeg},\dQuote{png},
   \dQuote{pdf} or \dQuote{ps}.}
  \item{outfile}{the name of the output file.}
  \item{w}{the width of the device in pixels.}
  \item{h}{the height of the device in pixels.}
  \item{\dots}{optional arguments to be passed to \code{plot}.}
}
\details{
 Produces a false color display, i.e. heatmap, of between array distances for slot \code{data}
 for an object of class \code{\link{ExprTreeSet}}, computed as the MAD of the M-values of each
 pair of arrays.

 For \code{names=NULL} full column names of slot \code{data} will be displayed while for 
 \code{names="namepart"} column names will be displayed without name extension. If \code{names}
 is a vector of column names, only these columns will displayed as mdaplot.

 For \code{bmar=NULL} the default list \code{bmar = list(b=6, cex.axis=1.0)} will be used initially.
 However, both bottom margin and axis label magnification will be adjusted depending on the number
 of label characters and the number of smaples.  
}
\author{Christian Stratowa}
\seealso{\code{\link{madplot}}}
\keyword{device}
