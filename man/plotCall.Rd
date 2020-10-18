\name{plotCall}
\alias{plotCall}
\title{Barplot of Percent Present and Absent Calls for Device}
\description{
 Creates a barplot of percent Present/Marginal/Absent calls for the selected device.
}
\usage{
plotCall(x,
         beside  = TRUE,
         names   = "namepart",
         col     = c("red","green","blue"),
         legend  = c("P","M","A"),
         ylim    = c(0,100),
         ylab    = "detection call [\%]",
         las     = 2,
         dev     = "screen",
         outfile = "CallPlot",
         w       = 800,
         h       = 540,
         ...) 
}
\arguments{
  \item{x}{object of class \code{\link{CallTreeSet}}.}
  \item{beside}{logical. If \code{FALSE}, the columns of height are portrayed as stacked bars,
   and if \code{TRUE} the columns are portrayed as juxtaposed bars.}
  \item{names}{optional vector of sample names.}
  \item{col}{color for P/M/A bars}
  \item{legend}{legend for the plot, defaults to P/M/A.}
  \item{ylim}{the y limits of the plot.}
  \item{ylab}{a label for the y axis.}
  \item{las}{the style of axis labels.}
  \item{dev}{graphics device to plot to, i.e. one of \dQuote{screen}, \dQuote{jpeg},\dQuote{png},
   \dQuote{pdf} or \dQuote{ps}.}
  \item{outfile}{the name of the output file.}
  \item{w}{the width of the device in pixels.}
  \item{h}{the height of the device in pixels.}
  \item{\dots}{optional arguments to be passed to \code{barplot}.}
}
\details{
 Creates a barplot of percent Present/Marginal/Absent calls.

 For \code{names=NULL} full column names of slot \code{data} will be displayed while for 
 \code{names="namepart"} column names will be displayed without name extension. If \code{names}
 is a vector of column names, only these columns will displayed as callplot.
}
\author{Christian Stratowa}
\seealso{\code{\link{callplot}}}
\keyword{device}
