\name{plotBorder}
\alias{plotBorder}
\title{Plots of Border Elements for Device}
\description{
 Produce box-and-whisker plot(s) of the positive and negative feature intensities for the
 selected device.
}
\usage{
plotBorder(x,
           type    = c("pos", "neg"),
           qualopt = "raw",
           transfo = log2,
           range   = 0,
           names   = "namepart",
           ylim    = NULL,
           bmar    = NULL,
           las     = 2,
           dev     = "screen",
           outfile = "BorderPlot",
           w       = 800,
           h       = 540,
           ...) 
}
\arguments{
  \item{x}{object of class \code{\link{QualTreeSet}}.}
  \item{type}{type of border elements to be used, one of \dQuote{pos}, \dQuote{neg}, or both.}
  \item{qualopt}{character string specifying whether to draw boxplots for \dQuote{raw},
  \dQuote{adjusted}, or \dQuote{normalized} border intensities.}
  \item{transfo}{a valid function to transform the data, usually \dQuote{log2}, or \dQuote{0}.}
  \item{range}{determines how far the plot whiskers extend out from the box.}
  \item{names}{optional vector of sample names.}
  \item{ylim}{the y limits of the plot.}
  \item{bmar}{optional \code{list} for bottom margin and axis label magnification \code{cex.axis}.}
  \item{las}{the style of axis labels.}
  \item{dev}{graphics device to plot to, i.e. one of \dQuote{screen}, \dQuote{jpeg},\dQuote{png},
   \dQuote{pdf} or \dQuote{ps}.}
  \item{outfile}{the name of the output file.}
  \item{w}{the width of the device in pixels.}
  \item{h}{the height of the device in pixels.}
  \item{\dots}{optional arguments to be passed to \code{borderplot}.}
}
\details{
 Creates a boxplot of the positive and negative feature intensities for an object of class
 \code{\link{QualTreeSet}}.

 For \code{names=NULL} full tree names will be displayed while for  \code{names="namepart"} tree
 names will be displayed without name extension. If \code{names} is a vector of tree names, only
 these columns will displayed as boxplot.

 For \code{bmar=NULL} the default list \code{bmar = list(b=6, cex.axis=1.0)} will be used initially.
 However, both bottom margin \code{b} and axis label magnification \code{cex.axis} will be adjusted
 depending on the number of label characters and the number of samples.  
}
\author{Christian Stratowa}
\seealso{\code{\link{borderplot}}}
\keyword{device}
