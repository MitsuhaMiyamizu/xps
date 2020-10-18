\name{plotRLE}
\alias{plotRLE}
\title{Box Plots of Relative Log Expression (RLE) for Device}
\description{
 Produce boxplot of Relative Log Expression (RLE) for the set of arrays and
 the selected device.
}
\usage{
plotRLE(x,
        which    = "UnitName",
        size     = 0,
        range    = 0,
        names    = "namepart",
        main     = "RLE Plot",
        ylim     = c(-1.0, 1.0),
        las      = 2,
        add.line = TRUE,
        outline  = FALSE,
        dev      = "screen",
        outfile  = "RLEPlot",
        w        = 800,
        h        = 540,
        verbose  = TRUE,
        ...) 
}
\arguments{
  \item{x}{object of class \code{\link{ExprTreeSet}} or \code{\link{QualTreeSet}}.}
  \item{which}{type of probes to be used, for details see \code{\link{validData}}.}
  \item{size}{length of sequence to be generated as subset.}
  \item{range}{determines how far the plot whiskers extend out from the box.}
  \item{names}{optional vector of sample names.}
  \item{main}{the main title for the plot.}
  \item{ylim}{range for the plotted y values.}
  \item{las}{the style of axis labels.}
  \item{add.line}{\code{logical}, if TRUE a horizontal line is drawn.}
  \item{outline}{if \code{outline} is not true, the outliers are not drawn.}
  \item{dev}{graphics device to plot to, i.e. one of \dQuote{screen}, \dQuote{jpeg},\dQuote{png},
   \dQuote{pdf} or \dQuote{ps}.}
  \item{outfile}{the name of the output file.}
  \item{w}{the width of the device in pixels.}
  \item{h}{the height of the device in pixels.}
  \item{verbose}{logical, if \code{TRUE} print status information.}
  \item{\dots}{optional arguments to be passed to \code{boxplot}.}
}
\details{
 Create boxplots of Relative Log Expression (RLE) values for the set of arrays, i.e. of M plots,
 where M is determined relative to a pseudo-median reference chip.

 For \code{names=NULL} full column names of slot \code{data} will be displayed while for 
 \code{names="namepart"} column names will be displayed without name extension. If \code{names}
 is a vector of column names, only these columns will displayed as boxplot.
}
\author{Christian Stratowa}
\seealso{\code{\link{rleplot}}}
\keyword{device}
