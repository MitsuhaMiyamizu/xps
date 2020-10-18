\name{plotNUSE}
\alias{plotNUSE}
\title{Box Plots of Normalized Unscaled Standard Errors (NUSE) for Device}
\description{
 Produce boxplot of Normalized Unscaled Standard Errors (NUSE) for the set of arrays and
 the selected device.
}
\usage{
plotNUSE(x,
         which    = "UnitName",
         size     = 0,
         range    = 0,
         names    = "namepart",
         main     = "NUSE Plot",
         ylim     = c(0.8,1.2),
         las      = 2,
         add.line = TRUE,
         outline  = FALSE,
         dev      = "screen",
         outfile  = "NUSEPlot",
         w        = 800,
         h        = 540,
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
  \item{\dots}{optional arguments to be passed to \code{boxplot}.}
}
\details{
 Create boxplots of Normalized Unscaled Standard Errors (NUSE) for the set of arrays.

 For \code{names=NULL} full column names of slot \code{data} will be displayed while for 
 \code{names="namepart"} column names will be displayed without name extension. If \code{names}
 is a vector of column names, only these columns will displayed as boxplot.
}
\author{Christian Stratowa}
\seealso{\code{\link{nuseplot}}}
\keyword{device}
