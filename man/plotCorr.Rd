\name{plotCorr}
\alias{plotCorr}
\title{Array-Array Expression Level Correlation Plot for Device}
\description{
 A heat map of the array-array Spearman rank correlation coefficients for the selected device.
}
\usage{
plotCorr(x,
         which      = "UnitName",
         transfo    = log2,
         method     = "spearman",
         col        = NULL,
         names      = "namepart",
         sort       = FALSE,
         reverse    = TRUE,
         bmar       = NULL,
         add.legend = FALSE,
         dev        = "screen",
         outfile    = "CorrelationPlot",
         w          = 540,
         h          = 540,
         ...) 
}
\arguments{
  \item{x}{object of class \code{\link{ExprTreeSet}}.}
  \item{which}{type of probes to be used, for details see \code{\link{validData}}.}
  \item{transfo}{a valid function to transform the data, usually \dQuote{log2}, or \dQuote{0}.}
  \item{method}{a character string indicating which correlation coefficient is to be computed.}
  \item{col}{vector of colors for plot, length is number of samples.}
  \item{names}{optional vector of sample names.}
  \item{sort}{logical, if TRUE the correlation matrix will be sorted decreasingly.}
  \item{reverse}{logical, if TRUE the correlation matrix will be replaced by \code{1 - cor()}.}
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
 Produces a heat map of the array-array Spearman rank correlation coefficients for slot \code{data}
 for an object of class \code{\link{ExprTreeSet}}.

 For \code{names=NULL} full column names of slot \code{data} will be displayed while for 
 \code{names="namepart"} column names will be displayed without name extension. If \code{names}
 is a vector of column names, only these columns will displayed as \code{corplot}.

 For \code{bmar=NULL} the default list \code{bmar = list(b=6, cex.axis=1.0)} will be used initially.
 However, both bottom margin and axis label magnification will be adjusted depending on the number
 of label characters and the number of samples.  
}
\author{Christian Stratowa}
\seealso{\code{\link{corplot}}}
\keyword{device}
