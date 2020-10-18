\name{plotDensity}
\alias{plotDensity}
\title{Plot Density Estimate for Device}
\description{
 Plot the density estimates for each sample for the selected device.
}
\usage{
plotDensity(x,
            which      = "",
            size       = 0,
            transfo    = log2,
            ylab       = "density",
            xlab       = "log intensity",
            names      = "namepart",
            type       = "l",
            col        = 1:6,
            lty        = 1:5,
            add.legend = FALSE,
            dev        = "screen",
            outfile    = "DensityPlot",
            w          = 540,
            h          = 540,
            verbose    = TRUE,
            ...) 
}
\arguments{
  \item{x}{object of class \code{\link{DataTreeSet}} or \code{\link{ExprTreeSet}}.}
  \item{which}{type of probes to be used, for details see \code{\link{validData}}.}
  \item{size}{length of sequence to be generated as subset.}
  \item{transfo}{a valid function to transform the data, usually \dQuote{log2}, or \dQuote{0}.}
  \item{xlab}{a title for the x axis.}
  \item{ylab}{a title for the y axis.}
  \item{names}{optional vector of sample names.}
  \item{type}{type for the plot.}
  \item{col}{colors to use for the different arrays.}
  \item{lty}{line types to use for the different arrays.}
  \item{add.legend}{logical, if TRUE then a legend will be drawn.}
  \item{dev}{graphics device to plot to, i.e. one of \dQuote{screen}, \dQuote{jpeg},\dQuote{png},
   \dQuote{pdf} or \dQuote{ps}.}
  \item{outfile}{the name of the output file.}
  \item{w}{the width of the device in pixels.}
  \item{h}{the height of the device in pixels.}
  \item{verbose}{logical, if \code{TRUE} print status information.}
  \item{\dots}{optional arguments to be passed to \code{plot}.}
}
\details{
 Plots the non-parametric density estimates for each sample.

 For \code{names=NULL} full column names of slot \code{data} will be displayed while for 
 \code{names="namepart"} column names will be displayed without name extension. If \code{names}
 is a vector of column names, only these columns will displayed as callplot.
}
\author{Christian Stratowa}
\seealso{\code{\link{hist}}}
\keyword{device}
