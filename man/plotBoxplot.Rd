\name{plotBoxplot}
\alias{plotBoxplot}
\title{Box Plots for Device}
\description{
  Produce box-and-whisker plot(s) of the samples for the selected device.
}
\usage{
plotBoxplot(x,
            which   = "",
            size    = 0,
            transfo = log2,
            range   = 0,
            names   = "namepart",
            mar     = NULL,
            las     = 2,
            cex     = 1.0,
            dev     = "screen",
            outfile = "BoxPlot",
            w       = 800,
            h       = 540,
            ...) 
}
\arguments{
  \item{x}{object of class \code{\link{DataTreeSet}} or \code{\link{ExprTreeSet}}.}
  \item{which}{type of probes to be used, for details see \code{\link{validData}}.}
  \item{size}{length of sequence to be generated as subset.}
  \item{transfo}{a valid function to transform the data, usually \code{log2}, or 0.}
  \item{range}{determines how far the plot whiskers extend out from the box.}
  \item{names}{optional vector of sample names.}
  \item{mar}{plot margin.}
  \item{las}{style of axis labels.}
  \item{cex}{amount by which plotting text and symbols should be magnified.}
  \item{dev}{graphics device to plot to, i.e. one of \dQuote{screen}, \dQuote{jpeg},\dQuote{png},
   \dQuote{pdf} or \dQuote{ps}.}
  \item{outfile}{the name of the output file.}
  \item{w}{the width of the device in pixels.}
  \item{h}{the height of the device in pixels.}
  \item{\dots}{optional arguments to be passed to \code{boxplot}.}
}
\details{
 Produces a boxplot for slot \code{data} for an object of class \code{\link{DataTreeSet}},
 \code{\link{ExprTreeSet}} or \code{\link{QualTreeSet}} for the selected graphics device.
}
\author{Christian Stratowa}
\seealso{\code{\link{boxplot}}, \code{\link{plotBorder}}, \code{\link{plotNUSE}}, \code{\link{plotRLE}}}
\keyword{device}
