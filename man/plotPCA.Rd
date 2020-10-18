\name{plotPCA}
\alias{plotPCA}
\title{PCA Plot for Device}
\description{
 This function produces a PCA plot of the first two principle components for the selected device.
}
\usage{
plotPCA(x,
        which      = "UnitName",
        transfo    = log2,
        method     = "none",
        groups     = NULL,
        screeplot  = FALSE,
        squarepca  = FALSE,
        pcs        = c(1,2),
        add.labels = FALSE,
        add.legend = FALSE,
        col        = NULL,
        names      = "namepart",
        as.list    = FALSE,
        dev        = "screen",
        outfile    = "PCAPlot",
        w          = 540,
        h          = 540,
        ...) 
}
\arguments{
  \item{x}{object of class \code{\link{ExprTreeSet}}.}
  \item{which}{type of probes to be used, for details see \code{\link{validData}}.}
  \item{transfo}{a valid function to transform the data, usually \dQuote{log2}, or \dQuote{0}.}
  \item{method}{a character string indicating which correlation coefficient is to be computed.
   One of \dQuote{pearson}, \dQuote{spearman}, \dQuote{kendall}, or \dQuote{none}.}
  \item{groups}{character \code{vector} listing the group names in order of the \code{names}.}
  \item{screeplot}{logical, if TRUE plot a \code{\link[stats]{screeplot}} instead of a PCA plot.}
  \item{squarepca}{logical, if TRUE make the y-axis of the PCA plot comparable to the x-axis.}
  \item{pcs}{a character \code{vector} of length two indicating which principal components to plot.}
  \item{add.labels}{logical, if TRUE then name labels will be added to the points.}
  \item{add.legend}{logical, if TRUE and \code{groups} are supplied then a legend indicating the groups
   will be drawn. Optionally, a character indicating the position of the legend, default is \dQuote{topleft}.}
  \item{col}{vector of colors for plot, length is number of samples.}
  \item{names}{optional vector of sample names.}
  \item{as.list}{logical, if TRUE then a \code{list} will be returned in addition to the plot.}
  \item{dev}{graphics device to plot to, i.e. one of \dQuote{screen}, \dQuote{jpeg},\dQuote{png},
   \dQuote{pdf} or \dQuote{ps}.}
  \item{outfile}{the name of the output file.}
  \item{w}{the width of the device in pixels.}
  \item{h}{the height of the device in pixels.}
  \item{\dots}{optional arguments to be passed to \code{plot}.}
}
\details{
 Function \code{plotPCA} produces a PCA plot of the first two principle components for slot
 \code{data} or the correlations between the columns of slot \code{data}, respectively, of an
 object of class \code{\link{ExprTreeSet}}.

 For \code{method="none"} function \code{[stats]prcomp} will be applied to slot \code{data} directly,
 otherwise \code{prcomp} will be applied to \code{(1 - cor(data))} with the respective \code{method}.

 For \code{screeplot=TRUE} a \code{\link[stats]{screeplot}} will be plotted instead of a PCA plot.

 For \code{names=NULL} full column names of slot \code{data} will be displayed while for 
 \code{names="namepart"} column names will be displayed without name extension. If \code{names}
 is a vector of column names, only these columns will displayed as mvaplot.
}
\author{Christian Stratowa}
\seealso{\code{\link{pcaplot}}}
\keyword{device}
